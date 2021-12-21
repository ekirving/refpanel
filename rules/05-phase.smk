#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import expand, touch, temp, unpack

from scripts.utils import list_source_samples, list_families, list_family_children

global workflow

"""
Rules to perform variant phasing on a joint-called reference panel  

see https://github.com/odelaneau/shapeit4/issues/17 
"""

# TODO need a version that does not rely on trios
# TODO build a `scaffold` using the trios from the extended 1000G NYGC dataset
# TODO fetch the 10x Genomics callset for the HGDP project and use that also
# TODO add benchmarks to all rules
# TODO implement the "best practices" criteria https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html


rule bcftools_subset_sample:
    """
    Subset a single sample from the join-callset, so we can efficiently parallelize the read-based phasing
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel.vcf.gz.tbi",
    output:
        vcf=temp("data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_subset.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_subset.vcf.gz.tbi"),
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools view --samples '{wildcards.sample}' -Oz -o {output.vcf} {input.vcf} && "
        "bcftools index --tbi {output.vcf}"


rule whatshap_read_based_phasing:
    """
    Annotate the VCF with read-based phase set blocks (PS), for use by shapeit4

    https://whatshap.readthedocs.io/en/latest/guide.html
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_subset.vcf.gz",
        tbi="data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_subset.vcf.gz.tbi",
        cram="data/source/{source}/cram/{sample}.cram",
        crai="data/source/{source}/cram/{sample}.cram.crai",
    output:
        vcf=temp("data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap.vcf.log",
    conda:
        "../envs/whatshap-1.2.1.yaml"
    shell:
        "whatshap phase"
        " --reference {input.ref}"
        " --tag PS"
        " --indels"
        " --ignore-read-groups"
        " --output {output.vcf}"
        " {input.vcf}"
        " {input.cram} 2> {log} && "
        "bcftools index --tbi {output.vcf}"


def bcftools_merge_phased_samples_input(wildcards):
    """
    Return a list of VCF/TBI files for each sample in the reference panel
    """
    panel = wildcards.panel
    chr = wildcards.chr

    vcf = []
    tbi = []

    for source, sample in list_source_samples(config, panel):
        vcf.append(f"data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap.vcf.gz")
        tbi.append(f"data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap.vcf.gz.tbi")

    file_list = f"data/panel/{panel}/vcf/sample/{panel}_{chr}.list"

    with open(file_list, "w") as fout:
        fout.write("\n".join(vcf) + "\n")

    return {"vcfs": vcf, "tbi": tbi, "list": file_list}


# noinspection PyUnresolvedReferences
rule bcftools_merge_phased_samples:
    """
    Merge the sample-level VCFs, with read-based phase set blocks, back into a single chromosome.

    Raise the `ulimit` to prevent exceeding the maximum number of open files. 
    """
    input:
        unpack(bcftools_merge_phased_samples_input),
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz.tbi",
    params:
        limit=lambda wildcards, input: len(input.vcf) + 10,
    conda:
        "../envs/htslib-1.14.yaml"
    resources:
        mem_mb=8 * 1024,
    shell:
        "ulimit -n {params.limit} && "
        "bcftools merge --file-list {input.list} -Oz -o {output.vcf}"


rule pedigree_family:
    """
    Extract a specific family from the pedigree.

    There are 602 trios in 1000G, but only 576 families (including quads and multi-generational families)  
    """
    input:
        ped=lambda wildcards: config["panel"][wildcards.panel]["pedigree"],
    output:
        ped="data/panel/{panel}/vcf/family/{panel}_{family}.ped",
    shell:
        """awk '$1=="{wildcards.family}"' {input.ped} > {output.ped}"""


rule bcftools_subset_family:
    """
    Subset a specific family from the join-callset, so we can efficiently parallelize the pedigree-based phasing
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel.vcf.gz.tbi",
        ped="data/panel/{panel}/vcf/family/{panel}_{family}.ped",
    output:
        vcf=temp("data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_subset.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_subset.vcf.gz.tbi"),
    params:
        samples=lambda wildcards, input: ",".join(
            set(pd.read_table(input.ped, header=None, sep=" ", usecols=[1, 2, 3]).values.flatten())
        ),
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools view --samples '{params.samples}' -Oz -o {output.vcf} {input.vcf} && "
        "bcftools index --tbi {output.vcf}"


rule whatshap_pedigree_phasing:
    """
    Build a pedigree based scaffold, for use by `shapeit4`.

    This takes the whole VCF as input, extracts a single family, and performs trio phasing.

    Raises a "Zero genetic distances encountered" warning due to high density of the genetic map.

    https://whatshap.readthedocs.io/en/latest/guide.html#phasing-pedigrees
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        map="data/reference/GRCh38/genetic_maps/whatshap/genetic_map_hg38_{chr}.map",
        ped="data/panel/{panel}/vcf/family/{panel}_{family}.ped",
        vcf="data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_subset.vcf.gz",
        tbi="data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_subset.vcf.gz.tbi",
    output:
        vcf=temp("data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family.vcf.log",
    conda:
        "../envs/whatshap-1.2.1.yaml"
    resources:
        mem_mb=6 * 1024,
    shell:
        "whatshap phase"
        " --reference {input.ref}"
        " --chromosome {wildcards.chr}"
        " --genmap {input.map}"
        " --ped {input.ped}"
        " --use-ped-samples"
        " --tag PS"
        " --indels"
        " --output {output.vcf}"
        " {input.vcf} 2> {log} && "
        "bcftools index --tbi {output.vcf}"


rule bcftools_extract_children:
    """
    Extract the children from the phased trios, as these will be merged to form the scaffold. 
    """
    input:
        vcf="data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family.vcf.gz",
        tbi="data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family.vcf.gz.tbi",
    output:
        vcf=temp("data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family_children.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family_children.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family_children.vcf.log",
    params:
        children=lambda wildcards: ",".join(list_family_children(config, wildcards.panel, wildcards.family)),
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools view --samples '{params.children}' -Oz -o {output.vcf} {input.vcf} && "
        "bcftools index --tbi {output.vcf}"


def bcftools_merge_phased_children_input(wildcards):
    """
    Return a list of VCF/TBI files for each family in the pedigree
    """
    panel = wildcards.panel
    chr = wildcards.chr

    vcf = []
    tbi = []

    for family in list_families(config, panel):
        vcf.append(f"data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family_children.vcf.gz"),
        tbi.append(f"data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family_children.vcf.gz.tbi"),

    return {"vcfs": vcf, "tbi": tbi}


# noinspection PyUnresolvedReferences
rule bcftools_merge_phased_children:
    """
    Merge the family-level trio phased VCFs back into a single chromosome.

    Produces a scaffold for `shapeit4` containing the 602 phased trios from 1000G.
    """
    input:
        unpack(bcftools_merge_phased_children_input),
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_trios.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_trios.vcf.gz.tbi",
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools merge -Oz -o {output.vcf} {input.vcfs}"


rule shapeit4_phase_vcf:
    """
    Phase the joint-callset, using the pedigree phased VCF as a scaffold and the read-based phase sets as input.

    See https://github.com/odelaneau/shapeit4/issues/17
    """
    input:
        map="data/reference/GRCh38/genetic_maps/shapeit4/{chr}.b38.gmap.gz",
        vcf1="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz",
        tbi1="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz.tbi",
        vcf2="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_trios.vcf.gz",
        tbi2="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_trios.vcf.gz.tbi",
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_phased.vcf.gz",
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_phased.vcf.log",
    threads: max(workflow.cores / 4, 8)
    conda:
        "../envs/shapeit-4.2.2.yaml"
    shell:
        "shapeit4"
        " --thread {threads}"
        " --input {input.vcf1}"
        " --scaffold {input.vcf2}"
        " --map {input.map}"
        " --region {wildcards.chr}"
        " --sequencing"
        " --output {output.vcf}"
        " --log {log}"


rule shapeit4_phase_all_vcfs:
    input:
        expand(
            "data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_phased.vcf.gz",
            chr=config["chroms"],
            allow_missing=True,
        ),
    output:
        touch("data/panel/{panel}/vcf/phase.done"),

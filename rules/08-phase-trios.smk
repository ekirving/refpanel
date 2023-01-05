#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import os, re

import pandas as pd
from snakemake.io import temp, unpack, expand, touch

from scripts.common import list_families

global workflow

"""
Rules to perform pedigree phasing of the joint-called reference panel

https://whatshap.readthedocs.io/en/latest/guide.html#phasing-pedigrees
"""


rule pedigree_family:
    """
    Extract a specific family from the pedigree.
    """
    input:
        ped=lambda wildcards: config["panel"][wildcards.panel]["pedigree"],
    output:
        ped="data/panel/{panel}/vcf/family/{panel}_{family}.ped",
    benchmark:
        "benchmarks/pedigree_family-{panel}-{family}.tsv"
    shell:
        """awk '$1=="{wildcards.family}"' {input.ped} > {output.ped}"""


rule bcftools_subset_family:
    """
    Subset a specific family from the joint-callset, so we can efficiently parallelize the pedigree-based phasing
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz.tbi",
        ped="data/panel/{panel}/vcf/family/{panel}_{family}.ped",
    output:
        vcf=temp("data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_subset.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_subset.vcf.gz.tbi"),
    params:
        samples=lambda wildcards, input: ",".join(
            sorted(set(pd.read_table(input.ped, header=None, sep=" ", usecols=[1, 2, 3]).values.flatten()))
        ),
    benchmark:
        "benchmarks/bcftools_subset_family-{panel}-{chr}-{family}.tsv"
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools view --samples '{params.samples}' -Oz -o {output.vcf} {input.vcf} && "
        "bcftools index --tbi {output.vcf}"


rule bcftools_subset_family_chrX:
    """
    Subset chrX for a family from the joint-callset, and split the PAR and non-PAR regions into separate files
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_chrX_vqsr_norm_annot_filter.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_chrX_vqsr_norm_annot_filter.vcf.gz.tbi",
        ped="data/panel/{panel}/vcf/family/{panel}_{family}.ped",
        bed1="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.M.1.bed",
        bed2="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.M.2.bed",
    output:
        vcf1=temp("data/panel/{panel}/vcf/family/{panel}_chrXm1_{family}_subset.vcf.gz"),
        tbi1=temp("data/panel/{panel}/vcf/family/{panel}_chrXm1_{family}_subset.vcf.gz.tbi"),
        vcf2=temp("data/panel/{panel}/vcf/family/{panel}_chrXm2_{family}_subset.vcf.gz"),
        tbi2=temp("data/panel/{panel}/vcf/family/{panel}_chrXm2_{family}_subset.vcf.gz.tbi"),
    params:
        samples=lambda wildcards, input: ",".join(
            sorted(set(pd.read_table(input.ped, header=None, sep=" ", usecols=[1, 2, 3]).values.flatten()))
        ),
    benchmark:
        "benchmarks/bcftools_subset_family_chrX-{panel}-{family}.tsv"
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools view --samples '{params.samples}' --regions-file {input.bed1} -Oz -o {output.vcf1} {input.vcf} && bcftools index --tbi {output.vcf1} && "
        "bcftools view --samples '{params.samples}' --regions-file {input.bed2} -Oz -o {output.vcf2} {input.vcf} && bcftools index --tbi {output.vcf2}"


rule whatshap_pedigree_phasing:
    """
    Build a pedigree based scaffold, for use by `shapeit4`.

    This takes a family-level VCF as input, and performs trio phasing of the children.

    Raises a "Zero genetic distances encountered" warning due to lack of recombination between some markers in the 
    genetic map.

    https://whatshap.readthedocs.io/en/latest/guide.html#phasing-pedigrees
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        map=lambda wildcards: "data/reference/GRCh38/genetic_maps/whatshap/genetic_map_hg38_{chr}.map".format(
            chr=re.sub("m[1-2]$", "", wildcards.chr)
        ),
        ped="data/panel/{panel}/vcf/family/{panel}_{family}.ped",
        vcf="data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_subset.vcf.gz",
        tbi="data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_subset.vcf.gz.tbi",
    output:
        vcf=temp("data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family.vcf.log",
    benchmark:
        "benchmarks/whatshap_pedigree_phasing-{panel}-{chr}-{family}.tsv"
    params:
        chr=lambda wildcards: re.sub("m[1-2]$", "", wildcards.chr),
    wildcard_constraints:
        # permit `chrXm1` and `chrXm2`, so we can handle the PAR regions of male X chromosomes
        chr="chr(\d+|Xm[1-2])",
    conda:
        "../envs/whatshap-1.2.1.yaml"
    shell:
        "whatshap phase"
        " --reference {input.ref}"
        " --chromosome {params.chr}"
        " --genmap {input.map}"
        " --ped {input.ped}"
        " --use-ped-samples"
        " --tag PS"
        " --indels"
        " --output {output.vcf}"
        " {input.vcf} 2> {log} && "
        "bcftools index --tbi {output.vcf}"


rule bcftools_concat_family_chrX:
    """
    Concatenate the PAR and non-PAR regions of chrX back together

    At present, `whatshap` cannot handle variable ploidy in a contig, so the PAR regions in chrX have to be subset
    (see https://github.com/whatshap/whatshap/issues/424)
    """
    input:
        vcf1="data/panel/{panel}/vcf/family/{panel}_chrXm1_{family}_family.vcf.gz",
        tbi1="data/panel/{panel}/vcf/family/{panel}_chrXm1_{family}_family.vcf.gz.tbi",
        vcf2="data/panel/{panel}/vcf/family/{panel}_chrXm2_{family}_family.vcf.gz",
        tbi2="data/panel/{panel}/vcf/family/{panel}_chrXm2_{family}_family.vcf.gz.tbi",
    output:
        vcf=temp("data/panel/{panel}/vcf/family/{panel}_chrX_{family}_family.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/family/{panel}_chrX_{family}_family.vcf.gz.tbi"),
    benchmark:
        "benchmarks/bcftools_concat_family_chrX-{panel}-{family}.tsv"
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools concat --allow-overlaps  -Oz -o {output.vcf} {input.vcf1} {input.vcf2} && "
        "bcftools index --tbi {output.vcf}"


# use the special case rule for trio phasing of chrX
ruleorder: bcftools_concat_family_chrX > whatshap_pedigree_phasing


def bcftools_merge_phased_families_input(wildcards):
    """
    Return a list of VCF/TBI files for each family in the pedigree
    """
    panel = wildcards.panel
    chr = wildcards.chr

    vcf = []
    tbi = []

    for family in list_families(config, panel):
        vcf.append(f"data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family.vcf.gz"),
        tbi.append(f"data/panel/{panel}/vcf/family/{panel}_{chr}_{family}_family.vcf.gz.tbi"),

    file_list = f"data/panel/{panel}/vcf/family/{panel}_{chr}.list"

    os.makedirs(os.path.dirname(file_list), exist_ok=True)

    with open(file_list, "w") as fout:
        fout.write("\n".join(vcf) + "\n")

    return {"vcfs": vcf, "tbi": tbi, "list": file_list}


# noinspection PyUnresolvedReferences
rule bcftools_merge_phased_families:
    """
    Merge the family-level trio phased VCFs back into a single chromosome.

    Produces a scaffold for `shapeit4` containing the 602 phased trios from 1000G.
    """
    input:
        unpack(bcftools_merge_phased_families_input),
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_trios.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_trios.vcf.gz.tbi",
    params:
        limit=MAX_OPEN_FILES,
    resources:
        mem_mb=8 * 1024,
    benchmark:
        "benchmarks/bcftools_merge_phased_families-{panel}-{chr}.tsv"
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "ulimit -n {params.limit} && "
        "bcftools merge --file-list {input.list} -Oz -o {output.vcf} && "
        "bcftools index --tbi {output.vcf}"


rule panel_pedigree_phasing:
    """
    Perform pedigree based phasing of all trios in a reference panel.
    """
    input:
        expand(
            "data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_trios.vcf.gz",
            chr=config["chroms"],
            allow_missing=True,
        ),
    output:
        touch("data/panel/{panel}/vcf/trios.done"),

#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import os

from snakemake.io import temp, unpack, expand, touch

from scripts.common import list_source_samples, MAX_OPEN_FILES

global workflow

"""
Rules to perform read-based phasing of the joint-called reference panel

https://whatshap.readthedocs.io/en/latest/guide.html 
"""


rule bcftools_subset_sample:
    """
    Subset a single sample from the join-callset, so we can efficiently parallelize the read-based phasing
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz.tbi",
    output:
        vcf=temp("data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_subset.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_subset.vcf.gz.tbi"),
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools view --samples '{wildcards.sample}' -Oz -o {output.vcf} {input.vcf} && "
        "bcftools index --tbi {output.vcf}"


# TODO fetch the 10x Genomics callset for the HGDP project and use that also
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

    os.makedirs(os.path.dirname(file_list), exist_ok=True)

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
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_whatshap.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_whatshap.vcf.gz.tbi",
    params:
        limit=MAX_OPEN_FILES,
    conda:
        "../envs/htslib-1.14.yaml"
    resources:
        mem_mb=8 * 1024,
    shell:
        "ulimit -n {params.limit} && "
        "bcftools merge --file-list {input.list} -Oz -o {output.vcf} && "
        "bcftools index --tbi {output.vcf}"


rule panel_read_based_phasing:
    """
    Perform read-based phasing of all samples in a reference panel
    """
    input:
        expand(
            "data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_whatshap.vcf.gz",
            chr=config["chroms"],
            allow_missing=True,
        ),
    output:
        touch("data/panel/{panel}/vcf/whatshap.done"),

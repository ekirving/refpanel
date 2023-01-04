#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import os

from snakemake.io import temp, unpack, expand, touch

from scripts.common import list_source_samples, sample_sex, MAX_OPEN_FILES

global workflow

"""
Rules to perform read-based phasing of the joint-called reference panel

https://whatshap.readthedocs.io/en/latest/guide.html 
"""


wildcard_constraints:
    # permit `chrXm1` and `chrXm2`, so we can handle the PAR regions of male X chromosomes
    chr="(chr(\d+|X(m[1-2])?|Y|M))|(others)",


rule bcftools_subset_sample:
    """
    Subset a specific sample from the joint-callset, so we can efficiently parallelize the read-based phasing
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz.tbi",
    output:
        vcf=temp("data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_subset.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_subset.vcf.gz.tbi"),
    benchmark:
        "benchmarks/bcftools_subset_sample-{panel}-{chr}-{source}-{sample}.tsv"
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools view --samples '{wildcards.sample}' -Oz -o {output.vcf} {input.vcf} && "
        "bcftools index --tbi {output.vcf}"


rule bcftools_subset_male_chrX:
    """
    Subset a male chrX sample from the joint-callset, and split the PAR and haploid regions into separate files
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_chrX_vqsr_norm_annot_filter.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_chrX_vqsr_norm_annot_filter.vcf.gz.tbi",
        bed1="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.M.1.bed",
        bed2="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.M.2.bed",
    output:
        vcf1=temp("data/panel/{panel}/vcf/sample/{panel}_chrXm1_{source}_{sample}_subset.vcf.gz"),
        tbi1=temp("data/panel/{panel}/vcf/sample/{panel}_chrXm1_{source}_{sample}_subset.vcf.gz.tbi"),
        vcf2=temp("data/panel/{panel}/vcf/sample/{panel}_chrXm2_{source}_{sample}_subset.vcf.gz"),
        tbi2=temp("data/panel/{panel}/vcf/sample/{panel}_chrXm2_{source}_{sample}_subset.vcf.gz.tbi"),
    benchmark:
        "benchmarks/bcftools_subset_male_chrX-{panel}-{source}-{sample}.tsv"
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools view --samples '{wildcards.sample}' --regions-file {input.bed1} -Oz -o {output.vcf1} {input.vcf} && bcftools index --tbi {output.vcf1} && "
        "bcftools view --samples '{wildcards.sample}' --regions-file {input.bed2} -Oz -o {output.vcf2} {input.vcf} && bcftools index --tbi {output.vcf2}"


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
    benchmark:
        "benchmarks/whatshap_read_based_phasing-{panel}-{chr}-{source}-{sample}.tsv"
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


rule whatshap_linked_read_phasing:
    """
    Incorporate phase information from linked-reads (e.g., 10x Genomics)

    At present, `whatshap` cannot directly phase from the BX tag used in 10x Genomics BAM files 
    (see https://github.com/whatshap/whatshap/pull/238).

    The suggested workflow is to use the phased VCFs output by `longranger` as input to `whatshap` via the `PHASEINPUT`
    positional argument, see:
    - https://bitbucket.org/whatshap/whatshap/issues/228/recomendation-about-10x-genomics
    - https://github.com/whatshap/whatshap/issues/243
    - https://doi.org/10.1038/s41467-017-01389-4#Sec11
    - https://doi.org/10.1038/s41467-018-08148-z
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_subset.vcf.gz",
        tbi="data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_subset.vcf.gz.tbi",
        cram="data/source/{source}/cram/{sample}.cram",
        crai="data/source/{source}/cram/{sample}.cram.crai",
        phase10x="data/source/{source}/gVCF/phase10x/{sample}.phase10x.vcf.gz",
    output:
        vcf=temp("data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap_10x.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap_10x.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap_10x.vcf.log",
    benchmark:
        "benchmarks/whatshap_linked_read_phasing-{panel}-{chr}-{source}-{sample}.tsv"
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
        " {input.cram} {input.phase10x} 2> {log} && "
        "bcftools index --tbi {output.vcf}"


rule bcftools_concat_male_chrX:
    """
    Concatenate the male PAR and haploid regions of chrX back together

    At present, `whatshap` cannot handle variable ploidy in a contig, so the PAR regions in male chrX have to be subset
    (see https://github.com/whatshap/whatshap/issues/424)
    """
    input:
        vcf1="data/panel/{panel}/vcf/sample/{panel}_chrXm1_{source}_{sample}_subset.vcf.gz",
        tbi1="data/panel/{panel}/vcf/sample/{panel}_chrXm1_{source}_{sample}_subset.vcf.gz.tbi",
        vcf2="data/panel/{panel}/vcf/sample/{panel}_chrXm2_{source}_{sample}_{whatshap}.vcf.gz",
        tbi2="data/panel/{panel}/vcf/sample/{panel}_chrXm2_{source}_{sample}_{whatshap}.vcf.gz.tbi",
    output:
        vcf=temp("data/panel/{panel}/vcf/sample/{panel}_chrXm_{source}_{sample}_{whatshap}.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/sample/{panel}_chrXm_{source}_{sample}_{whatshap}.vcf.gz.tbi"),
    benchmark:
        "benchmarks/bcftools_concat_male_chrX-{panel}-{source}-{sample}-{whatshap}.tsv"
    wildcard_constraints:
        whatshap="whatshap|whatshap_10x",
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools concat --allow-overlaps  -Oz -o {output.vcf} {input.vcf1} {input.vcf2} && "
        "bcftools index --tbi {output.vcf}"


def bcftools_merge_phased_samples_input(wildcards):
    """
    Return a list of VCF/TBI files for each sample in the reference panel
    """
    panel = wildcards.panel

    vcf = []
    tbi = []

    for source, sample, prephased in list_source_samples(config, panel):
        # handle the PAR region in male chrX
        sex = sample_sex(config, source, sample)
        chr = "chrXm" if wildcards.chr == "chrX" and sex == "M" else wildcards.chr

        if prephased:
            vcf.append(f"data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap_10x.vcf.gz")
            tbi.append(f"data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap_10x.vcf.gz.tbi")
        else:
            vcf.append(f"data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap.vcf.gz")
            tbi.append(f"data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap.vcf.gz.tbi")

    file_list = f"data/panel/{panel}/vcf/sample/{panel}_{wildcards.chr}.list"

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
    resources:
        mem_mb=8 * 1024,
    benchmark:
        "benchmarks/bcftools_merge_phased_samples-{panel}-{chr}.tsv"
    conda:
        "../envs/htslib-1.14.yaml"
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

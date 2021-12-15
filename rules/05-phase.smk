#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand, touch, temp

from scripts.utils import list_source_samples

global workflow

"""
Rules to perform variant phasing on a joint-called reference panel  

see https://github.com/odelaneau/shapeit4/issues/17 
"""

# TODO need a version that does not rely on trios
# TODO build a `scaffold` using the trios from the extended 1000G NYGC dataset
# TODO fetch the 10x Genomics callset for the HGDP project and use that also


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
        "bcftools view --samples '{wildcards.sample}' -Oz -o {output.vcf} {input.vcf}"


rule whatshap_phase_set_read_based:
    """
    Annotate the VCF with read-based phase set blocks, for use by shapeit4

    https://whatshap.readthedocs.io/en/latest/guide.html
    """
    input:
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
        " --tag=PS"
        " -o {output.vcf}"
        " {input.vcf}"
        " {input.cram} 2> {log}"


def bcftools_merge_samples_input(wildcards):
    panel = wildcards.panel
    chr = wildcards.chr
    samples = list_source_samples(config, panel)
    return [
        f"data/panel/{panel}/vcf/sample/{panel}_{chr}_{source}_{sample}_whatshap.vcf.gz" for source, sample in samples
    ]


rule bcftools_merge_samples:
    """
    Merge the sample-level VCFs, with read-based phase set blocks, back into a single chromosome.
    """
    input:
        bcftools_merge_samples_input,
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz.tbi",
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "bcftools merge -Oz -o {output.vcf} {input}"


# rule whatshap_phase_set_pedigree:
#     """
#     Build a pedigree based scaffold, for use by shapeit4
#
#     https://whatshap.readthedocs.io/en/latest/guide.html
#     https://whatshap.readthedocs.io/en/latest/guide.html#using-a-phased-vcf-instead-of-a-bam-cram-file
#     """
#     input:
#         ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
#         map="data/reference/GRCh38/genetic_maps/{chr}.b38.gmap.gz",
#         ped="data/source/1000g/1000g-trios.ped",
#         vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel.vcf.gz",
#         tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel.vcf.gz.tbi",
#     output:
#         vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz",
#         tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz.tbi",
#     log:
#         log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.log",
#     conda:
#         "../envs/whatshap-1.2.1.yaml"
#     shell:
#         "whatshap phase"
#         " --reference={input.ref}"
#         " --genmap {input.map}"
#         " --ped {input.ped}"
#         " --tag=PS"
#         " -o {output.vcf}"
#         " {input.vcf}"
#         " input.bam"


rule shapeit4_phase_vcf_trios:
    """
    Phase the joint-callset, using trios
    """
    input:
        map="data/reference/GRCh38/genetic_maps/{chr}.b38.gmap.gz",
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz.tbi",
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
        " --input {input.vcf}"
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

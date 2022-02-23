#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import math

from snakemake.io import protected, unpack, temp, expand, touch

from scripts.common import list_samples, MAX_MEM_MB, GATK_BATCH_SIZE, JAVA_TEMP_DIR

"""
Rules to merge all sample-level gVCF files for a data source
"""


wildcard_constraints:
    chr="(chr(\d+|X|Y|M))|(others)",
    type="SNP|INDEL",


def gatk3_batch_sample_chrom_gvcfs_input(wildcards):
    """Split the samples into batches so we don't use too much RAM"""
    source = wildcards.source
    samples = list_samples(config, source)
    start = (int(wildcards.batch) - 1) * GATK_BATCH_SIZE
    stop = start + GATK_BATCH_SIZE

    return {
        "ref": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "chr": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.{chr}.bed",
        "gvcfs": [f"data/source/{source}/gVCF/{sample}.g.vcf.gz" for sample in samples[start:stop]],
        "tbi": [f"data/source/{source}/gVCF/{sample}.g.vcf.gz.tbi" for sample in samples[start:stop]],
    }


# noinspection PyUnresolvedReferences
rule gatk3_batch_sample_chrom_gvcfs:
    """
    Combine all gVCFs in batches, from one chromosome in one datasource into a multi-sample gVCF
    """
    input:
        unpack(gatk3_batch_sample_chrom_gvcfs_input),
    output:
        vcf=temp("data/source/{source}/gVCF/merged/{source}_{chr}_{batch}.g.vcf.gz"),
        tbi=temp("data/source/{source}/gVCF/merged/{source}_{chr}_{batch}.g.vcf.gz.tbi"),
    log:
        log="data/source/{source}/gVCF/merged/{source}_{chr}_{batch}.g.vcf.log",
    params:
        gvcfs=lambda wildcards, input: [f"--variant {gvcf}" for gvcf in input.gvcfs],
    resources:
        mem_mb=min(28 * 1024, MAX_MEM_MB),  # ~3.71%
        tmpdir=JAVA_TEMP_DIR,
    benchmark:
        "benchmarks/gatk3_batch_sample_chrom_gvcfs-{source}-{chr}-{batch}.tsv"
    conda:
        # a bug in gatk v3.5 causes excessive memory usage when combining large numbers of samples
        "../envs/gatk-3.8.yaml"
    shell:
        "gatk3"
        " -XX:ConcGCThreads=1"
        " -Xms{resources.mem_mb}m"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " -T CombineGVCFs"
        " -R {input.ref}"
        " -L {input.chr}"
        " {params.gvcfs}"
        " -o {output.vcf} &> {log}"


def gatk3_multisample_chrom_gvcf_input(wildcards):
    """Split the samples into batches"""
    source = wildcards.source
    chr = wildcards.chr
    num_samples = len(list_samples(config, source))
    batches = range(1, math.ceil(num_samples / GATK_BATCH_SIZE) + 1)

    return {
        "ref": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "chr": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.{chr}.bed",
        "gvcfs": [f"data/source/{source}/gVCF/merged/{source}_{chr}_{batch}.g.vcf.gz" for batch in batches],
        "tbi": [f"data/source/{source}/gVCF/merged/{source}_{chr}_{batch}.g.vcf.gz.tbi" for batch in batches],
    }


# noinspection PyUnresolvedReferences
rule gatk3_multisample_chrom_gvcf:
    """
    Combine all gVCFs batches into one multi-sample gVCF for each chromosome
    """
    input:
        unpack(gatk3_multisample_chrom_gvcf_input),
    output:
        vcf=protected("data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.gz"),
        tbi=protected("data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.gz.tbi"),
    log:
        log="data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.log",
    params:
        gvcfs=lambda wildcards, input: [f"--variant {gvcf}" for gvcf in input.gvcfs],
    resources:
        mem_mb=min(50 * 1024, MAX_MEM_MB),  # ~6.67%
        tmpdir=JAVA_TEMP_DIR,
    benchmark:
        "benchmarks/gatk3_multisample_chrom_gvcf-{source}-{chr}.tsv"
    conda:
        # a bug in gatk v3.5 causes excessive memory usage when combining large numbers of samples
        "../envs/gatk-3.8.yaml"
    shell:
        "gatk3"
        " -XX:ConcGCThreads=1"
        " -Xms{resources.mem_mb}m"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " -T CombineGVCFs"
        " -R {input.ref}"
        " -L {input.chr}"
        " {params.gvcfs}"
        " -o {output.vcf} &> {log}"


rule source_merge_gvcfs:
    """
    Merge all `gVCF` files in a data source.
    """
    input:
        expand("data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.gz", chr=config["chroms"], allow_missing=True),
    output:
        touch("data/source/{source}/gVCF/merge.done"),

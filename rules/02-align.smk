#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import protected, temp, touch

from scripts.common import (
    fastq_path,
    read_group,
    list_accessions,
    GATK_NUM_THREADS,
    JAVA_MEMORY_MB,
    JAVA_TEMP_DIR,
    gatk_flags,
)

global workflow

"""
Rules to perform short-read alignment for the IGSR pipeline
"""


wildcard_constraints:
    source="[^_]+",


rule fastp_trim_adapters_se:
    """
    Pre-process single-end FASTQ files. 

    Detect and trim sequencing adapters, and remove reads that are low quality, have too many Ns, or are too short. 
    """
    input:
        fastq=lambda wildcards: fastq_path(config, wildcards.source, wildcards.accession),
    output:
        fastq=temp("data/source/{source}/fastq/{accession}_trim.fastq.gz"),
        json="data/source/{source}/fastq/{accession}_se_trim.json",
        html="data/source/{source}/fastq/{accession}_se_trim.html",
    log:
        log="data/source/{source}/fastq/{accession}_se_trim.log",
    threads: 4
    conda:
        "../envs/fastp-0.23.2.yaml"
    shell:
        "fastp "
        " --in1 {input.fastq}"
        " --out1 {output.fastq}"
        " --thread {threads}"
        " --json {output.json}"
        " --html {output.html} 2> {log}"


rule fastp_trim_adapters_pe:
    """
    Pre-process paired-end FASTQ files.

    Detect and trim sequencing adapters, and remove reads that are low quality, have too many Ns, or are too short.
    """
    input:
        fastq_r1=lambda wildcards: fastq_path(config, wildcards.source, wildcards.accession, "r1"),
        fastq_r2=lambda wildcards: fastq_path(config, wildcards.source, wildcards.accession, "r2"),
    output:
        fastq_r1=temp("data/source/{source}/fastq/{accession}_trim.r1.fastq.gz"),
        fastq_r2=temp("data/source/{source}/fastq/{accession}_trim.r2.fastq.gz"),
        json="data/source/{source}/fastq/{accession}_pe_trim.json",
        html="data/source/{source}/fastq/{accession}_pe_trim.html",
    log:
        log="data/source/{source}/fastq/{accession}_pe_trim.log",
    threads: 4
    conda:
        "../envs/fastp-0.23.2.yaml"
    shell:
        "fastp "
        " --in1 {input.fastq_r1}"
        " --in2 {input.fastq_r2}"
        " --out1 {output.fastq_r1}"
        " --out2 {output.fastq_r2}"
        " --thread {threads}"
        " --json {output.json}"
        " --html {output.html} 2> {log}"


rule bwa_mem_se:
    """
    Align a single-end FASTQ file to the reference.
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        bwt="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt",
        fastq="data/source/{source}/fastq/{accession}_trim.fastq.gz",
    output:
        bam=temp("data/source/{source}/bam/{accession}.se.bam"),
    log:
        log="data/source/{source}/bam/{accession}.se.bam.log",
    params:
        rg=lambda wildcards: read_group(config, wildcards.source, wildcards.accession),
    threads: max(workflow.cores / 4, 8)
    conda:
        "../envs/bwa-0.7.15.yaml"
    shell:
        "( bwa mem -Y "
        "   -K 100000000 "
        "   -t {threads} "
        "   -R '{params.rg}' "
        "   {input.ref} "
        "   {input.fastq} | "
        "  samtools view -Shb -o {output.bam} - "
        ") 2> {log}"


rule bwa_mem_pe:
    """
    Align a paired-end FASTQ file to the reference.
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        bwt="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt",
        fastq_r1="data/source/{source}/fastq/{accession}_trim.r1.fastq.gz",
        fastq_r2="data/source/{source}/fastq/{accession}_trim.r2.fastq.gz",
    output:
        bam=temp("data/source/{source}/bam/{accession}.pe.bam"),
    log:
        log="data/source/{source}/bam/{accession}.pe.bam.log",
    params:
        rg=lambda wildcards: read_group(config, wildcards.source, wildcards.accession),
    threads: max(workflow.cores / 4, 8)
    conda:
        "../envs/bwa-0.7.15.yaml"
    shell:
        "( bwa mem -Y "
        "   -K 100000000 "
        "   -t {threads} "
        "   -R '{params.rg}' "
        "   {input.ref} "
        "   {input.fastq_r1} "
        "   {input.fastq_r1} | "
        "  samtools view -Shb -o {output.bam} - "
        ") 2> {log}"


rule picard_fix_mate_info:
    """
    Fix mate information in the BAM
    """
    input:
        bam="data/source/{source}/bam/{accession}.pe.bam",
    output:
        bam=temp("data/source/{source}/bam/{accession}_fixedmate.bam"),
    log:
        log="data/source/{source}/bam/{accession}_fixedmate.bam.log",
    resources:
        mem_mb=JAVA_MEMORY_MB,
    conda:
        "../envs/picard-2.5.0.yaml"
    shell:
        "picard"
        " -Xmx{resources.mem_mb}m"
        " FixMateInformation"
        " MAX_RECORDS_IN_RAM=2000000"
        " VALIDATION_STRINGENCY=SILENT"
        " ADD_MATE_CIGAR=True"
        " ASSUME_SORTED=true"
        " I={input.bam}"
        " O={output.bam} 2> {log}"


def picard_merge_accessions_input(wildcards):
    """
    List all the accession BAMs for the current sample
    """
    source = wildcards.source

    bam = []
    for accession, layout, unpaired in list_accessions(config, wildcards.source, wildcards.sample):
        if layout == "PAIRED":
            bam.append(f"data/source/{source}/bam/{accession}_fixedmate.bam")

        if layout == "SINGLE" or unpaired:
            bam.append(f"data/source/{source}/bam/{accession}.se.bam")

    return bam


rule picard_merge_accessions:
    """
    Merging lane level bam files to Sample level bam files
    """
    input:
        picard_merge_accessions_input,
    output:
        bam=temp("data/source/{source}/bam/{sample}_merged.bam"),
    log:
        log="data/source/{source}/bam/{sample}_merged.bam.log",
    params:
        bams=lambda wildcards, input: [f"INPUT={bam}" for bam in input],
    resources:
        mem_mb=JAVA_MEMORY_MB,
        tmpdir=JAVA_TEMP_DIR,
    conda:
        "../envs/picard-2.5.0.yaml"
    shell:
        "picard"
        " -XX:ConcGCThreads=1"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " MergeSamFiles"
        " USE_THREADING=true"
        " MAX_RECORDS_IN_RAM=2000000"
        " VALIDATION_STRINGENCY=SILENT"
        " SORT_ORDER=queryname"
        " {params.bams}"
        " OUTPUT={output.bam} 2> {log}"


rule picard_sort_bam:
    """
    Coordinate sort BAM
    """
    input:
        bam="data/source/{source}/bam/{sample}_merged.bam",
    output:
        bam=temp("data/source/{source}/bam/{sample}_merged_sorted.bam"),
        bai=temp("data/source/{source}/bam/{sample}_merged_sorted.bai"),
    log:
        log="data/source/{source}/bam/{sample}_merged_sorted.bam.log",
    resources:
        mem_mb=JAVA_MEMORY_MB,
        tmpdir=JAVA_TEMP_DIR,
    conda:
        "../envs/picard-2.5.0.yaml"
    shell:
        "picard"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " SortSam"
        " MAX_RECORDS_IN_RAM=2000000"
        " VALIDATION_STRINGENCY=SILENT"
        " SORT_ORDER=coordinate"
        " CREATE_INDEX=true"
        " I={input.bam}"
        " O={output.bam} 2> {log}"


rule picard_mark_duplicates:
    """
    Mark duplicates

    https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md#duplicate-marking
    """
    input:
        bam="data/source/{source}/bam/{sample}_merged_sorted.bam",
    output:
        bam=temp("data/source/{source}/bam/{sample}_merged_sorted_dedup.bam"),
        bai=temp("data/source/{source}/bam/{sample}_merged_sorted_dedup.bam.bai"),
        met="data/source/{source}/bam/{sample}_merged_sorted_dedup.metrics",
    log:
        log="data/source/{source}/bam/{sample}_merged_sorted_dedup.bam.log",
    resources:
        mem_mb=JAVA_MEMORY_MB,
        tmpdir=JAVA_TEMP_DIR,
    conda:
        "../envs/picard-2.5.0.yaml"
    shell:
        "picard"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " MarkDuplicates"
        " MAX_RECORDS_IN_RAM=2000000"
        " VALIDATION_STRINGENCY=SILENT"
        " M={output.met}"
        " I={input.bam}"
        " O={output.bam} 2> {log} && "
        "picard BuildBamIndex I={output.bam} O={output.bai} 2> /dev/null"


ruleorder: picard_mark_duplicates > bwa_mem_pe


rule gatk3_base_recalibrator:
    """
    Generate the base recalibration table

    https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md#base-quality-score-recalibration
    """
    input:
        bam="data/source/{source}/bam/{sample}_merged_sorted_dedup.bam",
        bai="data/source/{source}/bam/{sample}_merged_sorted_dedup.bam.bai",
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        snps="data/reference/GRCh38/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz",
        indels="data/reference/GRCh38/other_mapping_resources/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels_fixed.vcf.gz",
        mills="data/reference/GRCh38/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz",
    output:
        tbl=temp("data/source/{source}/bam/{sample}_merged_sorted_dedup_recal.table"),
    log:
        log="data/source/{source}/bam/{sample}_merged_sorted_dedup_recal.table.log",
    threads: GATK_NUM_THREADS
    resources:
        mem_mb=JAVA_MEMORY_MB,
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -Xmx{resources.mem_mb}m"
        " -T BaseRecalibrator"
        " --downsample_to_fraction 0.1"
        " --num_cpu_threads_per_data_thread {threads}"
        " --preserve_qscores_less_than 6"
        " -R {input.ref}"
        " -o {output.tbl}"
        " -I {input.bam}"
        " -knownSites {input.snps}"
        " -knownSites {input.indels}"
        " -knownSites {input.mills} 2> {log}"


rule gatk3_recalibrator_print_reads:
    """
    Recalibrate base quality scores using known SNPs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        bam="data/source/{source}/bam/{sample}_merged_sorted_dedup.bam",
        bai="data/source/{source}/bam/{sample}_merged_sorted_dedup.bam.bai",
        tbl="data/source/{source}/bam/{sample}_merged_sorted_dedup_recal.table",
    output:
        bam=temp("data/source/{source}/bam/{sample}_merged_sorted_dedup_recal.bam"),
    log:
        log="data/source/{source}/bam/{sample}_merged_sorted_dedup_recal.bam.log",
    params:
        # fetch any sample specific flags
        flags=lambda wildcards: gatk_flags(config, wildcards.source, wildcards.sample),
    threads: GATK_NUM_THREADS
    resources:
        mem_mb=JAVA_MEMORY_MB,
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -Xmx{resources.mem_mb}m"
        " -T PrintReads"
        " --num_cpu_threads_per_data_thread {threads}"
        " --disable_indel_quals"
        " --preserve_qscores_less_than 6"
        " {params.flags}"
        " -SQQ 10"
        " -SQQ 20"
        " -SQQ 30"
        " -rf BadCigar"
        " -R {input.ref}"
        " -o {output.bam}"
        " -I {input.bam}"
        " -BQSR {input.tbl} 2> {log}"


rule samtools_cram:
    """
    Creating CRAM files
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        bam="data/source/{source}/bam/{sample}_merged_sorted_dedup_recal.bam",
    output:
        cram=protected("data/source/{source}/cram/{sample}.cram"),
        crai=protected("data/source/{source}/cram/{sample}.cram.crai"),
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "samtools view"
        " --cram"
        " --reference {input.ref}"
        " --write-index"
        " --output {output.cram}"
        " {input.bam}"


def source_list_all_crams(wildcards):
    """List all CRAM files for the given data source"""
    source = wildcards.source
    samples = pd.read_table(config["source"][source]["samples"])

    return [f"data/source/{source}/cram/{sample}.cram" for sample in samples["sample"]]


rule source_align_samples:
    """
    Align all samples in a data source.
    """
    input:
        source_list_all_crams,
    output:
        touch("data/source/{source}/cram/align.done"),

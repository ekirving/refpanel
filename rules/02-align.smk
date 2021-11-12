#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand, protected, temp

global workflow

"""
Rules to perform short-read alignment for the IGSR pipeline

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/1000G_README_2019April10_NYGCjointcalls.pdf 
"""

from scripts.utils import fastq_path, read_group, list_accessions


rule bwa_mem_pe:
    """
    Alignment at accession level

    https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md#alignment
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        bwt="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt",
        fastq_r1=lambda wildcards: fastq_path(config, wildcards.source, wildcards.accession, "r1"),
        fastq_r2=lambda wildcards: fastq_path(config, wildcards.source, wildcards.accession, "r2"),
    output:
        bam=temp("data/source/{source}/bam/{accession}.bam"),
    log:
        log="data/source/{source}/bam/{accession}.log",
    params:
        rg=lambda wildcards: read_group(config, wildcards.source, wildcards.accession),
    threads: workflow.cores / 4
    shell:
        "( bwa mem -Y "
        "   -K 100000000 "
        "   -t {threads} "
        "   -R '{params.rg}' "
        "   {input.ref} "
        "   {input.fastq_r1} "
        "   {input.fastq_r1} | "
        "  samtools view -Shb -o {output.bam} -"
        ") 2> {log}"


rule picard_fix_mate_info:
    """
    Fix mate information in the BAM
    """
    input:
        bam="data/source/{source}/bam/{accession}.bam",
    output:
        bam=temp("data/source/{source}/bam/{accession}_fixedmate.bam"),
    log:
        log="data/source/{source}/bam/{accession}_fixedmate.log",
    shell:
        "picard"
        " FixMateInformation"
        " MAX_RECORDS_IN_RAM=2000000"
        " VALIDATION_STRINGENCY=SILENT"
        " ADD_MATE_CIGAR=True"
        " ASSUME_SORTED=true"
        " I={input.bam}"
        " O={output.bam} 2> {log}"


def picard_merge_bam_files_input(wildcards):
    return expand(
        "data/source/{source}/bam/{accession}_fixedmate.bam",
        source=wildcards.source,
        accession=list_accessions(config, wildcards.source, wildcards.sample),
    )


rule picard_merge_bam_files:
    """
    Merging lane level bam files to Sample level bam files
    """
    input:
        picard_merge_bam_files_input,
    output:
        bam=temp("data/source/{source}/bam/{sample}_merged.bam"),
    log:
        log="data/source/{source}/bam/{sample}_merged.log",
    params:
        bams=lambda wildcards, input: [f"INPUT={bam}" for bam in input],
    shell:
        "picard"
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
        log="data/source/{source}/bam/{sample}_merged_sorted.log",
    shell:
        "picard"
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
        log="data/source/{source}/bam/{sample}_merged_sorted_dedup.log",
    shell:
        "picard"
        " MarkDuplicates"
        " MAX_RECORDS_IN_RAM=2000000"
        " VALIDATION_STRINGENCY=SILENT"
        " M={output.met}"
        " I={input.bam}"
        " O={output.bam} 2> {log} && "
        "samtools index {output.bam}"


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
        log="data/source/{source}/bam/{sample}_merged_sorted_dedup_recal_table.log",
    shell:
        "gatk3"
        " -T BaseRecalibrator"
        " --downsample_to_fraction 0.1"
        " -nct 4"
        " --preserve_qscores_less_than 6"
        " -R {input.ref}"
        " -o {output.tbl}"
        " -I {input.bam}"
        " -knownSites {input.snps}"
        " -knownSites {input.indels}"
        " -knownSites {input.mills} 2> {log}"


rule gatk3_print_reads:
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
        log="data/source/{source}/bam/{sample}_merged_sorted_dedup_recal.log",
    shell:
        "gatk3"
        " -T PrintReads"
        " -nct 4"
        " --disable_indel_quals"
        " --preserve_qscores_less_than 6"
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
    shell:
        "samtools view -C -T {input.ref} -o {output.cram} {input.bam} && "
        "samtools index {output.cram}"
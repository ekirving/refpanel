#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import protected, unpack, temp, touch

from scripts.common import sample_sex, JAVA_MEMORY_MB

"""
Rules to perform sample-level genotype calling for the IGSR pipeline

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/1000G_README_2019April10_NYGCjointcalls.pdf 
"""


rule gatk3_haplotype_caller:
    """
    Raw variant calls using HaplotypeCaller on a single sample, within a given region and sex-based ploidy

    see https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls_updated/README_2021November05_NYGCrawcalls_updated.docx
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cram="data/source/{source}/cram/{sample}.cram",
        region="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.{sex}.{ploidy}.bed",
    output:
        vcf=temp("data/source/{source}/gVCF/{sample}.{sex}.{ploidy}.g.vcf.gz"),
        tbi=temp("data/source/{source}/gVCF/{sample}.{sex}.{ploidy}.g.vcf.gz.tbi"),
    log:
        log="data/source/{source}/gVCF/{sample}.{sex}.{ploidy}.g.vcf.log",
    resources:
        mem_mb=JAVA_MEMORY_MB,
    benchmark:
        "benchmarks/gatk3_haplotype_caller-{source}-{sample}-{sex}-{ploidy}.tsv"
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -XX:ConcGCThreads=1"
        " -Xmx{resources.mem_mb}m"
        " -T HaplotypeCaller"
        " --genotyping_mode DISCOVERY"
        " -A AlleleBalanceBySample"
        " -A DepthPerAlleleBySample"
        " -A DepthPerSampleHC"
        " -A InbreedingCoeff"
        " -A MappingQualityZeroBySample"
        " -A StrandBiasBySample"
        " -A Coverage"
        " -A FisherStrand"
        " -A HaplotypeScore"
        " -A MappingQualityRankSumTest"
        " -A MappingQualityZero"
        " -A QualByDepth"
        " -A RMSMappingQuality"
        " -A ReadPosRankSumTest"
        " -A VariantType"
        " --logging_level INFO"
        " --emitRefConfidence GVCF"
        " -rf BadCigar"
        " --variant_index_parameter 128000"
        " --variant_index_type LINEAR"
        " --sample_ploidy {wildcards.ploidy}"
        " --intervals {input.region}"
        " -R {input.ref}"
        " --num_cpu_threads_per_data_thread 1"
        " -I {input.cram}"
        " -o {output.vcf} 2> {log}"


def gatk3_combine_ploidy_regions_input(wildcards):
    """Handle sex-dependent ploidy"""
    source = wildcards.source
    sample = wildcards.sample
    sex = sample_sex(config, source, sample)

    return {
        "ref": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "vcf1": f"data/source/{source}/gVCF/{sample}.{sex}.1.g.vcf.gz",
        "tbi1": f"data/source/{source}/gVCF/{sample}.{sex}.1.g.vcf.gz.tbi",
        "vcf2": f"data/source/{source}/gVCF/{sample}.{sex}.2.g.vcf.gz",
        "tbi2": f"data/source/{source}/gVCF/{sample}.{sex}.2.g.vcf.gz.tbi",
    }


# noinspection PyUnresolvedReferences
rule gatk3_combine_ploidy_regions:
    """
    Merge the separate gVCFs called with sex-specific ploidy in chrX and chrY
    """
    input:
        unpack(gatk3_combine_ploidy_regions_input),
    output:
        vcf=protected("data/source/{source}/gVCF/{sample}.g.vcf.gz"),
        tbi=protected("data/source/{source}/gVCF/{sample}.g.vcf.gz.tbi"),
    log:
        log="data/source/{source}/gVCF/{sample}.g.vcf.log",
    resources:
        mem_mb=JAVA_MEMORY_MB,
    benchmark:
        "benchmarks/gatk3_combine_ploidy_regions-{source}-{sample}.tsv"
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -Xmx{resources.mem_mb}m"
        " -T CombineGVCFs"
        " -R {input.ref}"
        " --variant {input.vcf1}"
        " --variant {input.vcf2}"
        " -o {output.vcf} 2> {log}"


def source_list_all_gvcfs(wildcards):
    """List all gVCF files for the given data source"""
    source = wildcards.source
    samples = pd.read_table(config["source"][source]["samples"])

    return [f"data/source/{source}/gVCF/{sample}.g.vcf.gz" for sample in samples["sample"]]


rule source_call_all_gvcfs:
    """
    Call genotype likelihoods in all samples in a data source.
    """
    input:
        source_list_all_gvcfs,
    output:
        touch("data/source/{source}/gVCF/call.done"),

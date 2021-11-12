#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import protected, unpack

from scripts.utils import sample_sex

global workflow

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
        vcf=protected("data/source/{source}/gVCF/{sample}.{sex}.{ploidy}.g.vcf.gz"),
        tbi=protected("data/source/{source}/gVCF/{sample}.{sex}.{ploidy}.g.vcf.gz.tbi"),
    log:
        log="data/source/{source}/gVCF/{sample}.{sex}.{ploidy}.g.vcf.log",
    resources:
        mem_mb=8*1024
    shell:
        "gatk3"
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
        " -nct 1"
        " -I {input.cram}"
        " -o {output.vcf} 2> {log}"


def gatk3_combine_gvcfs_input(wildcards):
    """Handle sex-dependent ploidy"""
    sex = sample_sex(config, wildcards.source, wildcards.sample)

    if sex == "M":
        gvcfs = ["data/source/{source}/gVCF/{sample}.M.1.g.vcf.gz", "data/source/{source}/gVCF/{sample}.M.2.g.vcf.gz"]
    else:
        gvcfs = ["data/source/{source}/gVCF/{sample}.F.2.g.vcf.gz"]

    return {"ref": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa", "gvcfs": gvcfs}


# noinspection PyUnresolvedReferences
rule gatk3_combine_gvcfs:
    """
    Merge the separate gVCFs called with sex-specific ploidy in chrX and chrY
    """
    input:
        unpack(gatk3_combine_gvcfs_input),
    output:
        vcf=protected("data/source/{source}/gVCF/{sample}.g.vcf.gz"),
        tbi=protected("data/source/{source}/gVCF/{sample}.g.vcf.gz.tbi"),
    log:
        log="data/source/{source}/gVCF/{sample}.g.vcf.log",
    shell:
        "gatk3"
        " -T CombineGVCFs"
        " -R {input.ref}"
        " --variant {input.gvcfs}"
        " -o {output.vcf}"
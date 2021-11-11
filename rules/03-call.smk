#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import protected

global workflow

"""
Rules to perform sample-level genotype calling for the IGSR pipeline

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/1000G_README_2019April10_NYGCjointcalls.pdf 
"""


rule gatk3_haplotype_caller:
    """
    Raw variant calls using HaplotypeCaller on single sample

    https://htmlpreview.github.io/?https://raw.githubusercontent.com/broadinstitute/gatk-docs/master/gatk3-tooldocs/3.5-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.html
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cram="data/source/{source}/cram/{sample}.cram",
    output:
        vcf=protected("data/source/{source}/gVCF/{sample}.g.vcf.gz"),
        tbi=protected("data/source/{source}/gVCF/{sample}.g.vcf.gz.tbi"),
    log:
        log="data/source/{source}/vcf/{sample}.log",
    shell:
        "gatk3"
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
        " -l INFO"
        " --emitRefConfidence GVCF"
        " -rf BadCigar"
        " --variant_index_parameter 128000"
        " --variant_index_type LINEAR"
        " -R {input.ref}"
        " -nct 1"
        " -I {input.cram}"
        " -o {output.vcf}"

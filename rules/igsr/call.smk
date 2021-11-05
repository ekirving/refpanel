#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand, protected

global workflow

"""
Rules to implement the International Genome Sample Resource (IGSR) genotype calling pipeline

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/1000G_README_2019April10_NYGCjointcalls.pdf 
"""


rule gatk3_haplotype_caller:
    """
    Raw variant calls using HaplotypeCaller on single sample

    https://htmlpreview.github.io/?https://raw.githubusercontent.com/broadinstitute/gatk-docs/master/gatk3-tooldocs/3.5-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.html
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cram="data/samples/bam/{sample}.cram",
    output:
        gvcf=protected("data/samples/vcf/{sample}.g.vcf"),
    log:
        log="data/samples/vcf/{sample}.log",
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
        " -o {output.gvcf}"


rule gatk3_genotype_gvcf:
    """
    Jointly call genotypes in all samples
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        gvcfs=expand("data/samples/vcf/{sample}.g.vcf", sample=["TODO"]),
    output:
        vcf=protected("data/samples/vcf/{panel}.vcf"),
    log:
        log="data/samples/vcf/{panel}.log",
    shell:
        "gatk3"
        " -T GenotypeGVCFs"
        " -R {input.ref}"
        " -nt 5"
        " --disable_auto_index_creation_and_locking_when_reading_rods"
        " --variant {input.gvcfs}"
        " -o {output.vcf}"


rule gatk3_variant_recalibrator_snp:
    """
    Variant Quality Score Recalibration (VQSR) to assign FILTER status for SNPs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/samples/vcf/{panel}.vcf",
        hap="data/reference/GRCh38/other_mapping_resources/hapmap_3.3.hg38.vcf.gz",
        omni="data/reference/GRCh38/other_mapping_resources/1000G_omni2.5.hg38.vcf.gz",
        snps="data/reference/GRCh38/other_mapping_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="data/reference/GRCh38/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz",
    output:
        recal="data/samples/vcf/{panel}_vqsr_snp.recal",
        tranche="data/samples/vcf/{panel}_vqsr_snp.tranches",
        plot="data/samples/vcf/{panel}_vqsr_snp_plots.R",
    log:
        log="data/samples/vcf/{panel}.log",
    shell:
        "gatk3"
        " -T VariantRecalibrator"
        " -R {input.ref}"
        " -nt 5"
        " -input {input.vcf}"
        " -mode SNP"
        " -recalFile {output.recal}"
        " -tranchesFile {output.tranche}"
        " -rscriptFile {output.plot}"
        " -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hap}"
        " -resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni}"
        " -resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.snps}"
        " -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp}"
        " -an QD"
        " -an MQ"
        " -an FS"
        " -an MQRankSum"
        " -an ReadPosRankSum"
        " -an SOR"
        " -an DP"
        " -tranche 100.0"
        " -tranche 99.8"
        " -tranche 99.6"
        " -tranche 99.4"
        " -tranche 99.2"
        " -tranche 99.0"
        " -tranche 95.0"
        " -tranche 90.0"


rule gatk3_variant_recalibrator_indel:
    """
    Variant Quality Score Recalibration (VQSR) to assign FILTER status for INDELs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/samples/vcf/{panel}.vcf",
        mills="data/reference/GRCh38/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz",
        dbsnp="data/reference/GRCh38/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz",
    output:
        recal="data/samples/vcf/{panel}_vqsr_indel.recal",
        tranche="data/samples/vcf/{panel}_vqsr_indel.tranches",
        plot="data/samples/vcf/{panel}_vqsr_indel_plots.R",
    log:
        log="data/samples/vcf/{panel}.log",
    shell:
        "gatk3"
        " -T VariantRecalibrator"
        " -R {input.ref}"
        " -nt 5"
        " -input {input.vcf}"
        " -mode INDEL"
        " -recalFile {output.recal}"
        " -tranchesFile {output.tranche}"
        " -rscriptFile {output.plot}"
        " -resource:mills,known=true,training=true,truth=true,prior=12.0 {input.mills}"
        " -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp}"
        " -an QD"
        " -an FS"
        " -an ReadPosRankSum"
        " -an MQRankSum"
        " -an SOR"
        " -an DP"
        " -tranche 100.0"
        " -tranche 99.0"
        " -tranche 95.0"
        " -tranche 92.0"
        " -tranche 90.0"
        " --maxGaussians 4"


rule gatk3_apply_recalibration_snp:
    """
    Apply the VQSR to the SNPs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/samples/vcf/{panel}.vcf",
        recal="data/samples/vcf/{panel}_vqsr_snp.recal",
        tranche="data/samples/vcf/{panel}_vqsr_snp.tranches",
    output:
        vcf="data/samples/vcf/{panel}_vqsr_snp.vcf",
    log:
        log="data/samples/vcf/{panel}_vqsr_snp.log",
    shell:
        "gatk3"
        " -T ApplyRecalibration"
        " -R {input.ref}"
        " -nt 5"
        " -input {input.vcf}"
        " -mode SNP"
        " --ts_filter_level 99.80"
        " -recalFile {input.recal}"
        " -tranchesFile {input.tranche}"
        " -o {output.vcf}"


rule gatk3_apply_recalibration_indel:
    """
    Apply the VQSR to the INDELs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/samples/vcf/{panel}.vcf",
        recal="data/samples/vcf/{panel}_vqsr_indel.recal",
        tranche="data/samples/vcf/{panel}_vqsr_indel.tranches",
    output:
        vcf="data/samples/vcf/{panel}_vqsr_indel.vcf",
    log:
        log="data/samples/vcf/{panel}_vqsr_indel.log",
    shell:
        "gatk3"
        " -T ApplyRecalibration"
        " -R {input.ref}"
        " -nt 5"
        " -input {input.vcf}"
        " -mode INDEL"
        " --ts_filter_level 99.0"
        " -recalFile {input.recal}"
        " -tranchesFile {input.tranche}"
        " -o {output.vcf}"

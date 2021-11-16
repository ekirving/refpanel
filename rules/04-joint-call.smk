#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import protected, unpack, temp, expand

from scripts.utils import list_samples

"""
Rules to perform joint genotype calling for the IGSR pipeline

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/1000G_README_2019April10_NYGCjointcalls.pdf 
"""

GATK_NUM_THREADS = 4


wildcard_constraints:
    chr="(chr(\d+|X|Y|M))|(others)",
    type="SNP|INDEL",


def gatk3_genotype_gvcf_input(wildcards):
    return {
        "ref": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "chr": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.{chr}.bed",
        "gvcfs": [
            f"data/source/{source}/gVCF/{sample}.g.vcf.gz" for source, sample in list_samples(config, wildcards.panel)
        ],
    }


# noinspection PyUnresolvedReferences
rule gatk3_genotype_gvcf:
    """
    Jointly call genotypes in all samples
    """
    input:
        unpack(gatk3_genotype_gvcf_input),
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_{chr}.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_{chr}.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}.vcf.log",
    params:
        gvcfs=lambda wildcards, input: " ".join([f"--variant {gvcf}" for gvcf in input.gvcfs]),
    threads: GATK_NUM_THREADS
    shell:
        "gatk3"
        " -T GenotypeGVCFs"
        " -R {input.ref}"
        " -L {input.chr}"
        " --num_threads {threads}"
        " --disable_auto_index_creation_and_locking_when_reading_rods"
        " {params.gvcfs}"
        " -o {output.vcf} 2> {log}"


rule picard_merge_chromosome_vcfs:
    """
    Merge the chromosomes back together so we can do variant recalibration over the whole genome
    """
    input:
        expand("data/panel/{panel}/vcf/{panel}_{chr}.vcf.gz", chr=config["chroms"], allow_missing=True),
    output:
        vcf=protected("data/panel/{panel}/vcf/{panel}_chrALL.vcf.gz"),
        tbi=protected("data/panel/{panel}/vcf/{panel}_chrALL.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_chrALL.vcf.log",
    params:
        vcfs=lambda wildcards, input: [f"INPUT={vcf}" for vcf in input],
    shell:
        "picard"
        " MergeVcfs"
        " {params.vcfs}"
        " OUTPUT={output.vcf} 2> {log}"


rule gatk3_split_variants:
    """
    Split variants into SNPs and INDELs to perform the variant recalibration in parallel
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_chrALL.vcf.gz",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_chrALL_{type}.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_chrALL_{type}.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_chrALL_{type}.vcf.log",
    shell:
        "gatk3 "
        " -T SelectVariants"
        " -R {input.ref}"
        " --variant {input.vcf}"
        " --selectTypeToInclude {wildcards.type}"
        " --out {output.vcf} 2> {log}"


rule gatk3_variant_recalibrator_snp:
    """
    Variant Quality Score Recalibration (VQSR) to assign FILTER status for SNPs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_chrALL_SNP.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_chrALL_SNP.vcf.gz.tbi",
        hap="data/reference/GRCh38/other_mapping_resources/hapmap_3.3.hg38.vcf.gz",
        omni="data/reference/GRCh38/other_mapping_resources/1000G_omni2.5.hg38.vcf.gz",
        snps="data/reference/GRCh38/other_mapping_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="data/reference/GRCh38/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz",
    output:
        recal="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_SNP.recal",
        tranche="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_SNP.tranches",
        plot="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_SNP_plots.R",
    log:
        log="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_SNP.log",
    threads: GATK_NUM_THREADS
    shell:
        "gatk3"
        " -T VariantRecalibrator"
        " -R {input.ref}"
        " --num_threads {threads}"
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
        " -tranche 90.0 2> {log}"


rule gatk3_variant_recalibrator_indel:
    """
    Variant Quality Score Recalibration (VQSR) to assign FILTER status for INDELs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_chrALL_INDEL.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_chrALL_INDEL.vcf.gz.tbi",
        mills="data/reference/GRCh38/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz",
        dbsnp="data/reference/GRCh38/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz",
    output:
        recal="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_INDEL.recal",
        tranche="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_INDEL.tranches",
        plot="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_INDEL_plots.R",
    log:
        log="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_INDEL.log",
    threads: GATK_NUM_THREADS
    shell:
        "gatk3"
        " -T VariantRecalibrator"
        " -R {input.ref}"
        " --num_threads {threads}"
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
        " --maxGaussians 4 2> {log}"


rule gatk3_apply_recalibration_snp:
    """
    Apply the VQSR to the SNPs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_chrALL_SNP.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_chrALL_SNP.vcf.gz.tbi",
        recal="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_SNP.recal",
        tranche="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_SNP.tranches",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_chrALL_vqsr_SNP.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_chrALL_vqsr_SNP.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_SNP.vcf.log",
    threads: GATK_NUM_THREADS
    shell:
        "gatk3"
        " -T ApplyRecalibration"
        " -R {input.ref}"
        " --num_threads {threads}"
        " -input {input.vcf}"
        " -mode SNP"
        " --ts_filter_level 99.80"
        " -recalFile {input.recal}"
        " -tranchesFile {input.tranche}"
        " -o {output.vcf} 2> {log}"


rule gatk3_apply_recalibration_indel:
    """
    Apply the VQSR to the INDELs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_chrALL_INDEL.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_chrALL_INDEL.vcf.gz.tbi",
        recal="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_INDEL.recal",
        tranche="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_INDEL.tranches",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_chrALL_vqsr_INDEL.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_chrALL_vqsr_INDEL.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_INDEL.vcf.log",
    threads: GATK_NUM_THREADS
    shell:
        "gatk3"
        " -T ApplyRecalibration"
        " -R {input.ref}"
        " --num_threads {threads}"
        " -input {input.vcf}"
        " -mode INDEL"
        " --ts_filter_level 99.0"
        " -recalFile {input.recal}"
        " -tranchesFile {input.tranche}"
        " -o {output.vcf} 2> {log}"


rule picard_merge_variant_vcfs:
    """
    Merge the SNP and INDEL VCF files
    """
    input:
        snp="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_SNP.vcf.gz",
        indel="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_INDEL.vcf.gz",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_chrALL_vqsr.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_chrALL_vqsr.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_chrALL_vqsr.vcf.log",
    shell:
        "picard"
        " MergeVcfs"
        " INPUT={input.snp}"
        " INPUT={input.indel}"
        " OUTPUT={output.vcf} 2> {log}"


# noinspection PyTypeChecker
rule bcftools_super_populations:
    """
    Make a bcftools samples file for calculating INFO tags at the super population level

    https://m.ensembl.org/Help/Faq?id=532
    """
    input:
        tsv=lambda wildcards: config["panel"][wildcards.panel]["samples"],
    output:
        tsv="data/panel/{panel}/{panel}-super_populations.tsv",
    params:
        col1=lambda wildcards, input: open(input.tsv).readline().split("\t").index("sample") + 1,
        col2=lambda wildcards, input: open(input.tsv).readline().split("\t").index("super_population") + 1,
    shell:
        r"""awk -v FS="\t" 'NR>1 {{ print ${params.col1} FS ${params.col2} }}' {input.tsv} > {output.tsv}"""


rule bcftools_fill_tags:
    """
    Calculate and fill missing tags and annotations
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_chrALL_vqsr.vcf.gz",
        tsv="data/panel/{panel}/{panel}-super_populations.tsv",
    output:
        vcf=protected("data/panel/{panel}/vcf/{panel}_chrALL_vqsr_annot.vcf.gz"),
        tbi=protected("data/panel/{panel}/vcf/{panel}_chrALL_vqsr_annot.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_annot.vcf.log",
    shell:
        "bcftools +fill-tags {input.vcf} -Oz -o {output.vcf} -- --tags all --samples-file {input.tsv} && "
        "bcftools index --tbi {output.vcf}"


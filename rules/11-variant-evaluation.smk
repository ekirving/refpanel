#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import temp

from scripts.common import JAVA_MEMORY_MB, JAVA_TEMP_DIR

"""
Evaluation of small variant calls (SNPs and small INDELs)
"""


rule bcftools_atomize:
    """
    Split multi-allelic variants into multiple rows (INDELs are already left-normalized)
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz.tbi",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_atom.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_atom.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_atom.vcf.log",
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "( bcftools norm --atomize -Oz -o {output.vcf} {input.vcf} && "
        "  bcftools index --tbi {output.vcf} ) 2> {log}"


rule gatk3_variant_eval:
    """
    Evaluate a callset with GATK VariantEval, comparing to both dbSNP and the NYGC 30x 1000G callset

    https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Evaluate_a_callset_with_VariantEval.md
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_atom.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_atom.vcf.gz.tbi",
        dbsnp="data/reference/GRCh38/dbsnp/GRCh38.dbSNP155.vcf.gz",  # TODO 1000G uses v.151
        comp="data/evaluation/1000g_nygc/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chr}.recalibrated_variants.vcf.gz",
    output:
        evl="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_atom.vcf.eval",
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_atom.vcf.eval.log",
    resources:
        mem_mb=JAVA_MEMORY_MB,
        tmpdir=JAVA_TEMP_DIR,
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -XX:ConcGCThreads=1"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " -T VariantEval"
        " -R {input.ref}"
        " --eval {input.vcf}"
        " --dbsnp {input.dbsnp}"
        " --comp {input.comp}"
        " --known_names '1000G'"
        " --out {output.evl} 2> {log}"


rule giab_genome_stratifications:
    """
    Download the Genome in a Bottle (GIAB) genome stratifications
    """
    output:
        "data/reference/GRCh38/giab-stratifications/v3.0/v3.0-GRCh38-all-stratifications.tsv",
        "data/reference/GRCh38/giab-stratifications/v3.0/v3.0-GRCh38-CMRG-stratifications.tsv",
        "data/reference/GRCh38/giab-stratifications/v3.0/v3.0-GRCh38-stratifications-all-except-genome-specific-stratifications.tsv",
        "data/reference/GRCh38/giab-stratifications/v3.0/v3.0-GRCh38-v4.2.1-stratifications.tsv",
        "data/reference/GRCh38/giab-stratifications/v3.0/v3.0-stratifications-GRCh38-md5s.txt",
    shell:
        "wget "
        " --mirror"
        " --quiet"
        " --no-host-directories"
        " --cut-dirs=6"
        " --directory-prefix=data/reference/GRCh38/giab-stratifications/v3.0/"
        " -o /dev/null"
        " https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v3.0/GRCh38/"


rule giab_na12878_hg001:
    """
    Download the Genome in a Bottle (GIAB) truth-set for NA12878 (HG001)
    """
    output:
        "data/evaluation/GIAB/NA12878_HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.bed",
        "data/evaluation/GIAB/NA12878_HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        "data/evaluation/GIAB/NA12878_HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
    shell:
        "wget "
        " --mirror"
        " --quiet"
        " --no-host-directories"
        " --cut-dirs=6"
        " --directory-prefix=data/evaluation/GIAB/NA12878_HG001/"
        " -o /dev/null"
        " ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/"


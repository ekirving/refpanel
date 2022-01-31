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
        "  bcftools index --tbi {output.vcf} )2> {log}"


rule gatk3_variant_eval:
    """
    Raw variant calls using HaplotypeCaller on a single sample, within a given region and sex-based ploidy

    https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Evaluate_a_callset_with_VariantEval.md
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_atom.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_atom.vcf.gz.tbi",
        dbsnp="data/reference/GRCh38/dbsnp/GRCh38.dbSNP155.vcf.gz",  # TODO 1000G uses v.151
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
        " --out {output.evl} 2> {log}"

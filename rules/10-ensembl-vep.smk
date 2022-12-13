#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import directory, touch

global workflow

"""
Ensembl Variant Effect Predictor (VEP) annotations

https://www.ensembl.org/info/docs/tools/vep 
"""


rule ensembl_vep_install_cache:
    """
    Download and install the VEP cache file for GRCh38
    """
    output:
        dir=directory("data/ensembl/vep/"),
    log:
        log="data/ensembl/vep/vep_install_cache.log",
    benchmark:
        "benchmarks/ensembl_vep_install_cache.tsv"
    conda:
        "../envs/ensembl-vep-105.0.yaml"
    shell:
        "vep_install"
        " --AUTO c"
        " --SPECIES homo_sapiens"
        " --ASSEMBLY GRCh38"
        " --CACHEDIR {output.dir}"
        " --NO_UPDATE &> {log}"


rule ensembl_vep_annotate_vcf:
    """
    Annotate the VCF with predicted effect of each variant on genes, transcripts, and protein sequence, as well as 
    regulatory regions.

    https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        dir="data/ensembl/vep/",
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_{whatshap}_phased.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_{whatshap}_phased.vcf.gz.tbi",
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_{whatshap}_phased_vep.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_{whatshap}_phased_vep.vcf.gz.tbi",
        htm="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_{whatshap}_phased_vep.html",
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_{whatshap}_phased_vep.vcf.log",
    threads: 4
    benchmark:
        "benchmarks/ensembl_vep_annotate_vcf-{panel}-{chr}-{whatshap}.tsv"
    wildcard_constraints:
        whatshap="whatshap|whatshap_trio",
    conda:
        "../envs/ensembl-vep-105.0.yaml"
    shell:
        "vep"
        " --offline"
        " --everything"
        " --species homo_sapiens"
        " --assembly GRCh38"
        " --dir_cache {input.dir}"
        " --vcf"
        " --fork {threads}"
        " --input_file {input.vcf}"
        " --fasta {input.ref}"
        " --stats_file {output.htm}"
        " --compress_output bgzip"
        " --output_file {output.vcf} &> {log} && "
        "tabix -p vcf {output.vcf}"


def panel_predict_variant_effects_input(wildcards):
    """
    Check if the current panel has a pedigree, or not.
    """
    panel = wildcards.panel
    chroms = [chr for chr in config["chroms"] if chr not in ["chrY", "chrM", "others"]]

    if config["panel"][wildcards.panel].get("pedigree") is None:
        vcf = [
            f"data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_whatshap_phased_vep.vcf.gz" for chr in chroms
        ]
    else:
        vcf = [
            f"data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_whatshap_trio_phased_vep.vcf.gz"
            for chr in chroms
        ]

    return vcf


rule panel_predict_variant_effects:
    """
    Perform variant effect prediction for all variants in the panel. 
    """
    input:
        panel_predict_variant_effects_input,
    output:
        touch("data/panel/{panel}/vcf/vep.done"),

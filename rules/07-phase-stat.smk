#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand, touch


global workflow

"""
Rules to perform statistical phasing of the joint-called reference panel

https://whatshap.readthedocs.io/en/latest/guide.html 
"""


# TODO add a version that does not rely on trios
rule shapeit4_phase_vcf:
    """
    Phase the joint-callset, using the pedigree phased VCF as a scaffold and the read-based phase sets as input.

    See https://github.com/odelaneau/shapeit4/issues/17
    """
    input:
        map="data/reference/GRCh38/genetic_maps/shapeit4/{chr}.b38.gmap.gz",
        vcf1="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz",
        tbi1="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_whatshap.vcf.gz.tbi",
        vcf2="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_trios.vcf.gz",
        tbi2="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_trios.vcf.gz.tbi",
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_phased.vcf.gz",
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_phased.vcf.log",
    threads: max(workflow.cores / 4, 8)
    conda:
        "../envs/shapeit-4.2.2.yaml"
    shell:
        "shapeit4"
        " --thread {threads}"
        " --input {input.vcf1}"
        " --scaffold {input.vcf2}"
        " --map {input.map}"
        " --region {wildcards.chr}"
        " --sequencing"
        " --output {output.vcf}"
        " --log {log}"


rule panel_statistical_phasing:
    """
    Perform statistical phasing of all samples a reference panel.
    """
    input:
        expand(
            "data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter_mendel_phased.vcf.gz",
            chr=config["chroms"],
            allow_missing=True,
        ),
    output:
        touch("data/panel/{panel}/vcf/phase.done"),

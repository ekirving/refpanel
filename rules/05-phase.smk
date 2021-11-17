#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

global workflow

"""
Rules to perform variant phasing on a joint-called reference panel  

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/1000G_2020Oct26_NYGC_Phasing_README.pdf 
"""

# use 25% of the total cores
SHAPEIT4_NUM_THREADS = workflow.cores / 4


rule shapeit4_phase_vcf:
    """
    Phase the joint-callset
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_chrALL_vqsr_annot_mendel_filter.vcf.gz",
        map="data/reference/GRCh38/genetic_maps.b38.tar.gz",
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_annot_filter_phased.vcf.gz",
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_annot_filter_phased.vcf.log",
    threads: SHAPEIT4_NUM_THREADS
    conda:
        "../envs/shapeit4.yaml"
    shell:
        "shapeit4"
        " --thread {threads}"
        " --input {input.vcf}"
        " --map {input.map}"
        " --region {wildcards.chr}"
        " --sequencing"
        " --output {output.vcf}"
        " --log {log}"

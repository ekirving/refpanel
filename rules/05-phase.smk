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
BEAGLE_NUM_THREADS = workflow.cores / 4


rule beagle_phase_vcf:
    """
    Phase the joint-callset using Beagle
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr.vcf.gz",
        map="data/reference/GRCh38/plink.GRCh38.map.zip",
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_phased.vcf.gz",
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_phased.vcf.log",
    threads: BEAGLE_NUM_THREADS
    shell:
        "beagle "
        " nthreads={threads}"
        " gt={input.vcf}"
        " map={input.map}"
        " out={output.vcf} 2> {log}"

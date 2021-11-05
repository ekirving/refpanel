#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand, protected, temp

global workflow

"""
Rules to download data files from the International Genome Sample Resource (IGSR)

http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/ 
"""

FTP_1000G_RAW_CALLS = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls/"


rule download_1000g_nygc_gvcf:
    """
    Download GATK HaplotypeCaller gVCFs for each high-coverage NYGC 1000G sample
    """
    output:
        vcf="data/1000G_NYGC/gVCF/{population}/{sample}.g.vcf.gz",
        tbi="data/1000G_NYGC/gVCF/{population}/{sample}.g.vcf.gz.tbi",
    shell:
        "wget --quiet -O {output.vcf} {FTP_1000G_RAW_CALLS}/{wildcards.population}/Sample_{wildcards.sample}/analysis/{wildcards.sample}.haplotypeCalls.er.raw.vcf.gz && "
        "wget --quiet -O {output.tbi} {FTP_1000G_RAW_CALLS}/{wildcards.population}/Sample_{wildcards.sample}/analysis/{wildcards.sample}.haplotypeCalls.er.raw.vcf.gz.tbi"

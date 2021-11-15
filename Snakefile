#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand


configfile: "config.yaml"


# rules to download the data sources
include: "rules/data/1000g.smk"
include: "rules/data/hgdp.smk"
include: "rules/data/sgdp.smk"
include: "rules/data/ggvp.smk"
#
# rules to apply the IGSR genotyping pipeline
include: "rules/01-reference.smk"
include: "rules/02-align.smk"
include: "rules/03-call.smk"
include: "rules/04-joint-call.smk"


# preference download rules, when they are available
ruleorder: sgdp_download_cram > ggvp_download_cram > samtools_cram
ruleorder: tgp_nygc_download_gvcf > hgdp_download_gvcf > gatk3_haplotype_caller


rule all:
    input:
        # "data/source/example/gVCF/HGDP00094.g.vcf.gz",
        # "data/source/example/gVCF/HGDP01356.g.vcf.gz",
        "data/panel/example-panel/vcf/example-panel_chr22_vqsr.vcf.gz",


rule refpanel:
    input:
        expand("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_phased.vcf.gz", panel=config["refpanel"], chr=["22"]),


rule download_data:
    input:
        "data/source/1000g/gVCF/download.done",
        "data/source/hgdp/gVCF/download.done",
        "data/source/sgdp/cram/download.done",
        "data/source/ggvp/cram/download.done",

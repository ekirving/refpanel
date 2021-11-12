#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


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
        # "data/source/sgdp/cram/armenia293.cram",
        # "data/source/ggvp/cram/SC_GMFUL5309722.cram",
        # "data/source/example/bam/ERR1373257.bam",
        # "data/source/example/cram/HGDP00094.cram",
        "data/source/example/bam/HGDP00094_merged_sorted_dedup_recal.table",


rule download_data:
    input:
        "data/source/1000g/gVCF/download.done",
        "data/source/hgdp/gVCF/download.done",
        "data/source/sgdp/cram/download.done",
        "data/source/ggvp/cram/download.done",

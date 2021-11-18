#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
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
include: "rules/05-phase.smk"


# preference download rules, when they are available
ruleorder: sgdp_download_cram > ggvp_download_cram > samtools_cram
ruleorder: tgp_nygc_download_gvcf > hgdp_download_gvcf > gatk3_combine_gvcfs


rule all:
    input:
        "",


rule refpanel:
    input:
        # expand("data/panel/{panel}/vcf/phase.done", panel=config["refpanel"]),
        expand("data/panel/{panel}/vcf/{panel}_chrALL_vqsr_annot_filter_mendel.vcf.gz", panel=config["refpanel"]),


rule download_data:
    input:
        "data/source/1000g/gVCF/download.done",
        "data/source/hgdp/gVCF/download.done",
        "data/source/sgdp/cram/download.done",
        "data/source/ggvp/cram/download.done",


def list_all_gvcfs():
    """List all the gVCF files that need generating from the downloaded CRAM files"""
    files = []

    for source in ["sgdp", "ggvp"]:
        samples = pd.read_table(config["source"][source]["samples"])

        files += [f"data/source/{source}/gVCF/{sample}.g.vcf.gz" for sample in samples["sample"]]

    batch = config.get("batch", 0)
    size = int(len(files) / 3)

    start = (batch - 1) * size
    finish = start + size if batch != 3 else None

    return files[start:finish]


rule generate_gvcfs:
    input:
        list_all_gvcfs(),

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
ruleorder: hgdp_download_cram > sgdp_download_cram > ggvp_filter_iupac_base_codes > samtools_cram
ruleorder: tgp_nygc_download_gvcf > gatk3_combine_ploidy_regions


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
        "data/source/hgdp/cram/download.done",
        "data/source/sgdp/cram/download.done",
        "data/source/ggvp/cram/download.done",


rule merge_1000g:
    input:
        expand("data/source/1000g/gVCF/merged/1000g_{chr}.g.vcf.gz", chr=config["chroms"]),


rule merge_hgdp:
    input:
        expand("data/source/hgdp/gVCF/merged/hgdp_{chr}.g.vcf.gz", chr=config["chroms"]),


rule merge_sgdp:
    input:
        expand("data/source/sgdp/gVCF/merged/sgdp_{chr}.g.vcf.gz", chr=config["chroms"]),


rule merge_ggvp:
    input:
        expand("data/source/ggvp/gVCF/merged/ggvp_{chr}.g.vcf.gz", chr=config["chroms"]),


def list_all_gvcfs():
    """List all the gVCF files that need generating from the downloaded CRAM files"""
    files = []

    for source in ["hgdp"]:  # ["sgdp", "ggvp"]
        samples = pd.read_table(config["source"][source]["samples"])

        files += [f"data/source/{source}/gVCF/{sample}.g.vcf.gz" for sample in samples["sample"]]

    batch = config.get("batch", 0)

    if batch == 0:
        return files

    size = int(len(files) / 3)

    start = (batch - 1) * size
    finish = start + size if batch != 3 else None

    return files[start:finish]


rule generate_gvcfs:
    input:
        list_all_gvcfs(),


rule generate_multisample_chrom_gvcfs:
    input:
        expand("data/source/1000g/gVCF/merged/1000g_{chr}.g.vcf.gz", chr=config["chroms"]),

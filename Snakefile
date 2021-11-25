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
ruleorder: sgdp_download_cram > ggvp_filter_iupac_base_codes > samtools_cram
ruleorder: tgp_nygc_download_gvcf > hgdp_download_gvcf > gatk3_combine_ploidy_regions


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

    import random

    running = ["abh100", "ALB212", "Ayodo_430C", "Ayodo_502C", "Ayodo_81S", "Bu5", "ch113", "DNK05", "DNK07",
               "HGDP00058", "HGDP00090", "HGDP00157", "HGDP00195", "HGDP00208", "HGDP00428", "HGDP00554", "HGDP00597",
               "HGDP00616", "HGDP00656", "HGDP00702", "HGDP00706", "HGDP00737", "HGDP00783", "HGDP00903", "HGDP00915",
               "HGDP00928", "HGDP00936", "HGDP00956", "HGDP01018", "HGDP01028", "HGDP01030", "HGDP01032", "HGDP01095",
               "HGDP01153", "HGDP01198", "HGDP01203", "HGDP01211", "HGDP01242", "HGDP01246", "HGDP01297", "HGDP01306",
               "HGDP01315", "HGDP01335", "HGDP01345", "HGDP01355", "HGDP01401", "HGDP01417", "I3", "Igor20",
               "Jordan445", "KD4", "Kor82", "mg27", "mixa0105", "mixe0007", "NA11201", "NA15202", "Nesk_25", "Nlk1",
               "R6", "SA0342", "SA0722", "Sir19", "Sir40", "TZ-11", "Ul5", "Y4"]

    for source in ["sgdp", "ggvp"]:
        samples = pd.read_table(config["source"][source]["samples"])

        files += [
            f"data/source/{source}/gVCF/{sample}.g.vcf.gz" for sample in samples["sample"] if sample not in running
        ]

    batch = config.get("batch", 0)
    size = int(len(files) / 3)

    start = (batch - 1) * size
    finish = start + size if batch != 3 else None

    random.seed(4)
    random.shuffle(files)

    return files[start:finish]


rule generate_gvcfs:
    input:
        list_all_gvcfs(),

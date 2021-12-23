#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import expand


configfile: "config.yaml"


# TODO add benchmarks to all rules
# TODO implement the "best practices" criteria https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html


# rules to download the data sources
include: "rules/data/1000g.smk"
include: "rules/data/hgdp.smk"
include: "rules/data/sgdp.smk"
include: "rules/data/ggvp.smk"
#
# rules for the genotyping and phasing pipeline
include: "rules/01-reference.smk"
include: "rules/02-align.smk"
include: "rules/03-qc-metrics.smk"
include: "rules/04-call.smk"
include: "rules/05-merge-calls.smk"
include: "rules/06-joint-call.smk"
include: "rules/07-phase-reads.smk"
include: "rules/08-phase-trios.smk"
include: "rules/09-phase-stat.smk"
include: "rules/10-ensembl-vep.smk"


# preference download rules, when they are available
ruleorder: tgp_nygc_download_cram > hgdp_filter_iupac_base_codes > sgdp_download_cram > ggvp_standardise_sample_names > samtools_cram
ruleorder: tgp_nygc_download_gvcf > gatk3_combine_ploidy_regions


rule all:
    input:
        "",


rule refpanel:
    input:
        # run the pipeline end-to-end
        expand("data/panel/{panel}/vcf/phase.done", panel=config["refpanel"]),


rule download_data:
    input:
        # download CRAMs for all data sources, and gVCFs for 1000g
        expand("data/source/{source}/cram/download.done", source=["1000g", "hgdp", "sgdp", "ggvp"]),
        expand("data/source/{source}/gVCF/download.done", source=["1000g"]),

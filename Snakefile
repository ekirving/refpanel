#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand

from scripts.common import list_sources


configfile: "config.yaml"


# rules to download the data sources
include: "rules/data/1000g.smk"
include: "rules/data/1000g_nygc.smk"
include: "rules/data/hgdp.smk"
include: "rules/data/sgdp.smk"
include: "rules/data/ggvp.smk"
include: "rules/data/ena.smk"
include: "rules/data/phase10x.smk"
#
# rules for the genotyping and phasing pipeline
include: "rules/01-reference.smk"
include: "rules/02-align.smk"
include: "rules/03-qc-metrics.smk"
include: "rules/04-call.smk"
include: "rules/05-merge-samples.smk"
include: "rules/06-joint-call.smk"
include: "rules/07-phase-reads.smk"
include: "rules/08-phase-trios.smk"
include: "rules/09-phase-stat.smk"
include: "rules/10-ensembl-vep.smk"
#
# rules for evaluating the reference panel
include: "rules/11-variant-evaluation.smk"


# preference download rules, when they are available
ruleorder: tgp_nygc_download_cram > hgdp_standardise_sample_names > sgdp_standardise_sample_names > ggvp_standardise_sample_names > samtools_cram
ruleorder: tgp_nygc_download_gvcf > gatk3_combine_ploidy_regions


rule all:
    input:
        "",


rule refpanel:
    input:
        # run the pipeline end-to-end
        expand("data/panel/{panel}/vcf/phase.done", panel=config["refpanel"]),


# get the list of data sources in the default reference panel
SOURCES = list_sources(config, config["refpanel"])

# the list of data sources that are hosted on the European Nucleotide Archive (ENA) -- excluding GGVP, which has CRAMs
ENA_SOURCES = [source for source in SOURCES if config["source"][source].get("ena_ftp", False) and source != "ggvp"]


rule download_data:
    input:
        # download sample data for all data sources (the starting point varies between sources)
        expand("data/source/{source}/fastq/download.done", source=ENA_SOURCES),
        expand("data/source/{source}/cram/download.done", source=["1000g", "hgdp", "sgdp", "ggvp", "PRJNA76713"]),
        expand("data/source/{source}/gVCF/download.done", source=["1000g"]),
        expand("data/source/{source}/gVCF/phase10x/download.done", source=["hgdp", "appg"]),


rule align_sources:
    input:
        # align all samples in all data sources
        expand("data/source/{source}/cram/align.done", source=SOURCES),


rule metrics_sources:
    input:
        # calculate alignment metrics for all samples in all data sources
        expand("data/source/{source}/cram/metrics.done", source=SOURCES),


rule call_sources:
    input:
        # call all samples in all data sources
        expand("data/source/{source}/gVCF/call.done", source=SOURCES),


rule merge_sources:
    input:
        # merge all samples in each data source
        expand("data/source/{source}/gVCF/merge.done", source=SOURCES),


rule merge_panel:
    input:
        # merge all data sources in the reference panel
        expand("data/panel/{panel}/gVCF/merge.done", panel=config["refpanel"]),


rule call_panel:
    input:
        # joint-call all samples in the reference panel
        expand("data/panel/{panel}/vcf/joint-call.done", panel=config["refpanel"]),


rule phase_panel:
    input:
        # phase all samples in the reference panel
        expand("data/panel/{panel}/vcf/phase.done", panel=config["refpanel"]),


rule predict_effects:
    input:
        # predict variant effects for all variants in the panel
        expand("data/panel/{panel}/vcf/vep.done", panel=config["refpanel"]),

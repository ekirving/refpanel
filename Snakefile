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
include: "rules/data/ena.smk"
include: "rules/data/phase10x.smk"
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
ruleorder: tgp_nygc_download_cram > hgdp_standardise_sample_names > sgdp_standardise_sample_names > ggvp_standardise_sample_names > samtools_cram
ruleorder: tgp_nygc_download_gvcf > gatk3_combine_ploidy_regions


rule all:
    input:
        "",


rule refpanel:
    input:
        # run the pipeline end-to-end
        expand("data/panel/{panel}/vcf/phase.done", panel=config["refpanel"]),


ENA_SOURCES = [
    source for source in config["source"] if config["source"][source].get("ena_ftp", False) and source != "ggvp"
]


rule download_data:
    input:
        # download CRAMs for 1000G, HGDP, SGDP, GGVP; gVCFs for 1000G; 10x gVCFs for HGDP and APPG; and FASTQs for everything else
        expand("data/source/{source}/cram/download.done", source=["1000g", "hgdp", "sgdp", "ggvp", "PRJNA76713"]),
        expand("data/source/{source}/gVCF/download.done", source=["1000g"]),
        expand("data/source/{source}/gVCF/phase10x/download.done", source=["hgdp", "appg"]),
        expand("data/source/{source}/fastq/download.done", source=ENA_SOURCES),


rule align_sources:
    input:
        # align all samples in all data sources
        expand("data/source/{source}/cram/align.done", source=config["source"]),


rule metrics_sources:
    input:
        # calculate alignment metrics for all samples in all data sources
        expand("data/source/{source}/cram/metrics.done", source=config["source"]),


rule call_sources:
    input:
        # call all samples in all data sources
        expand("data/source/{source}/gVCF/call.done", source=config["source"]),


rule merge_sources:
    input:
        # merge all samples in each data source
        expand("data/source/{source}/gVCF/merge.done", source=config["source"]),

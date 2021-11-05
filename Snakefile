#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


# configfile: "config.yaml"


include: "rules/igsr/ref.smk"
include: "rules/igsr/align.smk"
include: "rules/igsr/call.smk"
include: "rules/igsr/download.smk"


def download_1000g_nygc(_):
    files = []
    with open("data/1000G_samples.txt ") as fin:
        for sample in fin:
            population, sample = sample.split()
            files.append(f"data/1000G_NYGC/gVCF/{population}/{sample}.g.vcf.gz")

    return files


rule all:
    input:
        download_1000g_nygc,

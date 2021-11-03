#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


# configfile: "config.yaml"


include: "rules/igsr/ref.smk"
include: "rules/igsr/align.smk"


rule all:
    input:
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",

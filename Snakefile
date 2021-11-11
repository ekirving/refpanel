#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


# include: "rules/01-reference.smk"
# include: "rules/02-align.smk"
# include: "rules/03-call.smk"
include: "rules/download/hgdp.smk"
include: "rules/download/sgdp.smk"
include: "rules/download/tgp_nygc.smk"


rule all:
    input:
        "data/hgdp/gVCF/download.done",
        "data/sgdp/cram/download.done",
        "data/tpg_nygc/gVCF/download.done",

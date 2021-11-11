#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


include: "rules/01-reference.smk"
include: "rules/02-align.smk"
include: "rules/03-call.smk"
include: "rules/download/hgdp.smk"
include: "rules/download/nygc_1000g.smk"
include: "rules/download/sgdp.smk"


rule all:
    input:
        "data/1000g_nygc/gVCF/download.done",
        "data/hgdp/gVCF/download.done",
        "data/sgdp/cram/download.done",

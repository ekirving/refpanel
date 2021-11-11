#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


# include: "rules/01-reference.smk"
# include: "rules/02-align.smk"
# include: "rules/03-call.smk"
# include: "rules/04-joint-call.smk"
include: "rules/download/1000g.smk"
include: "rules/download/sgdp.smk"
include: "rules/download/hgdp.smk"
include: "rules/download/ggvp.smk"


rule all:
    input:
        "data/source/1000g/gVCF/download.done",
        "data/source/sgdp/cram/download.done",
        "data/source/hgdp/gVCF/download.done",


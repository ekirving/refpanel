#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


include: "rules/igsr/ref.smk"
include: "rules/igsr/align.smk"
include: "rules/igsr/call.smk"
include: "rules/igsr/nygc_1000g.smk"
include: "rules/igsr/hgdp.smk"


rule all:
    input:
        "data/1000G_NYGC/gVCF/download.done",
        "data/HGDP/gVCF/download.done",

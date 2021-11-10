#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch

"""
Rules to download data files for the Human Genome Diversity Project (HGDP)

https://www.internationalgenome.org/data-portal/data-collection/hgdp
"""


FTP_HGDP = "ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/gVCFs"


wildcard_constraints:
    ext="([a-z]+\.)*[a-z]+",


rule hgdp_download_gvcf:
    """
    Download GATK HaplotypeCaller gVCFs for each HGDP sample
    """
    output:
        vcf="data/HGDP/gVCF/{sample}.g.{ext}",
    resources:
        sanger_ftp=1,
    shell:
        "wget --quiet -O {output.vcf} {FTP_HGDP}/{wildcards.sample}.hgdp_wgs.20190516.{wildcards.ext}"


def hgdp_list_all_gvcf():
    samples = pd.read_table("data/HGDP/hgdp_wgs.20190516.metadata.txt")

    files = [
        [f"data/HGDP/gVCF/{sample}.g.vcf.gz", f"data/HGDP/gVCF/{sample}.g.vcf.gz.tbi"] for sample in samples["sample"]
    ]

    return files


rule hgdp_download_all_gvcf:
    input:
        hgdp_list_all_gvcf(),
    output:
        touch("data/HGDP/gVCF/download.done"),

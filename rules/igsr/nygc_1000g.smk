#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch

"""
Rules to download data files from the New York Genome Center (NYGC) high-coverage 1000 Genomes Project datset

https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
"""

FTP_NYGC_1000G = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls_updated"


wildcard_constraints:
    ext="vcf.gz(.tbi)?",


rule nygc_1000g_md5:
    """
    Make an md5 checksum file for validating the 1000G data
    """
    input:
        man="data/1000G_NYGC/20211105_NYGC_GATK_raw_calls_updated_manifest.txt",
    output:
        md5="data/1000G_NYGC/gVCF/{sample}.g.{ext}.md5",
    params:
        rgx=r"{sample}.haplotypeCalls.er.raw.{ext}\s",
        file="data/1000G_NYGC/gVCF/{sample}.g.{ext}",
    shell:
        """grep -P '{params.rgx}' {input.man} | awk '{{ print $3" {params.file}" }}' > {output.md5}"""


rule nygc_1000g_download_gvcf:
    """
    Download GATK HaplotypeCaller gVCFs for each high-coverage NYGC 1000G sample
    """
    input:
        md5="data/1000G_NYGC/gVCF/{sample}.g.{ext}.md5",
    output:
        vcf="data/1000G_NYGC/gVCF/{sample}.g.{ext}",
    resources:
        ebi_ftp=1,
    shell:
        "wget --quiet -O {output.vcf} {FTP_NYGC_1000G}/{wildcards.sample}.haplotypeCalls.er.raw.{wildcards.ext} && "
        "md5sum --status --check {input.md5}"


def nygc_1000g_list_all_gvcf():
    samples = pd.read_table("data/1000G_NYGC/igsr-1000_genomes_30x_on_grch38.tsv")

    files = [
        [f"data/1000G_NYGC/gVCF/{sample}.g.vcf.gz", f"data/1000G_NYGC/gVCF/{sample}.g.vcf.gz.tbi"]
        for sample in samples["Sample name"]
    ]

    return files


rule nygc_1000g_download_all_gvcf:
    input:
        nygc_1000g_list_all_gvcf(),
    output:
        touch("data/1000G_NYGC/gVCF/download.done"),

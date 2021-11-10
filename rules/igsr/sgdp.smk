#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch

"""
Rules to download data files for the Simons Genome Diversity Project (SGDP)

https://www.internationalgenome.org/data-portal/data-collection/sgdp
"""


wildcard_constraints:
    ext="([a-z]+\.)*[a-z]+",


rule sgdp_md5:
    """
    Make an md5 checksum file for validating the SGDP data
    """
    input:
        man="data/SGDP/simons_diversity_data.GRCh38DH.alignment.index",
    output:
        md5="data/SGDP/cram/{sample}.{ext}.md5",
    params:
        file="data/SGDP/gVCF/{sample}.g.{ext}",
    shell:
        """grep -P '/{wildcards.sample}/' {input.man} | awk '{{ print $3" {params.file}" }}' > {output.md5}"""


rule sgdp_download_gvcf:
    """
    Download GATK HaplotypeCaller gVCFs for each high-coverage NYGC 1000G sample
    """
    input:
        md5="data/SGDP/gVCF/{sample}.g.{ext}.md5",
    output:
        vcf="data/SGDP/gVCF/{sample}.g.{ext}",
    resources:
        ebi_ftp=1,
    shell:
        "wget --quiet -O {output.vcf} {FTP_sgdp}/{wildcards.sample}.haplotypeCalls.er.raw.{wildcards.ext} && "
        "md5sum --status --check {input.md5}"


def sgdp_list_all_gvcf():
    samples = pd.read_table("data/SGDP/igsr-1000_genomes_30x_on_grch38.tsv")

    files = [
        [f"data/SGDP/gVCF/{sample}.g.vcf.gz", f"data/SGDP/gVCF/{sample}.g.vcf.gz.tbi"]
        for sample in samples["Sample name"]
    ]

    return files


rule sgdp_download_all_gvcf:
    input:
        sgdp_list_all_gvcf(),
    output:
        touch("data/SGDP/gVCF/download.done"),

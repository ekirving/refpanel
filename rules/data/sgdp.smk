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
    ext="cram(.crai)?",


rule sgdp_md5:
    """
    Make an md5 checksum file for validating the SGDP data
    """
    input:
        man="data/source/sgdp/simons_diversity_data.GRCh38DH.alignment.tsv",
    output:
        md5="data/source/sgdp/cram/{sample}.{ext}.md5",
    params:
        file="data/source/sgdp/cram/{sample}.{ext}",
        col=lambda wildcards: 4 if "crai" in wildcards.ext else 2,
    shell:
        r"""grep -P '\t{wildcards.sample}\t' {input.man} | awk '{{ print ${params.col}" {params.file}" }}' > {output.md5}"""


rule sgdp_download_cram:
    """
    Download bwa-mem CRAM files for each fully-public SGDP sample
    """
    input:
        man="data/source/sgdp/simons_diversity_data.GRCh38DH.alignment.tsv",
        md5="data/source/sgdp/cram/{sample}.{ext}.md5",
    output:
        cram="data/source/sgdp/cram/{sample}.{ext}",
    params:
        col=lambda wildcards: 3 if "crai" in wildcards.ext else 1,
    resources:
        ebi_ftp=1,
    shell:
        r"grep -P '\t{wildcards.sample}\t' {input.man} | awk '{{ print ${params.col} }}' | "
        r"xargs wget --quiet -O {output.cram} -o /dev/null && "
        r"md5sum --status --check {input.md5}"


def sgdp_list_all_cram():
    samples = pd.read_table(config["source"]["sgdp"]["samples"])

    files = [
        [f"data/source/sgdp/cram/{sample}.cram", f"data/source/sgdp/cram/{sample}.cram.crai"]
        for sample in samples["sample"]
    ]

    return files


rule sgdp_download_all_cram:
    input:
        sgdp_list_all_cram(),
    output:
        touch("data/source/sgdp/cram/download.done"),

#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch

"""
Rules to download data files for the Gambian Genome Variation Project (GGVP)

https://www.internationalgenome.org/data-portal/data-collection/ggvp-grch38
"""


wildcard_constraints:
    ext="cram(.crai)?",


rule ggvp_md5:
    """
    Make an md5 checksum file for validating the GGVP data
    """
    input:
        man="data/source/ggvp/gambian_genome_variation_project.GRCh38DH.alignment.tsv",
    output:
        md5="data/source/ggvp/cram/{sample}.{ext}.md5",
    params:
        file="data/source/ggvp/cram/{sample}.{ext}",
        col=lambda wildcards: 4 if "crai" in wildcards.ext else 2,
    shell:
        r"""grep -P '\t{wildcards.sample}\t' {input.man} | awk '{{ print ${params.col}" {params.file}" }}' > {output.md5}"""


rule ggvp_download_cram:
    """
    Download bwa-mem CRAM files for each fully-public GGVP sample
    """
    input:
        man="data/source/ggvp/gambian_genome_variation_project.GRCh38DH.alignment.tsv",
        md5="data/source/ggvp/cram/{sample}.{ext}.md5",
    output:
        cram="data/source/ggvp/cram/{sample}.{ext}",
    params:
        col=lambda wildcards: 3 if "crai" in wildcards.ext else 1,
    resources:
        ebi_ftp=1,
    shell:
        r"grep -P '\t{wildcards.sample}\t' {input.man} | awk '{{ print ${params.col} }}' | "
        r"xargs wget --quiet -O {output.cram} && "
        r"md5sum --status --check {input.md5}"


def ggvp_list_all_cram():
    samples = pd.read_table("data/source/ggvp/gambian_genome_variation_project.GRCh38DH.alignment.tsv")

    files = [
        [f"data/source/ggvp/cram/{sample}.cram", f"data/source/ggvp/cram/{sample}.cram.crai"]
        for sample in samples["Sample"]
    ]

    return files


rule ggvp_download_all_cram:
    input:
        ggvp_list_all_cram(),
    output:
        touch("data/source/ggvp/cram/download.done"),

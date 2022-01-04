#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch, temp, ancient


"""
Rules to download data files for the Arabian Peninsula Population Genomic (APPG) study

https://www.ebi.ac.uk/ena/browser/view/PRJEB28504
"""


wildcard_constraints:
    sample="[\w-]+",


rule appg_md5:
    """
    Make md5 checksum files for validating the APPG FASTQs
    """
    input:
        man=ancient("data/source/appg/appg-accessions.tsv"),
    output:
        r1_md5=temp("data/source/appg/fastq/{accession}_r1.md5"),
        r2_md5=temp("data/source/appg/fastq/{accession}_r2.md5"),
    shell:
        r"""awk -v FS="\t" '$2=="{wildcards.accession}" {{ print $9,  $13 }}' {input.man} > {output.r1_md5} && """
        r"""awk -v FS="\t" '$2=="{wildcards.accession}" {{ print $10, $14 }}' {input.man} > {output.r2_md5}"""


rule appg_download_fastq:
    """
    Download FASTQ files for each APPG accession
    """
    input:
        man=ancient("data/source/appg/appg-accessions.tsv"),
        md5="data/source/appg/fastq/{accession}_{pair}.md5",
    output:
        fastq="data/source/appg/fastq/{accession}_{pair}.fastq.gz",
    params:
        col=lambda wildcards: 11 if wildcards.pair == "r1" else 12,
    resources:
        ebi_ftp=1,
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        r"""awk -v FS="\t" '$2=="{wildcards.accession}" {{ print ${params.col} }}' {input.man} | """
        "xargs wget --quiet -O {output.fastq} -o /dev/null && "
        "md5sum --status --check {input.md5}"


def appg_list_all_fastq():
    accessions = pd.read_table(config["source"]["appg"]["accessions"])

    files = [
        [f"data/source/appg/fastq/{accession}_r1.fastq.gz", f"data/source/appg/fastq/{accession}_r2.fastq.gz"]
        for accession in accessions["accession"]
    ]

    return files


rule appg_download_all_fastq:
    input:
        appg_list_all_fastq(),
    output:
        touch("data/source/appg/fastq/download.done"),

#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch, temp, ancient


"""
Rules to download FASTQ files from the ENA

https://www.ebi.ac.uk/ena/browser/
"""

# data sources that store their data in the ENA
ENA_SOURCES = [source for source in config["source"] if config["source"][source].get("ena_ftp", False)]


wildcard_constraints:
    source="[\w-]+",
    sample="[\w-]+",
    pair="r1|r2|se",


def get_column_indexes(file_path, columns, sep="\t", offset=1):
    """
    Peek at the first line of the input file to determine the indexes of column headers
    """
    with open(file_path) as fin:
        header = fin.readline().strip().split(sep)

    return [header.index(col) + offset for col in columns]


rule ena_fastq_md5:
    """
    Make an md5 checksum file for validating the downloaded FASTQ file
    """
    input:
        man=lambda wildcards: ancient(config["source"][wildcards.source]["accessions"]),
    output:
        md5=temp("data/source/{source}/fastq/{accession}_{pair}.md5"),
    params:
        idx=lambda wildcards, input: get_column_indexes(
            input.man, ["accession", f"fastq_{wildcards.pair}_md5", f"fastq_{wildcards.pair}"]
        ),
    benchmark:
        "benchmarks/ena_fastq_md5-{source}-{accession}-{pair}.tsv"
    shell:
        r"""awk -v FS="\t" '${params.idx[0]}=="{wildcards.accession}" {{ print ${params.idx[1]}, ${params.idx[2]} }}' {input.man} > {output.md5}"""


rule ena_fastq_download:
    """
    Download FASTQ files for each ENA accession
    """
    input:
        man=lambda wildcards: ancient(config["source"][wildcards.source]["accessions"]),
        md5="data/source/{source}/fastq/{accession}_{pair}.md5",
    output:
        fastq="data/source/{source}/fastq/{accession}_{pair}.fastq.gz",
    params:
        idx=lambda wildcards, input: get_column_indexes(input.man, ["accession", f"fastq_{wildcards.pair}_ftp"]),
    resources:
        ftp=1,
    wildcard_constraints:
        source="|".join(ENA_SOURCES),
    benchmark:
        "benchmarks/ena_fastq_download-{source}-{accession}-{pair}.tsv"
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        r"""awk -v FS="\t" '${params.idx[0]}=="{wildcards.accession}" {{ print ${params.idx[1]} }}' {input.man} | """
        "xargs wget --quiet -O {output.fastq} -o /dev/null && "
        "md5sum --status --check {input.md5}"


def ena_list_all_fastq(wildcards):
    """
    List all the FASTQ files for a given data source
    """
    accessions = pd.read_table(config["source"][wildcards.source]["accessions"])

    fastqs = []

    for _, accession in accessions.iterrows():
        fastqs += [accession.get("fastq_se"), accession.get("fastq_r1"), accession.get("fastq_r2")]

    return [fq for fq in fastqs if isinstance(fq, str)]


rule ena_download_all_fastq:
    input:
        ena_list_all_fastq,
    output:
        touch("data/source/{source}/fastq/download.done"),

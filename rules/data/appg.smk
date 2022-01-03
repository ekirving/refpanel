#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch, temp, ancient

from scripts.common import get_accession_sample

"""
Rules to download data files for the Arabian Peninsula Population Genomic (APPG) study

https://www.ebi.ac.uk/ena/browser/view/PRJEB28504
"""


wildcard_constraints:
    sample="[\w-]+",


rule appg_md5:
    """
    Make an md5 checksum file for validating the APPG CRAMs
    """
    input:
        man=ancient("data/source/appg/appg-accessions.tsv"),
    output:
        md5=temp("data/source/appg/cram/{accession}.raw.cram.md5"),
    params:
        file="data/source/appg/cram/{accession}.raw.cram",
    shell:
        r"""awk -v FS="\t" '$2=="{wildcards.accession}" {{ print $4" {params.file}" }}' {input.man} > {output.md5}"""


rule appg_download_cram:
    """
    Download CRAM files for each APPG accession
    """
    input:
        man=ancient("data/source/appg/appg-accessions.tsv"),
        md5="data/source/appg/cram/{accession}.raw.cram.md5",
    output:
        cram=temp("data/source/appg/cram/{accession}.raw.cram"),
        crai=temp("data/source/appg/cram/{accession}.raw.cram.crai"),
    resources:
        ebi_ftp=1,
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        r"""awk -v FS="\t" '$2=="{wildcards.accession}" {{ print $3 }}' {input.man} | """
        "xargs wget --quiet -O {output.cram} -o /dev/null && "
        "md5sum --status --check {input.md5} && "
        "smatools index {output.cram}"


rule appg_standardise_sample_names:
    """
    Standardise the sample naming in the read groups

    e.g. replace `ERS2705196` with `APPG7555879`
    """
    input:
        cram="data/source/appg/cram/{accession}.raw.cram",
        crai="data/source/appg/cram/{accession}.raw.cram.crai",
    output:
        cram="data/source/appg/cram/{accession}.cram",
        crai="data/source/appg/cram/{accession}.cram.crai",
    params:
        sample=lambda wildcards: get_accession_sample(config, "appg", wildcards.accession),
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        "samtools reheader --command \"sed 's/SM:[^\\t]*/SM:{params.sample}/g'\" {input.cram} > {output.cram} && "
        "samtools index {output.cram}"


def appg_list_all_crams():
    samples = pd.read_table(config["source"]["appg"]["samples"])

    files = [
        [f"data/source/appg/cram/{sample}.cram", f"data/source/appg/cram/{sample}.cram.crai"]
        for sample in samples["sample"]
    ]

    return files


rule appg_download_all_cram:
    input:
        appg_list_all_crams(),
    output:
        touch("data/source/appg/cram/download.done"),

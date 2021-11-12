#!/usr/bin/env python
# -*- coding",  # utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd


def fastq_path(config, source, accession, pair):
    """
    Get the path to the FASTQ file for the given accession
    """
    accessions = pd.read_table(config["source"][source]["accessions"]).set_index("accession", drop=False)

    path = accessions.loc[accession]["fastq_r1" if pair == "r1" else "fastq_r2"]

    return path


def read_group(config, source, accession):
    """
    Compose the read group headers by extracting the mapped fields from the sample metadata sheet

    For @RG tag definitions, see https://support.sentieon.com/appnotes/read_groups/
    """
    accessions = pd.read_table(config["source"][source]["accessions"]).set_index("accession", drop=False)

    # get the column name mapping for the read group tags
    mapping = config["source"][source]["readgroup"]

    # apply the mapping
    tags = [f"{tag}:{accessions.loc[accession][col]}" for tag, col in mapping.items()]

    return r"@RG\t" + r"\t".join(tags)


def list_accessions(config, source, sample):
    """
    Get the list of accession codes for the given sample
    """
    accessions = pd.read_table(config["source"][source]["accessions"]).set_index("sample", drop=False)

    return accessions.loc[[sample]]["accession"].tolist()


def list_samples(config, panel):
    """
    Get the list of samples contained in this reference panel
    """
    samples = pd.read_table(config["panel"][panel]["samples"]).set_index("sample", drop=False)

    return samples[["source", "sample"]].to_records(index=False).tolist()

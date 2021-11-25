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


def sample_sex(config, source, sample):
    """
    Get the sex of the sample from the metadata sheet (i.e., M or F)
    """
    samples = pd.read_table(config["source"][source]["samples"]).set_index("sample", drop=False)

    return samples.loc[sample]["sex"][0].upper()


def list_samples(config, source):
    """
    Get a list of samples for the given data source
    """
    metadata = pd.read_table(config["source"][source]["samples"]).set_index("sample", drop=False)

    return metadata["sample"].unique().tolist()


def list_sources(config, panel):
    """
    Get a list of sources for the given reference panel
    """
    metadata = pd.read_table(config["panel"][panel]["samples"]).set_index("sample", drop=False)

    return metadata["source"].unique().tolist()

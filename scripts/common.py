#!/usr/bin/env python
# -*- coding",  # utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import numpy as np
import pandas as pd
from psutil import virtual_memory

"""
Helper functions and common settings to support the main workflow
"""

# maximum available RAM
MAX_MEM_MB = int(virtual_memory().total / 1024 ** 2) - 1024

# GATK / JAVA default settings
JAVA_MEMORY_MB = 8 * 1024

# number of threads to assign GATK
GATK_NUM_THREADS = 4

# the maximum number of samples to merge together at one time with CombineGVCFs (to prevent massive memory usage)
GATK_BATCH_SIZE = 200

# the default `/tmp` partition is too small
JAVA_TEMP_DIR = "./tmp/"

# the hard limit for the maximum open file descriptors
MAX_OPEN_FILES = 99999


def fastq_path(config, source, accession, pair="r1"):
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

    # get the column name mapping for the read group tags, or use the default mapping
    mapping = config["source"][source].get("readgroup", config["readgroup"])

    # apply the mapping
    tags = [f"{tag}:{accessions.loc[accession][col]}" for tag, col in mapping.items()]

    return r"@RG\t" + r"\t".join(tags)


def list_accessions(config, source, sample):
    """
    Get the list of accession codes and their paired-end status
    """
    accessions = (
        pd.read_table(config["source"][source]["accessions"]).set_index("sample", drop=False).replace(np.nan, "")
    )

    # is the library paired-end or single-end?
    accessions["paired"] = accessions["fastq_r2"].map(lambda fq: fq != "")

    return accessions.loc[[sample]][["accession", "paired"]].to_records(index=False).tolist()


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


def list_source_samples(config, panel):
    """
    Get a list of (source, sample) for the given reference panel
    """
    metadata = pd.read_table(config["panel"][panel]["samples"]).set_index("sample", drop=False)

    return metadata[["source", "sample"]].to_records(index=False).tolist()


def list_families(config, panel):
    """
    Get a list of families from the reference panel pedigree
    """
    pedigree = pd.read_table(
        config["panel"][panel]["pedigree"],
        delimiter=" ",
        usecols=range(4),
        names=["family", "child", "father", "mother"],
        header=None,
    )

    return pedigree["family"].unique().tolist()


def list_family_children(config, panel, family):
    """
    Get a list of children belonging to a specific family from the pedigree
    """
    pedigree = pd.read_table(
        config["panel"][panel]["pedigree"],
        delimiter=" ",
        usecols=range(4),
        names=["family", "child", "father", "mother"],
        header=None,
    ).set_index("family", drop=False)

    return pedigree.loc[[family]]["child"].unique().tolist()

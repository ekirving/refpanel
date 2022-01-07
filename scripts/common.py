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


def fastq_path(config, source, accession, pair="se"):
    """
    Get the path to the FASTQ file for the given accession
    """
    accession = pd.read_table(config["source"][source]["accessions"]).set_index("accession", drop=False).loc[accession]

    return accession[f"fastq_{pair}"]


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
    meta = pd.read_table(config["source"][source]["accessions"]).set_index("sample", drop=False).replace(np.nan, "")

    if "layout" not in meta:
        # determine the library layout
        meta["layout"] = ["SINGLE" if row.get("fastq_r2", "") == "" else "PAIRED" for _, row in meta.iterrows()]

    # does the accession contain unpaired mates
    meta["unpaired"] = [row["layout"] == "PAIRED" and row.get("fastq_se", "") != "" for _, row in meta.iterrows()]

    return meta.loc[[sample]][["accession", "layout", "unpaired"]].to_records(index=False).tolist()


def sample_sex(config, source, sample):
    """
    Get the sex of the sample from the metadata sheet (i.e., M or F)
    """
    samples = pd.read_table(config["source"][source]["samples"]).set_index("sample", drop=False)

    sex = samples.loc[sample].get("sex")[0].upper()

    # sanity check that the metadata is well formed
    assert sex in {"M", "F"}

    return sex


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
    Get a list of (source, sample, prephased) for the given reference panel
    """
    metadata = pd.read_table(config["panel"][panel]["samples"]).set_index("sample", drop=False)

    prephased = []

    # do any of the data sources contain pre-phased linked-reads
    for source in metadata["source"].unique():
        if file_path := config["source"][source].get("prephased"):
            prephased += pd.read_table(file_path).set_index("sample", drop=False)["sample"].tolist()

    metadata["prephased"] = [sample in prephased for sample in metadata["sample"]]

    return metadata[["source", "sample", "prephased"]].to_records(index=False).tolist()


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

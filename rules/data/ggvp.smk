#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch, temp, ancient

"""
Rules to download data files for the Gambian Genome Variation Project (GGVP)

https://www.internationalgenome.org/data-portal/data-collection/ggvp-grch38
"""

# the CRAM files for these 6 samples are missing from the IGSR website
MISSING_CRAM = [
    "SC_GMFUL5306362",
    "SC_GMFUL5306370",
    "SC_GMFUL5306378",
    "SC_GMFUL5306403",
    "SC_GMFUL5306428",
    "SC_GMMAN5482315",
]


wildcard_constraints:
    sample=f"((?!{'|'.join(MISSING_CRAM)})[\w-])+",
    ext="cram(.crai)?",


rule ggvp_md5:
    """
    Make an md5 checksum file for validating the GGVP data
    """
    input:
        man=ancient("data/source/ggvp/gambian_genome_variation_project.GRCh38DH.alignment.tsv"),
    output:
        md5=temp("data/source/ggvp/cram/{sample}.raw.{ext}.md5"),
    params:
        file="data/source/ggvp/cram/{sample}.raw.{ext}",
        col=lambda wildcards: 4 if "crai" in wildcards.ext else 2,
    benchmark:
        "benchmarks/ggvp_md5-{sample}-{ext}.tsv"
    shell:
        r"""grep -P '\t{wildcards.sample}\t' {input.man} | awk '{{ print ${params.col}" {params.file}" }}' > {output.md5}"""


rule ggvp_download_cram:
    """
    Download bwa-mem CRAM files for each fully-public GGVP sample
    """
    input:
        man=ancient("data/source/ggvp/gambian_genome_variation_project.GRCh38DH.alignment.tsv"),
        md5="data/source/ggvp/cram/{sample}.raw.{ext}.md5",
    output:
        cram=temp("data/source/ggvp/cram/{sample}.raw.{ext}"),
    params:
        col=lambda wildcards: 3 if "crai" in wildcards.ext else 1,
    resources:
        ftp=1,
    benchmark:
        "benchmarks/ggvp_download_cram-{sample}-{ext}.tsv"
    shell:
        r"grep -P '\t{wildcards.sample}\t' {input.man} | awk '{{ print ${params.col} }}' | "
        r"xargs wget --quiet -O {output.cram} -o /dev/null && "
        r"md5sum --status --check {input.md5}"


rule ggvp_filter_iupac_base_codes:
    """
    Remove any reads with IUPAC base codes, as GATK 3.5 cannot handle them

    e.g. SC_GMWOF5428795 and SC_GMJOL5309875 contain `M` and `R` calls

    https://www.bioinformatics.org/sms/iupac.html
    """
    input:
        ref=ancient("data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa"),
        cram="data/source/ggvp/cram/{sample}.raw.cram",
        crai="data/source/ggvp/cram/{sample}.raw.cram.crai",
    output:
        cram=temp("data/source/ggvp/cram/{sample}.filtered.cram"),
        crai=temp("data/source/ggvp/cram/{sample}.filtered.cram.crai"),
    benchmark:
        "benchmarks/ggvp_filter_iupac_base_codes-{sample}.tsv"
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        "samtools view"
        " --expr 'seq=~\"^[ACGTN]+$\"'"
        " --cram"
        " --reference {input.ref}"
        " --write-index"
        " --output {output.cram}"
        " {input.cram}"


rule ggvp_standardise_sample_names:
    """
    Standardise the sample naming in the read groups

    e.g. replace `SC_GMFUL5306388-sc-2012-05-09T14:55:57Z-1371772` with `SC_GMFUL5306388`
    """
    input:
        cram="data/source/ggvp/cram/{sample}.filtered.cram",
        crai="data/source/ggvp/cram/{sample}.filtered.cram.crai",
    output:
        cram="data/source/ggvp/cram/{sample}.cram",
        crai="data/source/ggvp/cram/{sample}.cram.crai",
    benchmark:
        "benchmarks/ggvp_standardise_sample_names-{sample}.tsv"
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        "samtools reheader --command \"sed 's/SM:[^\\t]*/SM:{wildcards.sample}/g'\" {input.cram} > {output.cram} && "
        "samtools index {output.cram}"


def ggvp_list_all_cram():
    samples = pd.read_table(config["source"]["ggvp"]["samples"])

    files = [
        [f"data/source/ggvp/cram/{sample}.cram", f"data/source/ggvp/cram/{sample}.cram.crai"]
        for sample in samples["sample"]
        if sample not in MISSING_CRAM
    ]

    return files


rule ggvp_download_all_cram:
    input:
        ggvp_list_all_cram(),
    output:
        touch("data/source/ggvp/cram/download.done"),

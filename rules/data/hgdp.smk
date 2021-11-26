#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch

"""
Rules to download data files for the Human Genome Diversity Project (HGDP)

https://www.internationalgenome.org/data-portal/data-collection/hgdp
"""


wildcard_constraints:
    sample="[\w.]+",


rule hgdp_download_cram:
    """
    Download bwa-mem CRAM files for each fully-public HGDP sample
    """
    input:
        man="data/source/hgdp/links-to-read-alignments.txt",
    output:
        cram="data/source/hgdp/cram/{sample}.raw.cram",
        crai="data/source/hgdp/cram/{sample}.raw.cram.crai",
    resources:
        ebi_ftp=1,
    conda:
        "../../envs/htslib.yaml"
    shell:
        r"grep -P '^{wildcards.sample}\t' {input.man} | awk '{{ print $2 }}' | "
        r"xargs wget --quiet -O {output.cram} -o /dev/null && "
        r"samtools quickcheck {output.cram} && "
        r"samtools index {output.cram}"


rule hgdp_filter_iupac_base_codes:
    """
    Remove any reads with IUPAC base codes, as GATK 3.5 cannot handle them

    e.g. HGDP00726 and HGDP00466 contain `M` and `R` calls

    https://www.bioinformatics.org/sms/iupac.html
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cram="data/source/hgdp/cram/{sample}.raw.cram",
        crai="data/source/hgdp/cram/{sample}.raw.cram.crai",
    output:
        cram="data/source/hgdp/cram/{sample}.cram",
        crai="data/source/hgdp/cram/{sample}.cram.crai",
    conda:
        "../../envs/htslib.yaml"
    shell:
        "samtools view"
        " --expr 'seq=~\"^[ACGTN]+$\"'"
        " --cram"
        " --reference {input.ref}"
        " --write-index"
        " --output {output.cram}"
        " {input.cram}"


def hgdp_list_all_crams():
    samples = pd.read_table(config["source"]["hgdp"]["samples"])

    files = [
        [f"data/source/hgdp/cram/{sample}.cram", f"data/source/hgdp/cram/{sample}.cram.crai"]
        for sample in samples["sample"]
    ]

    return files


rule hgdp_download_all_cram:
    input:
        hgdp_list_all_crams(),
    output:
        touch("data/source/hgdp/cram/download.done"),

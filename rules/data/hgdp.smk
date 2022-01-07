#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch, temp, ancient

"""
Rules to download data files for the Human Genome Diversity Project (HGDP)

https://www.internationalgenome.org/data-portal/data-collection/hgdp
"""


wildcard_constraints:
    sample="[\w-]+",


rule hgdp_download_cram:
    """
    Download bwa-mem CRAM files for each fully-public HGDP sample
    """
    input:
        man="data/source/hgdp/links-to-read-alignments.txt",
    output:
        cram=temp("data/source/hgdp/cram/{sample}.raw.cram"),
        crai=temp("data/source/hgdp/cram/{sample}.raw.cram.crai"),
    resources:
        ftp=1,
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        r"grep -P '^{wildcards.sample}\t' {input.man} | awk '{{ print $2 }}' | "
        r"xargs wget --quiet -O {output.cram} -o /dev/null && "
        r"samtools quickcheck {output.cram} && "
        r"samtools index {output.cram} 2> /dev/null"


rule hgdp_filter_iupac_base_codes:
    """
    Remove any reads with IUPAC base codes, as GATK 3.5 cannot handle them

    e.g. HGDP00726 and HGDP00466 contain `M` and `R` calls

    https://www.bioinformatics.org/sms/iupac.html
    """
    input:
        ref=ancient("data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa"),
        cram="data/source/hgdp/cram/{sample}.raw.cram",
        crai="data/source/hgdp/cram/{sample}.raw.cram.crai",
    output:
        cram=temp("data/source/hgdp/cram/{sample}.filtered.cram"),
        crai=temp("data/source/hgdp/cram/{sample}.filtered.cram.crai"),
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


rule hgdp_standardise_sample_names:
    """
    Standardise the sample naming in the read groups

    e.g. replace `SAMEA3302828` with `HGDP00526`
    """
    input:
        cram="data/source/hgdp/cram/{sample}.filtered.cram",
        crai="data/source/hgdp/cram/{sample}.filtered.cram.crai",
    output:
        cram="data/source/hgdp/cram/{sample}.cram",
        crai="data/source/hgdp/cram/{sample}.cram.crai",
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        "samtools reheader --command \"sed 's/SM:[^\\t]*/SM:{wildcards.sample}/g'\" {input.cram} > {output.cram} && "
        "samtools index {output.cram}"


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


rule hgdp_download_10x_gvcfs:
    """
    Download the 10x Genomics gVCFs phased with linked-reads
    """
    input:
        man="data/source/hgdp/hgdp-10x-gVCFs.tsv",
    output:
        vcf="data/source/hgdp/gVCF/phase10x/{sample}.phase10x.vcf.gz",
        tbi="data/source/hgdp/gVCF/phase10x/{sample}.phase10x.vcf.gz.tbi",
    resources:
        ftp=1,
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        r"grep -P '^{wildcards.sample}\t' {input.man} | awk '{{ print $2 }}' | "
        r"xargs wget --quiet -O {output.vcf} -o /dev/null && "
        r"bcftools index --tbi {output.vcf}"


def hgdp_list_all_10x_gvcfs():
    samples = pd.read_table("data/source/hgdp/hgdp-10x-gVCFs.tsv")

    files = [f"data/source/hgdp/gVCF/phase10x/{sample}.phase10x.vcf.gz" for sample in samples["sample"]]

    return files


rule hgdp_download_all_10x_gvcfs:
    input:
        hgdp_list_all_10x_gvcfs(),
    output:
        touch("data/source/hgdp/gVCF/phase10x/download.done"),
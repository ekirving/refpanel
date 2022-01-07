#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch, ancient, temp

"""
Rules to download 10x Genomics pre-phased gVCFs
"""


rule phase10x_download_gvcfs:
    """
    Download gVCFs phased with 10x Genomics linked-reads
    
    Add a random wait time before starting the download, to prevent too many requests per second on the FTP service.
    """
    input:
        tsv=lambda wildcards: ancient(config["source"][wildcards.source]["prephased"]),
    output:
        vcf=temp("data/source/{source}/gVCF/phase10x/{sample}.phase10x.raw.vcf.gz"),
        tbi=temp("data/source/{source}/gVCF/phase10x/{sample}.phase10x.raw.vcf.gz.tbi"),
    resources:
        ftp=1,
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        r"sleep $(( $RANDOM % 10 + 1))s && "
        r"grep -P '^{wildcards.sample}\t' {input.tsv} | awk '{{ print $2 }}' | "
        r"xargs wget --quiet -O {output.vcf} -o /dev/null && "
        r"bcftools index --tbi {output.vcf}"


rule phase10x_standardise_sample_names:
    """
    Standardise the sample naming in the 10x gVFCs

    e.g. replace `longranger222_wgs_26958_APPG7555879_GRCh38_gatk` with `APPG7555879`
    """
    input:
        vcf="data/source/{source}/gVCF/phase10x/{sample}.phase10x.raw.vcf.gz",
        tbi="data/source/{source}/gVCF/phase10x/{sample}.phase10x.raw.vcf.gz.tbi",
    output:
        vcf="data/source/{source}/gVCF/phase10x/{sample}.phase10x.vcf.gz",
        tbi="data/source/{source}/gVCF/phase10x/{sample}.phase10x.vcf.gz.tbi",
        sample=temp("data/source/{source}/gVCF/phase10x/{sample}.phase10x.sample"),
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        "echo '{wildcards.sample}' > {output.sample} && "
        "bcftools reheader --samples {output.sample} {input.vcf} > {output.vcf} && "
        "bcftools index --tbi {output.vcf}"


def phase10x_list_all_gvcfs(wildcards):
    source = wildcards.source
    samples = pd.read_table(config["source"][source]["prephased"])

    files = [f"data/source/{source}/gVCF/phase10x/{sample}.phase10x.vcf.gz" for sample in samples["sample"]]

    return files


rule phase10x_download_all_gvcfs:
    input:
        phase10x_list_all_gvcfs,
    output:
        touch("data/source/{source}/gVCF/phase10x/download.done"),

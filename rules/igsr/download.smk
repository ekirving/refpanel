#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import ancient, touch

"""
Rules to download data files from the International Genome Sample Resource (IGSR)

http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/ 
"""

FTP_1000G_RAW_CALLS = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls/"


wildcard_constraints:
    ext="([a-z]+\.)*[a-z]+",


rule md5_1000g_nygc:
    """
    Make an md5 checksum file for validating the downloaded data
    """
    input:
        man="data/1000G_NYGC/20190425_NYGC_GATK_raw_calls_manifest.txt",
    output:
        md5="data/1000G_NYGC/gVCF/{population}/{sample}.g.{ext}.md5",
    params:
        rgx=r"/{population}/Sample_{sample}/.+{ext}\t",
        file="data/1000G_NYGC/gVCF/{population}/{sample}.g.{ext}",
    shell:
        """grep -P '{params.rgx}' {input.man} | awk '{{ print $3" {params.file}" }}' > {output.md5}"""


rule download_1000g_nygc_gvcf:
    """
    Download GATK HaplotypeCaller gVCFs for each high-coverage NYGC 1000G sample
    """
    input:
        md5=ancient("data/1000G_NYGC/gVCF/{population}/{sample}.g.{ext}.md5"),
    output:
        vcf="data/1000G_NYGC/gVCF/{population}/{sample}.g.{ext}",
    resources:
        ebi_ftp=1,
    shell:
        "wget --quiet -O {output.vcf} {FTP_1000G_RAW_CALLS}/{wildcards.population}/Sample_{wildcards.sample}/analysis/{wildcards.sample}.haplotypeCalls.er.raw.{wildcards.ext} && "
        "md5sum --status --check {input.md5}"


def download_1000g_nygc_2504_samples_inputs(_):
    files = []
    with open("data/1000G_NYGC/20190425_NYGC_2504_samples.txt") as fin:
        for record in fin:
            population, sample = record.split()
            files.append(f"data/1000G_NYGC/gVCF/{population}/{sample}.g.vcf.gz")
            files.append(f"data/1000G_NYGC/gVCF/{population}/{sample}.g.vcf.gz.tbi")

    return files


rule download_1000g_nygc_2504_samples:
    input:
        download_1000g_nygc_2504_samples_inputs,
    output:
        touch("data/1000G_NYGC/gVCF/download.done"),

#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch

"""
Rules to download data files from the New York Genome Center (NYGC) high-coverage version of the 1000 Genomes Project (TGP)

https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
"""

FTP_TGP_NYGC = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls_updated"


wildcard_constraints:
    sample="[\w-]+",
    ext="vcf.gz(.tbi)?",


rule tgp_nygc_md5:
    """
    Make an md5 checksum file for validating the 1000G data
    """
    input:
        man="data/source/1000g/20211105_NYGC_GATK_raw_calls_updated_manifest.txt",
    output:
        md5="data/source/1000g/gVCF/{sample}.g.{ext}.md5",
    params:
        rgx=r"{sample}.haplotypeCalls.er.raw.{ext}\s",
        file="data/source/1000g/gVCF/{sample}.g.{ext}",
    shell:
        """grep -P '{params.rgx}' {input.man} | awk '{{ print $3" {params.file}" }}' > {output.md5}"""


rule tgp_nygc_download_gvcf:
    """
    Download GATK HaplotypeCaller gVCFs for each high-coverage NYGC 1000G sample
    """
    input:
        md5="data/source/1000g/gVCF/{sample}.g.{ext}.md5",
    output:
        vcf="data/source/1000g/gVCF/{sample}.g.{ext}",
    resources:
        ebi_ftp=1,
    shell:
        "wget --quiet -O {output.vcf} -o /dev/null {FTP_TGP_NYGC}/{wildcards.sample}.haplotypeCalls.er.raw.{wildcards.ext} && "
        "md5sum --status --check {input.md5}"


def tgp_nygc_list_all_gvcf():
    samples = pd.read_table(config["source"]["1000g"]["samples"])

    files = [
        [f"data/source/1000g/gVCF/{sample}.g.vcf.gz", f"data/source/1000g/gVCF/{sample}.g.vcf.gz.tbi"]
        for sample in samples["sample"]
    ]

    return files


rule tgp_nygc_download_all_gvcf:
    input:
        tgp_nygc_list_all_gvcf(),
    output:
        touch("data/source/1000g/gVCF/download.done"),


rule tgp_nygc_plink_ped:
    """
    Make a PLINK pedigree file for the 602 trios.
    
    NB. `HG02567` appears in the pedigree, but was not sequenced as part of the 602 trios in the NYGC callset
    """
    input:
        tsv="data/source/1000g/20130606_g1k_3202_samples_ped_population.txt",
    output:
        ped="data/source/1000g/1000g-trios.ped",
    shell:
        """awk 'NR>1 && $3!=0 && $4!=0' {input.tsv} | grep -Pv 'HG02567' > {output.ped}"""

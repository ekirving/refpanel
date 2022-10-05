#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch, ancient

"""
Rules to download data files from the New York Genome Center (NYGC) high-coverage version of the 1000 Genomes Project (TGP)

https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
"""


wildcard_constraints:
    sample="[\w-]+",
    ext="vcf.gz(.tbi)?",


rule tgp_nygc_cram_md5:
    """
    Make an md5 checksum file for validating the 1000G data
    """
    input:
        man=ancient("data/source/1000g/igsr_30x_GRCh38.tsv"),
    output:
        md5="data/source/1000g/cram/{sample}.cram.md5",
    params:
        file="data/source/1000g/cram/{sample}.cram",
    benchmark:
        "benchmarks/tgp_nygc_cram_md5-{sample}.tsv"
    shell:
        r"""awk -v FS="\t" '$6=="{wildcards.sample}" {{ print $2" {params.file}" }}' {input.man} > {output.md5}"""


rule tgp_nygc_download_cram:
    """
    Download bwa-mem CRAM files for each high-coverage NYGC 1000G sample
    """
    input:
        man=ancient("data/source/1000g/igsr_30x_GRCh38.tsv"),
        md5="data/source/1000g/cram/{sample}.cram.md5",
    output:
        cram="data/source/1000g/cram/{sample}.cram",
        crai="data/source/1000g/cram/{sample}.cram.crai",
    resources:
        ftp=1,
    benchmark:
        "benchmarks/tgp_nygc_download_cram-{sample}.tsv"
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        r"""awk -v FS="\t" '$6=="{wildcards.sample}" {{ print $1 }}' {input.man} | """
        r"""xargs wget --quiet -O {output.cram} -o /dev/null && """
        r"""md5sum --status --check {input.md5} && """
        r"""samtools index {output.cram} 2> /dev/null"""


def tgp_nygc_list_all_cram():
    samples = pd.read_table(config["source"]["1000g"]["samples"])

    files = [
        [f"data/source/1000g/cram/{sample}.cram", f"data/source/1000g/cram/{sample}.cram.crai"]
        for sample in samples["sample"]
    ]

    return files


rule tgp_nygc_download_all_cram:
    input:
        tgp_nygc_list_all_cram(),
    output:
        touch("data/source/1000g/cram/download.done"),


rule tgp_nygc_gvcf_md5:
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
    benchmark:
        "benchmarks/tgp_nygc_gvcf_md5-{sample}-{ext}.tsv"
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
        ftp=1,
    params:
        ftp="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls_updated",
    benchmark:
        "benchmarks/tgp_nygc_download_gvcf-{sample}-{ext}.tsv"
    shell:
        "wget --quiet -O {output.vcf} -o /dev/null {params.ftp}/{wildcards.sample}.haplotypeCalls.er.raw.{wildcards.ext} && "
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

    There are 602 trios in 1000G, but only 576 families (including some quads and multi-generational families)

    NB. `HG02567` appears in the pedigree, but was not sequenced as part of the 602 trios in the NYGC callset

    NB. there are 8 individuals (i.e., HG00656, HG00657, NA19238, NA19660, NA19661, NA19678, NA19679, NA20282) who are 
    present in 2 families each, due to family naming that doesn't take into account some overlapping pedigrees 
    (see https://www.biorxiv.org/content/10.1101/078600v1.full#app-3)

    So we fix the problematic family codes here:
    * CHS2 = SH074 + SH089
    * YRI1 = Y028 + Y117
    * MXL1 = m004 + m009
    * MXL2 = m008 + m011
    * ASW4 = 2467 + 2469
    """
    input:
        tsv="data/source/1000g/20130606_g1k_3202_samples_ped_population.txt",
    output:
        ped="data/source/1000g/1000g-trios.ped",
    benchmark:
        "benchmarks/tgp_nygc_plink_ped.tsv"
    shell:
        "awk 'NR>1 && $3!=0 && $4!=0' {input.tsv} | "
        "sed -E 's/SH074|SH089/CHS2/' | "
        "sed -E 's/Y028|Y117/YRI1/' | "
        "sed -E 's/m004|m009/MXL1/' | "
        "sed -E 's/m008|m011/MXL2/' | "
        "sed -E 's/2467|2469/ASW4/' | "
        "grep -v 'HG02567' > {output.ped}"

#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import ancient


"""
Rules to download VCF files from the 1000G NYGC callset for evaluation and comparison 

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/
"""


rule tgp_nygc_vcf_md5:
    """
    Make an md5 checksum file for validating the 1000G VCF files
    """
    input:
        tsv=ancient("data/source/1000g/20201028_genotyped-manifest.tsv"),
    output:
        md5="data/evaluation/1000g_nygc/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chr}.recalibrated_variants.vcf.gz.md5",
    params:
        file="data/evaluation/1000g_nygc/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chr}.recalibrated_variants.vcf.gz",
    benchmark:
        "benchmarks/tgp_nygc_vcf_md5-{chr}.tsv"
    shell:
        r"""grep -P '{wildcards.chr}\b' {input.tsv} | awk '{{ print $3" {params.file}" }}' > {output.md5}"""


rule tgp_nygc_vcf_download:
    """
    Download the 1000G VCF files
    """
    input:
        tsv=ancient("data/source/1000g/20201028_genotyped-manifest.tsv"),
        md5="data/evaluation/1000g_nygc/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chr}.recalibrated_variants.vcf.gz.md5",
    output:
        vcf="data/evaluation/1000g_nygc/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chr}.recalibrated_variants.vcf.gz",
        tbi="data/evaluation/1000g_nygc/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chr}.recalibrated_variants.vcf.gz.tbi",
    resources:
        ftp=1,
    benchmark:
        "benchmarks/tgp_nygc_vcf_download-{chr}.tsv"
    conda:
        "../../envs/htslib-1.14.yaml"
    shell:
        r"grep -P '{wildcards.chr}\b' {input.tsv} | awk '{{ print $1 }}' | "
        r"xargs wget --quiet -O {output.vcf} -o /dev/null && "
        r"md5sum --status --check {input.md5} && "
        r"bcftools index --tbi {output.vcf}"

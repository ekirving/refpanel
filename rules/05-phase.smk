#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand, touch

global workflow

"""
Rules to perform variant phasing on a joint-called reference panel  

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/1000G_2020Oct26_NYGC_Phasing_README.pdf 
"""

# TODO need a version that does not rely on trios
# TODO build a `scaffold` using the trios from the extended 1000G NYGC dataset
# TODO fetch the 10x Genomics callset for the HGDP project and use that also


rule shapeit4_phase_vcf_trios:
    """
    Phase the joint-callset, using trios
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_annot_filter_mendel.vcf.gz",
        map="data/reference/GRCh38/genetic_maps.b38.tar.gz",
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_annot_mendel_filter_phased.vcf.gz",
        # TODO does shapeit4 make it's own index?
        # tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_annot_mendel_filter_phased.vcf.gz.tbi",
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_annot_mendel_filter_phased.vcf.log",
    threads: max(workflow.cores / 4, 8)
    conda:
        "../envs/shapeit-4.2.2.yaml"
    shell:
        "shapeit4"
        " --thread {threads}"
        " --input {input.vcf}"
        " --map {input.map}"
        " --region {wildcards.chr}"
        " --sequencing"
        " --output {output.vcf}"
        " --log {log}"


rule shapeit4_phase_all_vcfs:
    input:
        expand(
            "data/panel/{panel}/vcf/{panel}_{chr}_vqsr_annot_mendel_filter_phased.vcf.gz",
            chr=config["chroms"],
            allow_missing=True,
        ),
    output:
        touch("data/panel/{panel}/vcf/phase.done"),

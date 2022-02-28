#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import protected, unpack, expand, touch

from scripts.common import MAX_MEM_MB, JAVA_TEMP_DIR, list_sources

"""
Rules to merge all data source level gVCF files for a reference panel
"""


def gatk3_panel_merge_chrom_gvcf_input(wildcards):
    sources = list_sources(config, wildcards.panel)
    chr = wildcards.chr
    return {
        "ref": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "chr": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.{chr}.bed",
        "gvcfs": [f"data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.gz" for source in sources],
        "tbi": [f"data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.gz.tbi" for source in sources],
    }


# noinspection PyUnresolvedReferences
rule gatk3_panel_merge_chrom_gvcf:
    """
    Combine all source level gVCFs into one reference panel gVCF for each chromosome
    """
    input:
        unpack(gatk3_panel_merge_chrom_gvcf_input),
    output:
        vcf=protected("data/panel/{panel}/gVCF/{panel}_{chr}.g.vcf.gz"),
        tbi=protected("data/panel/{panel}/gVCF/{panel}_{chr}.g.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/gVCF/{panel}_{chr}.g.vcf.log",
    params:
        gvcfs=lambda wildcards, input: [f"--variant {gvcf}" for gvcf in input.gvcfs],
    resources:
        mem_mb=min(50 * 1024, MAX_MEM_MB),  # ~6.67%
        tmpdir=JAVA_TEMP_DIR,
    benchmark:
        "benchmarks/gatk3_multisample_chrom_gvcf-{panel}-{chr}.tsv"
    conda:
        # a bug in gatk v3.5 causes excessive memory usage when combining large numbers of samples
        "../envs/gatk-3.8.yaml"
    shell:
        "gatk3"
        " -XX:ConcGCThreads=1"
        " -Xms{resources.mem_mb}m"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " -T CombineGVCFs"
        " -R {input.ref}"
        " -L {input.chr}"
        " {params.gvcfs}"
        " -o {output.vcf} &> {log}"


rule panel_merge_gvcfs:
    """
    Merge all `gVCF` files in a reference panel.
    """
    input:
        expand("data/panel/{panel}/gVCF/{panel}_{chr}.g.vcf.gz", chr=config["chroms"], allow_missing=True),
    output:
        touch("data/panel/{panel}/gVCF/merge.done"),

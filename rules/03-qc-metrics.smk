#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd
from snakemake.io import touch, multiext

from scripts.common import JAVA_MEMORY_MB

global workflow

"""
Rules to apply quality control (QC) metrics for sample alignments
"""


rule picard_collect_multiple_metrics:
    """
    Collect multiple classes of metrics

    https://gatk.broadinstitute.org/hc/en-us/articles/360036485252-CollectMultipleMetrics-Picard-
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cram="data/source/{source}/cram/{sample}.cram",
    output:
        multiext(
            "data/source/{source}/cram/{sample}",
            ".alignment_summary_metrics",
            ".base_distribution_by_cycle.pdf",
            ".base_distribution_by_cycle_metrics",
            ".insert_size_histogram.pdf",
            ".insert_size_metrics",
            ".quality_by_cycle.pdf",
            ".quality_by_cycle_metrics",
            ".quality_distribution.pdf",
            ".quality_distribution_metrics",
        ),
    log:
        log="data/source/{source}/cram/{sample}.multiple_metrics.log",
    params:
        prefix="data/source/{source}/cram/{sample}",
    resources:
        mem_mb=JAVA_MEMORY_MB,
    conda:
        "../envs/picard-2.5.0.yaml"
    shell:
        "picard"
        " -Xmx{resources.mem_mb}m"
        " CollectMultipleMetrics "
        " R={input.ref}"
        " I={input.cram}"
        " O={params.prefix} 2> {log}"


rule picard_collect_wgs_metrics:
    """
    Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.

    https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard-
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cram="data/source/{source}/cram/{sample}.cram",
    output:
        txt="data/source/{source}/cram/{sample}.wgs_metrics",
    log:
        log="data/source/{source}/cram/{sample}.wgs_metrics.log",
    resources:
        mem_mb=JAVA_MEMORY_MB,
    conda:
        "../envs/picard-2.5.0.yaml"
    shell:
        "picard"
        " -Xmx{resources.mem_mb}m"
        " CollectWgsMetrics "
        " R={input.ref}"
        " I={input.cram}"
        " O={output.txt} 2> {log}"


rule samtools_idxstats:
    """
    Check the coverage per contig in the alignments
    """
    input:
        cram="data/source/{source}/cram/{sample}.cram",
    output:
        idx="data/source/{source}/cram/{sample}.cram.idxstats",
    shell:
        "samtools idxstats {input.cram} > {output.idx}"


rule sample_sex_chr_ratio:
    """
    Calculate the X/Y coverage ratio
    """
    input:
        idx="data/source/{source}/cram/{sample}.cram.idxstats",
    output:
        xyr="data/source/{source}/cram/{sample}.cram.xy_ratio",
    shell:
        "awk '"
        ' $1=="chrX" {{ xcov=$3/$2 }}; '
        ' $1=="chrY" {{ ycov=$3/$2 }}; '
        " END {{ print xcov/ycov }}' {input.idx} > {output.xyr}"


def source_list_metrics(wildcards):
    """List all alignment metrics files for the given data source"""
    source = wildcards.source
    samples = pd.read_table(config["source"][source]["samples"])

    return [
        [
            f"data/source/{source}/cram/{sample}.multiple_metrics.log",
            f"data/source/{source}/cram/{sample}.wgs_metrics.log",
            f"data/source/{source}/cram/{sample}.cram.xy_ratio",
        ]
        for sample in samples["sample"]
    ]


rule source_alignment_metrics:
    """
    Calculate alignment metrics for the given data source.
    """
    input:
        source_list_metrics,
    output:
        touch("data/source/{source}/cram/metrics.done"),

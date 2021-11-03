#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

"""
The human reference genome used by the International Genome Sample Resource (IGSR)

https://www.internationalgenome.org/
"""


rule reference_grch38:
    """
    Download the GRCh38 assembly FASTA, and all the associated index files needed to perform the alignments
    """
    output:
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.dict",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.amb",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.pac",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.sa",
        "data/reference/GRCh38/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz",
        "data/reference/GRCh38/other_mapping_resources/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz",
        "data/reference/GRCh38/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz",
    shell:
        "wget "
        " --mirror"
        " --quiet"
        " --no-host-directories"
        " --cut-dirs=5"
        " --directory-prefix=data/reference/GRCh38/"
        " ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/"


rule reference_grch38_autosomes:
    """
    Make a Picard-style .interval_list of the autosomes

    https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists
    """
    input:
        dict="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.dict",
    output:
        list="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla_autosomes.interval_list",
    shell:
        "grep -P '@HD|SN:chr\d+\b' {input.dict} > {output.list}"


rule reference_grch38_hapmap:
    """Download HapMap 3.3"""
    output:
        "data/reference/GRCh38/other_mapping_resources/hapmap_3.3.hg38.vcf.gz",
    shell:
        "wget --quiet -O {output} https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz"


rule reference_grch38_1000G_omni:
    """Download 1000G Omni 2.5"""
    output:
        "data/reference/GRCh38/other_mapping_resources/1000G_omni2.5.hg38.vcf.gz",
    shell:
        "wget --quiet -O {output} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz"


rule reference_grch38_1000G_snps:
    """Download 1000G Phase 1 high confidence SNPs"""
    output:
        "data/reference/GRCh38/other_mapping_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
    shell:
        "wget --quiet -O {output} https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import temp

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
        " -o /dev/null"
        " ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/"


rule reference_grch38_hapmap:
    """Download HapMap 3.3"""
    output:
        vcf="data/reference/GRCh38/other_mapping_resources/hapmap_3.3.hg38.vcf.gz",
        tbi="data/reference/GRCh38/other_mapping_resources/hapmap_3.3.hg38.vcf.gz.tbi",
    shell:
        "wget --quiet -O {output.vcf} -o /dev/null https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz && "
        "wget --quiet -O {output.tbi} -o /dev/null https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi && "
        "gunzip --test {output.vcf}"


rule reference_grch38_1000G_omni25:
    """Download 1000G Omni 2.5"""
    output:
        vcf="data/reference/GRCh38/other_mapping_resources/1000G_omni2.5.hg38.vcf.gz",
        tbi="data/reference/GRCh38/other_mapping_resources/1000G_omni2.5.hg38.vcf.gz.tbi",
    shell:
        "wget --quiet -O {output.vcf} -o /dev/null https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz && "
        "wget --quiet -O {output.tbi} -o /dev/null https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi && "
        "gunzip --test {output.vcf}"


rule reference_grch38_1000G_phase1_snps:
    """Download 1000G Phase 1 high confidence SNPs"""
    output:
        vcf="data/reference/GRCh38/other_mapping_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        tbi="data/reference/GRCh38/other_mapping_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
    shell:
        "wget --quiet -O {output.vcf} -o /dev/null https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz && "
        "wget --quiet -O {output.tbi} -o /dev/null https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi && "
        "gunzip --test {output.vcf}"


rule reference_ncbi_remapper_fix_errors:
    """
    Fix the VCF validation errors that GATK complains about.

    1. INFO tags cannot contain `=`
    2. INFO column cannot contain whitespace
    """
    input:
        vcf="data/reference/GRCh38/other_mapping_resources/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz",
    output:
        vcf="data/reference/GRCh38/other_mapping_resources/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels_fixed.vcf.gz",
        tbi="data/reference/GRCh38/other_mapping_resources/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels_fixed.vcf.gz.tbi",
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "gunzip -c {input.vcf} | "
        " sed 's/POS=POS-1/POS_POS-1/g' | "
        " sed '/^#/!s/ /_/g' | "
        " bgzip > {output.vcf} && "
        "bcftools index --tbi {output.vcf}"


rule reference_grch38_fai_to_bed:
    """Convert a FASTA index into a BED file"""
    input:
        fai="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai",
    output:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.bed",
    shell:
        r"""awk -v FS="\t" -v OFS="\t" '{{print $1 FS "0" FS ($2-1)}}' {input.fai} > {output.bed}"""


rule reference_grch38_male_haploid_bed:
    """
    Make a bed file for the male haploid regions (i.e., chrX (minus PAR1 and PAR2), chrY and chrM)
    """
    input:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.bed",
    output:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.M.1.bed",
        sex=temp("data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.sex.bed"),
    params:
        par1=r"\t".join(["chrX", "10000", "2781479"]),
        par2=r"\t".join(["chrX", "155701382", "156030895"]),
    conda:
        "../envs/bedtools-2.30.yaml"
    shell:
        r"grep -P 'chrX|chrY|chrM' {input.bed} > {output.sex} && "
        r"printf '{params.par1}\n{params.par2}\n' | "
        r" bedtools subtract -a {output.sex} -b stdin > {output.bed}"


rule reference_grch38_male_diploid_bed:
    """
    Make a bed file for the male diploid regions (i.e., drop the haploid regions)
    """
    input:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.bed",
        hap="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.M.1.bed",
    output:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.M.2.bed",
    conda:
        "../envs/bedtools-2.30.yaml"
    shell:
        "bedtools subtract -nonamecheck -a {input.bed} -b {input.hap} > {output.bed}"


rule reference_grch38_female_haploid_bed:
    """
    Make a bed file for the female haploid regions (i.e., chrM)
    """
    input:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.bed",
    output:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.F.1.bed",
    shell:
        "grep chrM {input.bed} > {output.bed}"


rule reference_grch38_female_diploid_bed:
    """
    Make a bed file for female diploid regions (i.e., drop chrY and chrM)
    """
    input:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.bed",
    output:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.F.2.bed",
    shell:
        "grep -vP 'chrY|chrM' {input.bed} > {output.bed}"


rule reference_grch38_chrom_bed:
    """
    Make a bed file for a specific chromosome
    """
    input:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.bed",
    output:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.{chr}.bed",
    shell:
        r"grep -P '^{wildcards.chr}\t' {input.bed} > {output.bed}"


rule reference_grch38_chrom_others_bed:
    """
    Make a bed file for all the non-standard contigs (i.e., random, unknown, alt, decoy, and HLA)
    """
    input:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.bed",
    output:
        bed="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.others.bed",
    shell:
        r"grep -vP '^chr(\d+|X|Y|M)\t' {input.bed} > {output.bed}"


rule reference_grch38_genetic_map:
    """
    Fetch the GRCh38 genetic map
    """
    output:
        map="data/reference/GRCh38/genetic_maps.b38.tar.gz",
    shell:
        "wget --quiet -O {output.map} -o /dev/null https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz && "
        "gunzip --test {output.map}"

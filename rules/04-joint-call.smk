#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import math

from psutil import virtual_memory
from snakemake.io import protected, unpack, temp, expand, touch

from scripts.utils import list_samples, list_sources

"""
Rules to perform joint genotype calling for the IGSR pipeline

https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/1000G_README_2019April10_NYGCjointcalls.pdf 
"""

# GATK / JAVA default settings
GATK_NUM_THREADS = 4
JAVA_MEMORY_MB = 8 * 1024

# maximum available RAM
MAX_MEM_MB = int(virtual_memory().total / 1024 ** 2) - 1024

# the maximum number of samples to merge together at one time with CombineGVCFs (to prevent massive memory usage)
GATK_BATCH_SIZE = 200

# the default `/tmp` partition is too small
GATK_TEMP_DIR = "./tmp/"


wildcard_constraints:
    chr="(chr(\d+|X|Y|M))|(others)",
    type="SNP|INDEL",


def gatk3_batch_sample_chrom_gvcfs_input(wildcards):
    """Split the samples into batches so we don't use too much RAM"""
    source = wildcards.source
    samples = list_samples(config, source)
    start = (int(wildcards.batch) - 1) * GATK_BATCH_SIZE
    stop = start + GATK_BATCH_SIZE

    return {
        "ref": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "chr": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.{chr}.bed",
        "gvcfs": [f"data/source/{source}/gVCF/{sample}.g.vcf.gz" for sample in samples[start:stop]],
        "tbi": [f"data/source/{source}/gVCF/{sample}.g.vcf.gz.tbi" for sample in samples[start:stop]],
    }


# noinspection PyUnresolvedReferences
rule gatk3_batch_sample_chrom_gvcfs:
    """
    Combine all gVCFs in batches, from one chromosome in one datasource into a multi-sample gVCF
    """
    input:
        unpack(gatk3_batch_sample_chrom_gvcfs_input),
    output:
        vcf=temp("data/source/{source}/gVCF/merged/{source}_{chr}_{batch}.g.vcf.gz"),
        tbi=temp("data/source/{source}/gVCF/merged/{source}_{chr}_{batch}.g.vcf.gz.tbi"),
    log:
        log="data/source/{source}/gVCF/merged/{source}_{chr}_{batch}.g.vcf.log",
    params:
        gvcfs=lambda wildcards, input: [f"--variant {gvcf}" for gvcf in input.gvcfs],
    resources:
        mem_mb=min(28 * 1024, MAX_MEM_MB),  # ~3.71%
        tmpdir=GATK_TEMP_DIR,
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
        " -o {output.vcf} 2> {log}"


def gatk3_multisample_chrom_gvcf_input(wildcards):
    """Split the samples into batches"""
    source = wildcards.source
    chr = wildcards.chr
    num_samples = len(list_samples(config, source))
    batches = range(1, math.ceil(num_samples / GATK_BATCH_SIZE) + 1)

    return {
        "ref": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "chr": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.{chr}.bed",
        "gvcfs": [f"data/source/{source}/gVCF/merged/{source}_{chr}_{batch}.g.vcf.gz" for batch in batches],
        "tbi": [f"data/source/{source}/gVCF/merged/{source}_{chr}_{batch}.g.vcf.gz.tbi" for batch in batches],
    }


# noinspection PyUnresolvedReferences
rule gatk3_multisample_chrom_gvcf:
    """
    Combine all gVCFs batches into one multi-sample gVCF for each chromosome
    """
    input:
        unpack(gatk3_multisample_chrom_gvcf_input),
    output:
        vcf=protected("data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.gz"),
        tbi=protected("data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.gz.tbi"),
    log:
        log="data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.log",
    params:
        gvcfs=lambda wildcards, input: [f"--variant {gvcf}" for gvcf in input.gvcfs],
    resources:
        mem_mb=min(72 * 1024, MAX_MEM_MB),  # ~9.53%
        tmpdir=GATK_TEMP_DIR,
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
        " -o {output.vcf} 2> {log}"


rule source_merge_gvcfs:
    """
    Merge all `gVCF` files in a data source.
    """
    input:
        expand("data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.gz", chr=config["chroms"], allow_missing=True),
    output:
        touch("data/source/{source}/gVCF/merge.done"),


def gatk3_genotype_chrom_gvcf_input(wildcards):
    sources = list_sources(config, wildcards.panel)
    chr = wildcards.chr
    return {
        "ref": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "chr": "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.{chr}.bed",
        "gvcfs": [f"data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.gz" for source in sources],
        "tbi": [f"data/source/{source}/gVCF/merged/{source}_{chr}.g.vcf.gz.tbi" for source in sources],
    }


# noinspection PyUnresolvedReferences
rule gatk3_genotype_chrom_gvcf:
    """
    Jointly call genotypes in all samples for a specific chromosome

    NB. GATK does not honour the --num_threads flag and will use all available cores
    """
    input:
        unpack(gatk3_genotype_chrom_gvcf_input),
    output:
        vcf=protected("data/panel/{panel}/vcf/{panel}_{chr}.vcf.gz"),
        tbi=protected("data/panel/{panel}/vcf/{panel}_{chr}.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}.vcf.log",
    params:
        gvcfs=lambda wildcards, input: [f"--variant {gvcf}" for gvcf in input.gvcfs],
    resources:
        mem_mb=min(94 * 1024, MAX_MEM_MB),  # ~12.4%
        tmpdir=GATK_TEMP_DIR,
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -Xms{resources.mem_mb}m"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " -T GenotypeGVCFs"
        " -R {input.ref}"
        " -L {input.chr}"
        " --num_threads {GATK_NUM_THREADS}"
        " --disable_auto_index_creation_and_locking_when_reading_rods"
        " {params.gvcfs}"
        " -o {output.vcf} 2> {log}"


rule gatk3_split_variants:
    """
    Split variants into SNPs and INDELs to perform the variant recalibration in parallel
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_{chr}.vcf.gz",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_{chr}_{type}.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_{chr}_{type}.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_{type}.vcf.log",
    resources:
        mem_mb=JAVA_MEMORY_MB,
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -Xmx{resources.mem_mb}m"
        " -T SelectVariants"
        " -R {input.ref}"
        " --variant {input.vcf}"
        " --selectTypeToInclude {wildcards.type}"
        " --out {output.vcf} 2> {log}"


rule gatk3_variant_recalibrator_snp:
    """
    Variant Quality Score Recalibration (VQSR) to assign FILTER status for SNPs

    https://gatk.broadinstitute.org/hc/en-us/articles/360035531612?id=1259
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_SNP.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_SNP.vcf.gz.tbi",
        hap="data/reference/GRCh38/other_mapping_resources/hapmap_3.3.hg38.vcf.gz",
        omni="data/reference/GRCh38/other_mapping_resources/1000G_omni2.5.hg38.vcf.gz",
        snps="data/reference/GRCh38/other_mapping_resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="data/reference/GRCh38/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz",
    output:
        recal="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_SNP.recal",
        tranche="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_SNP.tranches",
        plot="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_SNP_plots.R",
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_SNP.log",
    threads: GATK_NUM_THREADS
    resources:
        mem_mb=JAVA_MEMORY_MB,
        tmpdir=GATK_TEMP_DIR,
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " -T VariantRecalibrator"
        " -R {input.ref}"
        " --num_threads {threads}"
        " -input {input.vcf}"
        " -mode SNP"
        " -recalFile {output.recal}"
        " -tranchesFile {output.tranche}"
        " -rscriptFile {output.plot}"
        " -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hap}"
        " -resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni}"
        " -resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.snps}"
        " -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp}"
        " -an QD"
        " -an MQ"
        " -an FS"
        " -an MQRankSum"
        " -an ReadPosRankSum"
        " -an SOR"
        " -an DP"
        " -tranche 100.0"
        " -tranche 99.8"
        " -tranche 99.6"
        " -tranche 99.4"
        " -tranche 99.2"
        " -tranche 99.0"
        " -tranche 95.0"
        " -tranche 90.0 2> {log}"


rule gatk3_variant_recalibrator_indel:
    """
    Variant Quality Score Recalibration (VQSR) to assign FILTER status for INDELs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_INDEL.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_INDEL.vcf.gz.tbi",
        mills="data/reference/GRCh38/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz",
        dbsnp="data/reference/GRCh38/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz",
    output:
        recal="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_INDEL.recal",
        tranche="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_INDEL.tranches",
        plot="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_INDEL_plots.R",
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_INDEL.log",
    threads: GATK_NUM_THREADS
    resources:
        mem_mb=JAVA_MEMORY_MB,
        tmpdir=GATK_TEMP_DIR,
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " -T VariantRecalibrator"
        " -R {input.ref}"
        " --num_threads {threads}"
        " -input {input.vcf}"
        " -mode INDEL"
        " -recalFile {output.recal}"
        " -tranchesFile {output.tranche}"
        " -rscriptFile {output.plot}"
        " -resource:mills,known=true,training=true,truth=true,prior=12.0 {input.mills}"
        " -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp}"
        " -an QD"
        " -an FS"
        " -an ReadPosRankSum"
        " -an MQRankSum"
        " -an SOR"
        " -an DP"
        " -tranche 100.0"
        " -tranche 99.0"
        " -tranche 95.0"
        " -tranche 92.0"
        " -tranche 90.0"
        " --maxGaussians 4 2> {log}"


rule gatk3_apply_recalibration_snp:
    """
    Apply the VQSR to the SNPs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_SNP.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_SNP.vcf.gz.tbi",
        recal="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_SNP.recal",
        tranche="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_SNP.tranches",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_SNP.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_SNP.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_SNP.vcf.log",
    threads: GATK_NUM_THREADS
    resources:
        mem_mb=JAVA_MEMORY_MB,
        tmpdir=GATK_TEMP_DIR,
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " -T ApplyRecalibration"
        " -R {input.ref}"
        " --num_threads {threads}"
        " -input {input.vcf}"
        " -mode SNP"
        " --ts_filter_level 99.80"
        " -recalFile {input.recal}"
        " -tranchesFile {input.tranche}"
        " -o {output.vcf} 2> {log}"


rule gatk3_apply_recalibration_indel:
    """
    Apply the VQSR to the INDELs
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_INDEL.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_INDEL.vcf.gz.tbi",
        recal="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_INDEL.recal",
        tranche="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_INDEL.tranches",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_INDEL.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_INDEL.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_INDEL.vcf.log",
    threads: GATK_NUM_THREADS
    resources:
        mem_mb=JAVA_MEMORY_MB,
        tmpdir=GATK_TEMP_DIR,
    conda:
        "../envs/gatk-3.5.yaml"
    shell:
        "gatk3"
        " -Xmx{resources.mem_mb}m"
        " -Djava.io.tmpdir='{resources.tmpdir}'"
        " -T ApplyRecalibration"
        " -R {input.ref}"
        " --num_threads {threads}"
        " -input {input.vcf}"
        " -mode INDEL"
        " --ts_filter_level 99.0"
        " -recalFile {input.recal}"
        " -tranchesFile {input.tranche}"
        " -o {output.vcf} 2> {log}"


rule picard_merge_variant_vcfs:
    """
    Merge the SNP and INDEL VCF files
    """
    input:
        snp="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_SNP.vcf.gz",
        indel="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_INDEL.vcf.gz",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr.vcf.log",
    resources:
        mem_mb=JAVA_MEMORY_MB,
    conda:
        "../envs/picard-2.5.0.yaml"
    shell:
        "picard"
        " -Xmx{resources.mem_mb}m"
        " MergeVcfs"
        " INPUT={input.snp}"
        " INPUT={input.indel}"
        " OUTPUT={output.vcf} 2> {log}"


# noinspection PyTypeChecker
rule bcftools_samples_file:
    """
    Make a bcftools samples file for calculating INFO tags at the super population level

    https://m.ensembl.org/Help/Faq?id=532
    """
    input:
        tsv=lambda wildcards: config["panel"][wildcards.panel]["samples"],
    output:
        tsv="data/panel/{panel}/{panel}-superpops.tsv",
    params:
        col1=lambda wildcards, input: open(input.tsv).readline().split().index("sample") + 1,
        col2=lambda wildcards, input: open(input.tsv).readline().split().index("superpopulation") + 1,
    shell:
        r"""awk -v FS="\t" 'NR>1 && ${params.col2}!="" {{ print ${params.col1} FS ${params.col2} }}' {input.tsv} > {output.tsv}"""


rule bcftools_norm:
    """
    Split polyallelic INDELs, left-align and normalise.

    NB. Polyallelic INDELs have separate rsIDs for each allele, whereas polyallelic SNPs do not
    """
    input:
        ref="data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr.vcf.gz.tbi",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm.vcf.log",
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "( bcftools norm"
        "   --fasta-ref {input.ref}"
        "   --atomize"
        "   --multiallelics +snps"
        "   -Oz -o {output.vcf} {input.vcf} && "
        "  bcftools index --tbi {output.vcf} "
        ")2> {log}"


rule bcftools_trio_file:
    """
    Convert a PLINK pedigree file into bcftools format (i.e., mother1,father1,child1)
    """
    input:
        ped=lambda wildcards: config["panel"][wildcards.panel].get("pedigree", "/dev/null"),
    output:
        tsv="data/source/{panel}/{panel}-trios.tsv",
    shell:
        """awk '{{ print $4","$3","$2 }}' {input.ped} > {output.tsv}"""


rule bcftools_annotate:
    """
    Annotate variants with dbSNP and fill all tags.

    NB. `+fill-tags --tags all` does not set all tags! 
    see https://github.com/samtools/bcftools/blob/develop/plugins/fill-tags.c#L404
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm.vcf.gz.tbi",
        super="data/panel/{panel}/{panel}-superpops.tsv",
        trios="data/source/{panel}/{panel}-trios.tsv",
        dbsnp="data/reference/GRCh38/dbsnp/GRCh38.dbSNP155.vcf.gz",
    output:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot.vcf.gz.tbi",
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot.vcf.log",
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "( bcftools annotate -a {input.dbsnp} -c ID -Ou {input.vcf} | "
        "  bcftools +mendelian {input.vcf} --mode a --trio-file {input.trios} -Ou | "
        "  bcftools +fill-tags -Oz -o {output.vcf} -- --tags all,F_MISSING --samples-file {input.super} && "
        "  bcftools index --tbi {output.vcf} "
        ") 2> {log} "


rule bcftools_filter_vcf:
    """
    Filter the VCF to only retain high-quality variants for downstream use

    Defined using the following criteria: 
    1) VQSR PASS;
    2) GT missingness < 5%; 
    3) HWE p-value > 1e-10 in at least one of the super-populations;
    4) Mendelian error rate < 5%
    5) MAC >= 2 (i.e., no singletons)
    """
    input:
        vcf="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot.vcf.gz",
        tbi="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot.vcf.gz.tbi",
        super="data/panel/{panel}/{panel}-superpops.tsv",
        trios="data/source/{panel}/{panel}-trios.tsv",
    output:
        vcf=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz"),
        tbi=temp("data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz.tbi"),
    log:
        log="data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.log",
    params:
        # compose the HWE filter for the super-populations found in this reference panel
        hwe=lambda wildcards, input: " | ".join(
            [f"HWE_{pop}>1e-10" for pop in set(line.split()[1] for line in open(input.super))]
        ),
        # use the count of trios to convert the Mendelian error count into a fraction
        max_merr=lambda wildcards, input: max(math.ceil(len(open(input.trios).readlines()) * 0.05), 1),
    conda:
        "../envs/htslib-1.14.yaml"
    shell:
        "( bcftools view --include 'FILTER=\"PASS\" & F_MISSING<0.05 & ({params.hwe}) & MERR<{params.max_merr} & MAC>=2' -Oz -o {output.vcf} {input.vcf} && "
        "  bcftools index --tbi {output.vcf} "
        ") 2> {log} "


rule panel_joint_call:
    """
    Joint-call a reference panel.
    """
    input:
        expand(
            "data/panel/{panel}/vcf/{panel}_{chr}_vqsr_norm_annot_filter.vcf.gz",
            chr=config["chroms"],
            allow_missing=True,
        ),
    output:
        touch("data/panel/{panel}/vcf/joint-call.done"),

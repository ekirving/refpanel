# refpanel

A fully automated and reproducible pipeline for building large reference panels of jointly-called and phased human
genomes.

This pipeline is designed to produce a callset compatible with the [alignment and SNP calling workflow](
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20190405_NYGC_b38_pipeline_description.pdf)
used by the New York Genome Center (NYGC) in their recent [1000G 30x paper](
https://www.biorxiv.org/content/10.1101/2021.02.06.430068); implemented here
using [snakemake](https://snakemake.readthedocs.io/en/stable/)
and [conda](https://docs.conda.io/projects/conda/en/latest/) for full reproducibility.

## Installation

Download the `refpanel` source code

```bash
git clone git@github.com:ekirving/refpanel.git && cd refpanel
```

This pipeline uses the [conda package manager](https://docs.conda.io/projects/conda/en/latest/index.html) (or the
faster [mamba](https://mamba.readthedocs.io/en/latest/index.html) front-end) to handle installation of all software
dependencies. If you do not already have `conda` or `mamba` installed, then please install one first.

Once `conda` is setup, build and activate a new environment for the `refpanel` pipeline

```bash
conda env create --name refpanel --file environment.yaml
```

```bash
conda activate refpanel
```

## Data sources

This pipeline comes preconfigured to build a joint-callset, called `igsr` (n=4,893), involving all publicly available
samples from:

* 1000 Genomes Project (1000G), 30x NYGC version (n=3,202; 2,504 unrelated + 698 related); \
  [doi:10.1101/2021.02.06.430068](https://doi.org/10.1101/2021.02.06.430068)
* Human Genome Diversity Project (HGDP) (n=929; including overlap with SGDP); \
  [doi:10.1126/science.aay5012](https://doi.org/10.1126/science.aay5012)
* Simons Genome Diversity Project (SGDP) (n=231; excluding overlap with 1000G and HGDP); \
  [doi:10.1038/nature18964](https://doi.org/10.1038/nature18964)
* Gambian Genome Variation Project (GGVP) (n=394; excluding overlap with 1000G); \
  [doi:10.1038/s41467-019-13480-z](https://doi.org/10.1038/s41467-019-13480-z)
* Arabian Peninsula Population Genomic (APPG) study (n=137); \
  [doi:10.1016/j.cell.2021.07.013](https://doi.org/10.1016/j.cell.2021.07.013)

The data from these projects is hosted by the
the [International Genome Sample Resource (IGSR) database ](https://www.internationalgenome.org/)
([doi:10.1093/nar/gkw829](https://doi.org/10.1093/nar/gkw829)) and the [European Nucleotide Archive (ENA)](
https://www.ebi.ac.uk/ena/browser/home).

If you wish to build a customised joint-callset (e.g., with non-public samples), please refer to
the [configuration docs](docs/config.md).

### Downloading data

To ensure all data is processed consistently, `refpanel` downloads `CRAM` files for all samples in 1000G, HGDP, SGDP and
GGVP, as well as `gVCF` files for 1000G and `fastq` files for APPG.

To (optionally) pre-fetch all the data dependencies, run:

```bash
./download_data.sh &
```

:warning: **These files are very large**: Please make sure you have [sufficient disk space](docs/diskspace.md) to store
them!

## Joint-calling pipeline

In brief, `refpanel` produces a jointly-called and phased callset via the following steps:

* [Adapter trimming](rules/02-align.smk) with `fastp` (v.0.23.2)
* [Alignment to `GRCh38`](rules/02-align.smk) with `bwa mem` (v0.7.15)
* [Fix-mate](rules/02-align.smk), [merge](rules/02-align.smk), [sort](rules/02-align.smk),
  and [mark duplicates](rules/02-align.smk) with `picard` (v2.5.0)
* [Base recalibration](rules/02-align.smk) with `gatk BaseRecalibrator` (v3.5)
* [Conversion to `cram`](rules/02-align.smk) with `samtools` (v1.14)
* [Per-sample calling of `gVCFs`](rules/04-call.smk) with `gatk HaplotypeCaller` (
  with [sex-dependent ploidy](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls_updated/README_2021November05_NYGCrawcalls_updated.docx))
* [Merging samples](rules/05-merge-calls.smk) with `gatk CombineGVCFs`
* [Joint-calling of all samples](rules/06-joint-call.smk) with `gatk GenotypeGVCFs`
* [Variant quality score recalibration](rules/06-joint-call.smk) with `gatk VariantRecalibrator`
* [Annotation with dbSNP build 155](rules/06-joint-call.smk) with `bcftools` (v1.14)
* [Hard-filtering of SNPs and INDELs](rules/06-joint-call.smk) with `bcftools`:
    1) VQSR PASS;
    2) GT missingness < 5%;
    3) HWE p-value > 1e-10 in at least one super-population;
    4) Mendelian error rate < 5% (using 602 trios from 1000G);
    5) MAC ≥ 2 (i.e., no singletons)
* [Read-based phasing](rules/07-phase-reads.smk) with `whatshap` (v1.2.1) using:
    * _Illumina_ paired-end reads from all projects;
    * _10x Genomics_ long-read sequencing from HGDP (n=XXX);
* [Pedigree phasing](rules/08-phase-trios.smk) with `whatshap` using:
    * Trios from 1000G (n=602);
* [Statistical phasing](rules/09-phase-stat.smk) with `shapeit4` (v4.2.2)
* [Variant effect prediction](rules/10-ensembl-vep.smk) with `ensembl-vep` (v105.0)

For more information, refer to the [DAG of the rule graph](docs/rulegraph.pdf) or the code itself.

## Running the pipeline

To execute the full pipeline, end-to-end, run:

```bash
./run_pipeline.sh &
```

:warning: **This will take a long time**: Please make sure you run this on a server with as many CPUs, and as much RAM,
as possible (e.g., this pipeline was developed and run on a cluster of nodes with 96 cores and 755Gb of RAM each).

The pipeline can also be broken down [into separate steps](docs/steps.md), for distribution across multiple nodes in a
cluster.
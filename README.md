# refpanel

A fully automated and reproducible pipeline for building large reference panels of jointly-called and phased human
genomes, aligned to `GRCh38`.

This pipeline was inspired by the [alignment and SNP calling workflow](
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20190405_NYGC_b38_pipeline_description.pdf)
used by the New York Genome Center (NYGC) in their recent [1000G 30x paper](
https://www.biorxiv.org/content/10.1101/2021.02.06.430068); implemented here, with improvements,
using [snakemake](https://snakemake.readthedocs.io/en/stable/) and
[conda](https://docs.conda.io/projects/conda/en/latest/) for full reproducibility.

:warning: **This pipeline is in active development and subject to ongoing improvements.**

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

This pipeline comes preconfigured to build a joint-callset, called `refpanel-v1` (n=4,833), involving all publicly 
available samples from:

* 1000 Genomes Project (1000G), 30x NYGC version (n=3,202; 2,504 unrelated + 698 related); \
  [doi:10.1101/2021.02.06.430068](https://doi.org/10.1101/2021.02.06.430068)
* Human Genome Diversity Project (HGDP) (n=816; excluding overlap with SGDP); \
  [doi:10.1126/science.aay5012](https://doi.org/10.1126/science.aay5012)
* Simons Genome Diversity Project (SGDP) (n=256; excluding overlap with 1000G); \
  [doi:10.1038/nature18964](https://doi.org/10.1038/nature18964)
* Gambian Genome Variation Project (GGVP) (n=400; excluding overlap with 1000G); \
  [doi:10.1038/s41467-019-13480-z](https://doi.org/10.1038/s41467-019-13480-z)
* Arabian Peninsula Population Genomic (APPG) study (n=137); \
  [doi:10.1016/j.cell.2021.07.013](https://doi.org/10.1016/j.cell.2021.07.013)

Plus additional public genomes from:
* Schuster et al. (2010) *Nature* (n=3); [doi:10.1038/nature08795](https://doi.org/10.1038/nature08795) 
* Meyer et al. (2012) *Science* (n=8); [doi:10.1126/science.1224344](https://doi.org/10.1126/science.1224344)
* Mondal et al. (2016) *Nature Genetics* (n=60); [doi:10.1038/ng.3621](https://doi.org/10.1038/ng.3621)
* McColl et al. (2018) *Science* (n=2); [doi:10.1126/science.aat3628](https://doi.org/10.1126/science.aat3628)
* Lorente-Galdos et al. (2019) *Genome Biology* (n=9); [doi:10.1186/s13059-019-1684-5](https://doi.org/10.1186/s13059-019-1684-5)

The data from these projects is hosted by the
the [International Genome Sample Resource (IGSR) database ](https://www.internationalgenome.org/)
([doi:10.1093/nar/gkw829](https://doi.org/10.1093/nar/gkw829)) and the [European Nucleotide Archive (ENA)](
https://www.ebi.ac.uk/ena/browser/home) ([doi:10.1093/nar/gkq967](https://doi.org/10.1093/nar/gkq967)).

If there is publicly available whole-genome sequencing data that you would like incorporated into `refpanel-v2` please
raise an issue on GitHub with the details of the publication.

If you wish to build a customised joint-callset (e.g., including non-public samples), please refer to
the [configuration docs](docs/config.md).

### Downloading data

To ensure all data is processed consistently, `refpanel` downloads `gVCF` files for 1000G; `CRAM` files for 1000G, HGDP, 
SGDP and GGVP;  and `fastq` files for all other data sources.

To (optionally) pre-fetch all the data dependencies, run:

```bash
./refpanel download_data &
```

All output will be automatically written to a log file `refpanel-<YYYY-MM-DD-HHMM.SS>.log` 

:warning: **These files are very large**: Please make sure you have [sufficient disk space](docs/diskspace.md) to store
them!

## Ancestry composition

Superpopulation assignments are based on the original 1000G, HGDP and SGDP metadata.

| Superpopulation                     | Code | Samples |
|-------------------------------------|------|--------:|
| African Ancestry                    | AFR  |   1,431 |
| American Ancestry                   | AMR  |     564 |
| Central Asian and Siberian Ancestry | CAS  |      25 |
| Central and South Asian Ancestry    | CSA  |     199 |
| East Asian Ancestry                 | EAS  |     826 |
| European Ancestry                   | EUR  |     788 |
| Middle Eastern Ancestry             | MEA  |     297 |
| Oceanian Ancestry                   | OCE  |      38 |
| South Asian Ancestry                | SAS  |     618 |
| West Eurasian Ancestry              | WEA  |      47 |

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
    4) Mendelian error rate < 5%, using trios from 1000G (n=602) and GGVP (n=133);
    5) MAC ≥ 2 (i.e., no singletons)
* [Read-based phasing](rules/07-phase-reads.smk) with `whatshap` (v1.2.1) using:
    * _Illumina_ paired-end reads from all projects;
    * _10x Genomics_ linked-read sequencing from HGDP (n=26) and APPG (n=137);
* [Pedigree phasing](rules/08-phase-trios.smk) with `whatshap` using:
    * Trios from 1000G (n=602) and GGVP (n=133);
* [Statistical phasing](rules/09-phase-stat.smk) with `shapeit4` (v4.2.2)
* [Variant effect prediction](rules/10-ensembl-vep.smk) with `ensembl-vep` (v105.0)

For more information, refer to the [DAG of the rule graph](docs/rulegraph.pdf) or the code itself.

## Running the pipeline

To execute the full pipeline, end-to-end, run:

```bash
./refpanel &
```

All output will be automatically written to a log file `refpanel-<YYYY-MM-DD-HHMM.SS>.log`

:warning: **This will take a long time**: Please make sure you run this on a server with as many CPUs, and as much RAM,
as possible (e.g., this pipeline was developed and run on a cluster of nodes with 96 cores and 755Gb of RAM each).

The pipeline can also be broken down [into separate steps](docs/steps.md), for distribution across multiple nodes in a
cluster.
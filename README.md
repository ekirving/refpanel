# refpanel

A fully automated and reproducible pipeline for building large reference panels of jointly-called and phased human
genomes, aligned to `GRCh38`.

This pipeline was inspired by the [alignment and SNP calling workflow](
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20190405_NYGC_b38_pipeline_description.pdf)
used by the New York Genome Center (NYGC) in their recent paper [High-coverage whole-genome sequencing of the expanded 
1000 Genomes Project cohort including 602 trios](
https://doi.org/10.1016/j.cell.2022.08.004); implemented here, with improvements,
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

This pipeline comes preconfigured to build a joint-callset, called `refpanel-v2` (n=5,100), involving all publicly 
available samples from:

* **1000 Genomes Project (1000G)**, 30x NYGC version (n=3,202; inc. 2,504 phase 3 samples + 698 additional related samples); \
  The 1000 Genomes Project Consortium. (2015) *Nature* [doi:10.1038/nature15393](https://doi.org/10.1038/nature15393) \
  Byrska-Bishop et al. (2022) *Cell* [doi:10.1016/j.cell.2022.08.004](https://doi.org/10.1016/j.cell.2022.08.004)
* **Simons Genome Diversity Project (SGDP)** (n=256; excluding overlap with 1000G); \
  Mallick et al. (2016) *Nature* [doi:10.1038/nature18964](https://doi.org/10.1038/nature18964)
* **Gambian Genome Variation Project (GGVP)** (n=400; excluding overlap with 1000G); \
  Band et al. (2019) *Nature Communications* [doi:10.1038/s41467-019-13480-z](https://doi.org/10.1038/s41467-019-13480-z)
* **Human Genome Diversity Project (HGDP)** (n=816; excluding overlap with SGDP); \
  Bergström et al. (2020) *Science* [doi:10.1126/science.aay5012](https://doi.org/10.1126/science.aay5012) 
* **Arabian Peninsula Population Genomic Study (APPG)**  (n=137); \
  Almarri et al. (2021) *Cell* [doi:10.1016/j.cell.2021.07.013](https://doi.org/10.1016/j.cell.2021.07.013)

Plus additional public genomes from:
* Meyer et al. (2012) *Science* (n=8); [doi:10.1126/science.1224344](https://doi.org/10.1126/science.1224344)
* Mondal et al. (2016) *Nature Genetics* (n=60); [doi:10.1038/ng.3621](https://doi.org/10.1038/ng.3621)
* Rodriguez-Flores et al. (2016) *Genome Research* (n=108); [doi:10.1101/gr.191478.115](https://doi.org/10.1101/gr.191478.115)
* de Barros Damgaard et al. (2018) *Science* (n=41); [10.1126/science.aar7711](https://doi.org/10.1126/science.aar7711)
* McColl et al. (2018) *Science* (n=2); [doi:10.1126/science.aat3628](https://doi.org/10.1126/science.aat3628)
* Lindo et al. (2018) *Science Advances* (n=25); [10.1126/sciadv.aau4921](https://doi.org/10.1126/sciadv.aau4921)
* Gelabert et al. (2019) *BMC Genomics* (n=12); [doi:10.1186/s12864-019-5529-0](https://doi.org/10.1186/s12864-019-5529-0)
* Serra-Vidal et al. (2019) *Current Biology* (n=21); [doi:10.1016/j.cub.2019.09.050](https://doi.org/10.1016/j.cub.2019.09.050)
* Lorente-Galdos et al. (2019) *Genome Biology* (n=9); [doi:10.1186/s13059-019-1684-5](https://doi.org/10.1186/s13059-019-1684-5)
* Crooks et al. (2020) *BMC Genetics* (n=3); [doi:10.1186/s12863-020-00917-4](https://doi.org/10.1186/s12863-020-00917-4)

The data from these projects is hosted by the
[International Genome Sample Resource (IGSR) database ](https://www.internationalgenome.org/)
([doi:10.1093/nar/gkw829](https://doi.org/10.1093/nar/gkw829)) and the [European Nucleotide Archive (ENA)](
https://www.ebi.ac.uk/ena/browser/home) ([doi:10.1093/nar/gkq967](https://doi.org/10.1093/nar/gkq967)).

If there are publicly available whole-genome sequencing data that you would like incorporated into `refpanel-v3` please
raise an issue on GitHub with the details of the publication and they will be considered for inclusion in future releases.

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
| African Ancestry                    | AFR  |   1,460 |
| American Ancestry                   | AMR  |     589 |
| Central Asian and Siberian Ancestry | CAS  |      66 |
| Central and South Asian Ancestry    | CSA  |     199 |
| East Asian Ancestry                 | EAS  |     826 |
| European Ancestry                   | EUR  |     790 |
| Middle Eastern Ancestry             | MEA  |     407 |
| Oceanian Ancestry                   | OCE  |      38 |
| South Asian Ancestry                | SAS  |     678 |
| West Eurasian Ancestry              | WEA  |      47 |

## Joint-calling pipeline

In brief, `refpanel` produces a jointly-called and phased callset via the following steps:

* [Adapter trimming](rules/02-align.smk) with `fastp` (v.0.23.2)
* [Alignment to `GRCh38`](rules/02-align.smk) with `bwa mem` (v0.7.15)
* [Fix-mate](rules/02-align.smk), [merge](rules/02-align.smk), [sort](rules/02-align.smk),
  and [mark duplicates](rules/02-align.smk) with `picard` (v2.5.0)
* [Base recalibration](rules/02-align.smk) with `gatk BaseRecalibrator` (v3.5)
* [Conversion to `cram`](rules/02-align.smk) with `samtools` (v1.14)
* [Per-sample calling of `gVCFs`](rules/04-call.smk) with `gatk HaplotypeCaller`
  (with [sex-dependent ploidy](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls_updated/README_2021November05_NYGCrawcalls_updated.docx))
* [Merging samples](rules/05-merge-samples.smk) with `gatk CombineGVCFs`
* [Joint-calling of all samples](rules/06-joint-call.smk) with `gatk GenotypeGVCFs`
* [Variant quality score recalibration](rules/06-joint-call.smk) with `gatk VariantRecalibrator`
* [Annotation with dbSNP build 155](rules/06-joint-call.smk) with `bcftools` (v1.14)
* [Hard-filtering of SNPs and INDELs](rules/06-joint-call.smk) with `bcftools`:
    1) VQSR PASS;
    2) GT missingness < 5%;
    3) HWE p-value > 1e-10 in at least one super-population;
    4) Mendelian error rate < 5%, using trios from 1000G (n=602) and GGVP (n=133);
    5) MAC ≥ 2 (i.e., no singletons); 
* [Read-based phasing](rules/07-phase-reads.smk) with `whatshap` (v1.2.1) using:
    * _Illumina_ paired-end reads from all projects;
    * _10x Genomics_ linked-read sequencing from HGDP (n=26) and APPG (n=137);
* [Pedigree phasing](rules/08-phase-trios.smk) with `whatshap` using:
    * Trios from 1000G (n=602) and GGVP (n=130);
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
as possible (e.g., this pipeline was developed and run on a cluster of nodes, each with 96 cores and 755Gb of RAM each).

The pipeline can also be broken down [into separate steps](docs/steps.md), for distribution across multiple nodes in a
cluster.
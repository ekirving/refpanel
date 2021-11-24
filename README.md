# refpanel
A fully automated and reproducible pipeline for building large reference panels of jointly-called and phased human 
genomes.

This pipeline is designed to produce a callset compatible with the [alignment and SNP calling workflow](
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20190405_NYGC_b38_pipeline_description.pdf) 
used by the New York Genome Center (NYGC) in their recent [1000G 30x paper](
https://www.biorxiv.org/content/10.1101/2021.02.06.430068); implemented here using [snakemake](https://snakemake.readthedocs.io/en/stable/)
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

This pipeline comes preconfigured to build a joint-callset, called `igsr` (n=4,756), involving all publicly available samples from:
* [1000 Genomes Project (1000G), 30x NYGC version](https://doi.org/10.1101/2021.02.06.430068) (n=3,202; inc. 2,504 unrelated + 698 trios)
* [Human Genome Diversity Project (HGDP) ](https://doi.org/10.1126/science.aay5012) (n=929; including overlap with SGDP)
* [Simons Genome Diversity Project (SGDP)](https://doi.org/10.1038/nature18964) (n=231; excluding overlap with 1000G and HGDP)
* [Gambian Genome Variation Project (GGVP)](https://doi.org/10.1038/s41467-019-13480-z) (n=394; excluding overlap with 1000G)

The data from these projects is hosted by the the [International Genome Sample Resource (IGSR)](
https://doi.org/10.1093/nar/gkw829) database at https://www.internationalgenome.org/  

If you wish to build a customised joint-callset (e.g., with additional samples), you will need to provide a minimal 
amount of [additional metadata](docs/config.md).

### Downloading data

For 1000G and HGDP, we download compatible `gVCF` files output by GATK `HaplotypeCaller`. For SGDP and GGVP,
we download `CRAM` files, as `gVCF` files are not publicly available.

To (optionally) pre-fetch all the data dependencies, run:
```bash
./download_data.sh &
```

:warning: **These files are very large**: Please make sure you have [sufficient disk space](docs/diskspace.md) to store them!

## Joint-calling pipeline

In brief, `refpanel` produces a phased joint-callset by:
* [Alignment](rules/02-align.smk#L25) to GRCh38 with `bwa mem` (v0.7.15)
* [Fix-mate](rules/02-align.smk#L57), [merge](rules/02-align.smk#L91), [sort](rules/02-align.smk#L119), and [mark duplicates](rules/02-align.smk#L146) with `picard` (v2.5.0)
* [Base recalibration](rules/02-align.smk#L176) with `gatk BaseRecalibrator` (v3.5)
* [Conversion to `cram`](rules/02-align.smk#L248) with `samtools` (v1.3.1)
* Per-sample calling of `gVCFs` with `gatk HaplotypeCaller` (with [sex-dependent ploidy](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls_updated/README_2021November05_NYGCrawcalls_updated.docx))
* Joint-calling of all samples with `gatk GenotypeGVCFs`
* Variant quality score recalibration ([VQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612?id=1259)) with `gatk VariantRecalibrator`
* Hard-filtering of SNPs and INDELs with `bcftools` (v1.14) on:
  1) VQSR PASS;
  2) GT missingness < 5%; 
  3) HWE p-value > 1e-10 in at least one super-population;
  4) Mendelian error rate < 5% (based on 1000G trios);
* Phasing with `shapeit4` (v4.2.2) using:
  * Trio data from 1000G; and
  * 10x Genomics long-reads from HGDP;

For more information, see the [DAG of the rule graph](docs/rulegraph.pdf) or refer to the relevant `snakemake` ruleset.

## Running the pipeline

To execute the entire pipeline, end-to-end, run:
```bash
./run_pipeline.sh &
```

:warning: **This will take a long time**: Please make sure you run this on a server with as many CPUs, and as much RAM, 
as possible (e.g., this pipeline was developed on a machine with 96 cores and 755Gb of RAM).

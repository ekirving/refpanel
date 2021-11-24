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

## Configuration

This pipeline comes preconfigured to build a joint-callset, called `igsr`, involving all publicly available samples from:
* [1000 Genomes Project (1000G), 30x NYGC version](https://doi.org/10.1101/2021.02.06.430068) (n=3,202; inc. 2,504 unrelated + 698 trios)
* [Human Genome Diversity Project (HGDP) ](https://doi.org/10.1126/science.aay5012) (n=929; including overlap with SGDP)
* [Simons Genome Diversity Project (SGDP)](https://doi.org/10.1038/nature18964) (n=231; excluding overlap with 1000G and HGDP)
* [Gambian Genome Variation Project (GGVP)](https://doi.org/10.1038/s41467-019-13480-z) (n=394; excluding overlap with 1000G)

The data from these projects is hosted by the the [International Genome Sample Resource (IGSR)](
https://doi.org/10.1093/nar/gkw829) database at https://www.internationalgenome.org/  

If you wish to build a customised joint-callset (e.g., one that contains additional samples), `refpanel` 
requires a [minimal amount of configuration](docs/config.md).

## Downloading data

For 1000G and HGDP, `refpanel` will download compatible `gVCF` files created by GATK `HaplotypeCaller`. For SGDP and GGVP,
`refpanel` will download `CRAM` files, and run `HaplotypeCaller` itself, before passing all samples to GATK `GenotypeGVCFs`

To (optionally) pre-fetch all the data dependencies, run:
```bash
./download_data.sh &
```

:warning: **These files are very large**: Please make sure you have sufficient space to store them!

| Folder                  | Size |
|-------------------------|------|
| data/source/1000g/gVCF/ | 20T  |
| data/source/hgdp/gVCF/  | 1.4T |
| data/source/sgdp/cram/  | 16T  |
| data/source/ggvp/cram/  | 4.8T |


You will also need a lot of additional space to store the intermediate outputs from running the pipeline.

## Running the pipeline

To execute the entire pipeline, end-to-end, run:
```bash
./run_pipeline.sh &
```

:warning: **This will take a long time**: Please make sure you run this on a server with at least 500Gb of RAM and as 
many cores as possible.


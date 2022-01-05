#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2021, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(argparser))
quiet(library(readr))
quiet(library(dplyr))
quiet(library(tidyr))

# get the command line arguments
p <- arg_parser("Convert ENA formattted metadata into `refpanel` format")
p <- add_argument(p, "--ena", help = "ENA metadata")
p <- add_argument(p, "--output", help = "Output file")

argv <- parse_args(p)

# parse the ENA metadata
meta <- read_tsv(argv$ena, na = c("", "NA", "unspecified"), col_types = cols()) %>%

  # split the delimited cols
  separate(col = fastq_md5, into = c("fastq_r1_md5", "fastq_r2_md5"), sep = ";") %>%
  separate(col = fastq_ftp, into = c("fastq_r1_ftp", "fastq_r2_ftp"), sep = ";") %>%

  # select and rename the ENA columns
  select(
    project = study_accession,
    sample = sample_title,
    alias = sample_alias,
    accession = run_accession,
    center = center_name,
    instrument = instrument_platform,
    library = library_name,
    layout = library_layout,
    description = study_title,
    fastq_r1_md5,
    fastq_r2_md5,
    fastq_r1_ftp,
    fastq_r2_ftp
  ) %>%

  # if no library code is defined, use the sample name
  mutate(library = coalesce(library, sample)) %>%

  # add the local file paths
  mutate(fastq_r1 = ifelse(fastq_r1_ftp != "", paste0("data/source/", project, "/fastq/", accession, "_r1.fastq.gz"))) %>%
  mutate(fastq_r2 = ifelse(fastq_r2_ftp != "", paste0("data/source/", project, "/fastq/", accession, "_r2.fastq.gz")))

write_tsv(meta, argv$output)

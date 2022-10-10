#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2021, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(argparser))
quiet(library(dplyr))
quiet(library(readr))
quiet(library(stringr))
quiet(library(tidyr))

# get the command line arguments
p <- arg_parser("Convert ENA formatted metadata into `refpanel` format")
p <- add_argument(p, "--source", help = "Data source name")
p <- add_argument(p, "--ena", help = "ENA metadata")
p <- add_argument(p, "--output", help = "Output file")

argv <- parse_args(p)

# argv$source <- "ggvp"
# argv$ena <- "filereport_read_run_ggvp.tsv"
# argv$output <- "data/source/ggvp/ggvp-accessions.tsv"

# parse the ENA metadata
meta <- read_tsv(argv$ena, na = c("", "NA", "unspecified"), col_types = cols(.default = "c")) %>%
  # select and rename the ENA columns
  select(
    project = study_accession,
    sample = sample_title,
    alias = sample_alias,
    accession = run_accession,
    center = center_name,
    platform = instrument_platform,
    instrument = instrument_model,
    library = library_name,
    layout = library_layout,
    description = study_title,
    fastq_md5,
    fastq_ftp,
  ) %>%
  # drop the alias if it's identical to the sample code
  mutate(alias = ifelse(sample == alias, NA, alias)) %>%
  # if no library code is defined, use the sample name
  mutate(library = coalesce(library, sample))

# flip the sample and alias if the sample code contains whitespace
if (sum(str_count(meta$sample, " "), na.rm = TRUE) > 0 && sum(str_count(meta$alias, " "), na.rm = TRUE) == 0) {
  meta <- rename(meta, alias = sample, sample = alias) %>% relocate(sample, .before = alias)
}

# handle missing sample names
meta <- mutate(meta, sample = coalesce(sample, alias))

# handle single-end libraries
meta_se <- filter(meta, layout == "SINGLE") %>%
  rename(fastq_se_md5 = fastq_md5, fastq_se_ftp = fastq_ftp)

# handle regular paired-end libraries
meta_pe <- filter(meta, layout == "PAIRED", str_count(fastq_ftp, ";") == 1) %>%
  # split the delimited cols
  separate(col = fastq_md5, into = c("fastq_r1_md5", "fastq_r2_md5"), sep = ";") %>%
  separate(col = fastq_ftp, into = c("fastq_r1_ftp", "fastq_r2_ftp"), sep = ";")

# handle special case of paired-end libraries that have unpaired mates
meta_pe_se <- meta %>%
  filter(layout == "PAIRED", str_count(fastq_ftp, ";") == 2) %>%
  # split the delimited cols
  separate(col = fastq_md5, into = c("fastq_se_md5", "fastq_r1_md5", "fastq_r2_md5"), sep = ";") %>%
  separate(col = fastq_ftp, into = c("fastq_se_ftp", "fastq_r1_ftp", "fastq_r2_ftp"), sep = ";")

# join the libraries back together
meta <- bind_rows(meta_se, meta_pe, meta_pe_se) %>%
  arrange(accession) %>%
  # add the explicit ftp:// protocol
  mutate(fastq_se_ftp = ifelse(!is.na(fastq_se_ftp), paste0("ftp://", fastq_se_ftp), NA)) %>%
  mutate(fastq_r1_ftp = ifelse(!is.na(fastq_r1_ftp), paste0("ftp://", fastq_r1_ftp), NA)) %>%
  mutate(fastq_r2_ftp = ifelse(!is.na(fastq_r2_ftp), paste0("ftp://", fastq_r2_ftp), NA)) %>%
  # add the local file paths
  mutate(fastq_se = ifelse(!is.na(fastq_se_ftp), paste0("data/source/", argv$source, "/fastq/", accession, "_se.fastq.gz"), NA)) %>%
  mutate(fastq_r1 = ifelse(!is.na(fastq_r1_ftp), paste0("data/source/", argv$source, "/fastq/", accession, "_r1.fastq.gz"), NA)) %>%
  mutate(fastq_r2 = ifelse(!is.na(fastq_r2_ftp), paste0("data/source/", argv$source, "/fastq/", accession, "_r2.fastq.gz"), NA)) %>%
  # drop all columns that are completely empty
  select_if(~ !all(is.na(.)))

# save the metadata
write_tsv(meta, argv$output, na = "")

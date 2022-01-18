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
p <- arg_parser("Convert IGSR formatted metadata into `refpanel` format")
p <- add_argument(p, "--source", help = "Data source name")
p <- add_argument(p, "--igsr", help = "IGSA metadata")
p <- add_argument(p, "--output", help = "Output file")

argv <- parse_args(p)

# parse the IGSR metadata
meta <- read_tsv(argv$igsr, col_types = cols(.default = "c")) %>%
  # add the data source name
  mutate(source = argv$source) %>%
  # select and rename the IGSR columns
  select(
    source,
    sample = "Sample name",
    sex = "Sex",
    population_name = "Population name",
    superpopulation_name = "Superpopulation name"
  ) %>%
  # some samples have more than one pop and super-pop name (because the sample is part of multiple data sources)
  separate_rows(c(population_name, superpopulation_name), sep = ",") %>%
  arrange(sample)

# drop duplicates from other data sources
if (argv$source == "1000g" || argv$source == "ggvp") {
  meta <- meta %>%
    filter(!str_detect(superpopulation_name, "HGDP|SGDP"))
} else if (argv$source == "hgdp") {
  meta <- meta %>%
    filter(str_detect(superpopulation_name, "HGDP"))
} else if (argv$source == "sgdp") {
  meta <- meta %>%
    filter(str_detect(superpopulation_name, "SGDP"))
}

# make sure that sample codes are unique
duplicates <- meta %>%
  group_by(sample) %>%
  tally() %>%
  filter(n > 1) %>%
  pull(sample)

stopifnot(length(duplicates) == 0)

# load the IGSR populations
pops <- read_tsv("data/source/igsr_populations.tsv", show_col_types = FALSE)

# join the population metadata
meta <- meta %>%
  select(-superpopulation_name) %>%
  left_join(pops, by = "population_name") %>%
  relocate(population, .before = population_name)

# save the metadata
write_tsv(meta, argv$output, na = "")

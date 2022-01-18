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
quiet(library(tidyr))
quiet(library(yaml))

# get the command line arguments
p <- arg_parser("Merge all the data sources into a single reference panel")
p <- add_argument(p, "--private", flag = TRUE, help = "Include private data sources")
p <- add_argument(p, "--output", help = "Output file")

argv <- parse_args(p)

config <- read_yaml("config.yaml")

# load all data sources
sources <- lapply(names(config$source), function(source) {
  if (is.null(config$source[[source]]$private) || config$source[[source]]$private == argv$private) {
    read_tsv(config$source[[source]]$samples, show_col_types = FALSE) %>%
      mutate(source = source) %>%
      select(source, sample, population, population_name, superpopulation, superpopulation_name)
  }
})

panel <- bind_rows(sources)

# make sure that sample codes are unique
duplicates <- panel %>%
  group_by(sample) %>%
  tally() %>%
  filter(n > 1) %>%
  pull(sample)

stopifnot(length(duplicates) == 0)

# make sure that population codes are consistent
duplicates <- panel %>%
  select(population, population_name) %>%
  unique() %>%
  group_by(population_name) %>%
  tally() %>%
  filter(n > 1) %>%
  pull(population_name)

stopifnot(length(duplicates) == 0)


# save the metadata
write_tsv(panel, argv$output, na = "")

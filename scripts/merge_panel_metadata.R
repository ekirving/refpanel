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
p <- add_argument(p, "--output", help = "Output file", default = "data/panel/refpanel-v1/refpanel-v1.tsv")

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
dup <- panel %>%
  group_by(sample) %>%
  tally() %>%
  filter(n > 1)

stopifnot(nrow(dup) == 0)

# make sure that population codes are consistent
dup <- panel %>%
  select(population, population_name) %>%
  unique() %>%
  group_by(population) %>%
  tally() %>%
  filter(n > 1)

panel %>% filter(population %in% dup$population)

stopifnot(nrow(dup) == 0)

# check the super-population pairings are consistent
dup <- panel %>%
  select(population_name, superpopulation_name) %>%
  unique() %>%
  group_by(population_name) %>%
  tally() %>%
  filter(n > 1)

stopifnot(nrow(dup) == 0)

# save the metadata
write_tsv(panel, argv$output, na = "")

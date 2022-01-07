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
sources <- lapply(config$source, function(source) {
  if (is.null(source$private) || source$private == argv$private) {
    read_tsv(source$panel, show_col_types = FALSE)
  }
})

panel <- bind_rows(sources)

# save the metadata
write_tsv(panel, argv$output, na = "")
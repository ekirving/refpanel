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
p <- add_argument(p, "--output", help = "Output file", default = "data/panel/refpanel-v2/refpanel-v2.tsv")
p <- add_argument(p, "--pedigree", help = "Pedigree file", default = "data/panel/refpanel-v2/refpanel-v2.ped")

argv <- parse_args(p)

# argv$private <- TRUE
# argv$output <- "data/panel/cgg-afr-v2/cgg-afr-v2.tsv"

config <- read_yaml("config.yaml")

# load all data sources
sources <- lapply(names(config$source), function(source) {
    if (is.null(config$source[[source]]$private) || config$source[[source]]$private == argv$private) {
        read_tsv(config$source[[source]]$samples, show_col_types = FALSE) %>%
            mutate(source = source) %>%
            select(source, sample, population, population_name, superpopulation, superpopulation_name) %>%
            arrange(sample)
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

stopifnot(nrow(dup) == 0)

# check the super-population pairings are consistent
dup <- panel %>%
    select(population_name, superpopulation_name) %>%
    unique() %>%
    group_by(population_name) %>%
    tally() %>%
    filter(n > 1)

stopifnot(nrow(dup) == 0)

# summarize the metadata
panel %>% nrow()
panel %>%
    group_by(source) %>%
    tally()
panel %>%
    group_by(superpopulation, superpopulation_name) %>%
    tally() %>%
    select(Superpopulation = superpopulation_name, Code = superpopulation, Samples = n)

# load all the pedigrees
pedigrees <- lapply(names(config$source), function(source) {
    if (is.null(config$source[[source]]$private) || config$source[[source]]$private == argv$private) {
        if (!is.null(config$source[[source]]$pedigree)) {
            read_table(config$source[[source]]$pedigree, col_types = cols(), col_names = c("family", "child", "father", "mother", "sex", "population", "superpopulation"))
        }
    }
})

pedigrees <- bind_rows(pedigrees)

# check that each sample only appears in one family
dup <- pedigrees %>%
    pivot_longer(cols = c("child", "father", "mother"), names_to = "relation", values_to = "sample") %>%
    group_by(sample) %>%
    summarise(n = n_distinct(family)) %>%
    filter(n > 1)

stopifnot(nrow(dup) == 0)

# save the metadata
write_tsv(panel, argv$output, na = "")
write_delim(pedigrees, argv$pedigree, col_names = FALSE, na = "")

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
quiet(library(ggplot2))
quiet(library(readr))
quiet(library(sf))
quiet(library(stringr))
quiet(library(tidyr))
quiet(library(purrr))

# get the command line arguments
p <- arg_parser("Convert IGSR formatted metadata into `refpanel` format")
p <- add_argument(p, "--igsr", help = "IGSA metadata", default = "igsr_populations.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "data/source/igsr_populations.tsv")

argv <- parse_args(p)

meta <- read_tsv(argv$igsr, col_types = cols(.default = "c")) %>%
    # select and rename the IGSR columns
    select(
        population = "Population code",
        population_name = "Population name",
        population_desc = "Population description",
        superpopulation = "Superpopulation code",
        superpopulation_name = "Superpopulation name",
        latitude = "Population latitude",
        longitude = "Population longitude"
    )

# map the super-population names to their 3-letter codes
super_codes <- list(
    "African Ancestry" = "AFR",
    "American Ancestry" = "AMR",
    "Central Asian and Siberian Ancestry" = "CAS",
    "Central and South Asian Ancestry" = "CSA",
    "East Asian Ancestry" = "EAS",
    "European Ancestry" = "EUR",
    "Middle Eastern Ancestry" = "MEA",
    "Oceanian Ancestry" = "OCE",
    "South Asian Ancestry" = "SAS",
    "West Eurasian Ancestry" = "WEA"
)

# map the HGDP and SGDP super-populations onto the 1000G names
super_map <- list(

    # HGDP
    "Africa (HGDP)" = "African Ancestry",
    "America (HGDP)" = "American Ancestry",
    "Central South Asia (HGDP)" = "Central and South Asian Ancestry",
    "East Asia (HGDP)" = "East Asian Ancestry",
    "Europe (HGDP)" = "European Ancestry",
    "Middle East (HGDP)" = "Middle Eastern Ancestry",
    "Oceania (HGDP)" = "Oceanian Ancestry",

    # SGDP
    "Africa (SGDP)" = "African Ancestry",
    "America (SGDP)" = "American Ancestry",
    "Central Asia and Siberia (SGDP)" = "Central Asian and Siberian Ancestry",
    "East Asia (SGDP)" = "East Asian Ancestry",
    "Oceania (SGDP)" = "Oceanian Ancestry",
    "South Asia (SGDP)" = "South Asian Ancestry",
    "West Eurasia (SGDP)" = "West Eurasian Ancestry"
)

# create 3-letter population codes based on the population names
new_codes <- abbreviate(sort(unique(meta$population_name)), minlength = 3, strict = TRUE) %>% toupper()

# old_codes[old_codes == "ESN"]
# new_codes[new_codes == "ESN"]

# manual override of some duplicate codes
new_codes["Aleut"] <- "ALE"
new_codes["Brahui"] <- "BRA"
new_codes["Eskimo Naukan"] <- "ESK"
new_codes["Mongola"] <- "MGA"
new_codes["Palestinian"] <- "PAL"

# get the codes that have already been defined
old_codes <- meta %>%
    filter(!is.na(population)) %>%
    select(population, population_name) %>%
    arrange(population_name)

old_codes <- setNames(as.list(old_codes$population), old_codes$population_name)

# keep the codes that have already been defined (e.g. YRI for Yoruba)
for (key in names(old_codes)) {
    new_codes[key] <- old_codes[key]
}

# codes must be unique
stopifnot(length(new_codes) == length(unlist(unique(new_codes))))

# new_codes[duplicated(new_codes)]

# apply the population codes
meta <- meta %>% mutate(population = as.character(new_codes[population_name]))

# check the population codes are unique
dup <- meta %>%
    select(population, population_name) %>%
    unique() %>%
    group_by(population) %>%
    tally() %>%
    filter(n > 1)

stopifnot(nrow(dup) == 0)

# check the pairings are consistent
dup <- meta %>%
    select(population, population_name) %>%
    unique() %>%
    group_by(population_name) %>%
    tally() %>%
    filter(n > 1)

stopifnot(nrow(dup) == 0)

# apply the super-population name and code mapping
meta <- meta %>%
    mutate(superpopulation_name = ifelse(superpopulation_name %in% names(super_map), as.character(super_map[superpopulation_name]), superpopulation_name)) %>%
    mutate(superpopulation = ifelse(superpopulation_name %in% names(super_codes), as.character(super_codes[superpopulation_name]), superpopulation))

# some populations exist in multiple sources, so prefer 1000G -> HGDP -> SGDP for the metadata
prefer_source <- function(pops) {
    unlist(lapply(pops, function(pop) {
        1 + str_detect(pop, "HGDP") * 2 + str_detect(pop, "SGDP") * 3
    }))
}

# only keep one entry per population
meta <- meta %>%
    group_by(population_name) %>%
    slice_min(order_by = prefer_source(population_desc)) %>%
    mutate(population_desc = str_replace(population_desc, " [(](SGDP|HGDP)[)]", ""))

# check the super-population pairings are consistent
dup <- meta %>%
    select(population_name, superpopulation_name) %>%
    unique() %>%
    group_by(population_name) %>%
    tally() %>%
    filter(n > 1)

stopifnot(nrow(dup) == 0)

# save the metadata
write_tsv(meta, argv$output, na = "")

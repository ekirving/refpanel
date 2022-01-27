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
quiet(library(purrr))

# get the command line arguments
p <- arg_parser("Convert EGDP formatted metadata into `refpanel` format")
p <- add_argument(p, "--egdp", help = "EGDP metadata", default = "data/source/egdp/41586_2016_BFnature19792_MOESM240_ESM.txt ")
p <- add_argument(p, "--publication", help = "EGDP publication", default = "this study")
p <- add_argument(p, "--output", help = "Output file", default = "data/source/egdp/egdp-samples.tsv")

argv <- parse_args(p)

igsr <- read_tsv("data/panel/refpanel-v1/refpanel-v1.tsv", col_types = cols(.default = "c"))

meta <- read_tsv(argv$egdp, col_types = cols(.default = "c")) %>%
  mutate(source = "egta", population = "", superpopulation = "") %>%
  unite(col = "population_desc", `Population`, `location`, `Country_of_origin`, remove = FALSE, na.rm = TRUE, sep = ", ") %>%
  # select and rename the EGDP columns
  select(
    source,
    sample = "ID",
    sex = "Sex",
    population,
    population_name = "Population",
    population_desc,
    superpopulation,
    superpopulation_name = "Continent_region_of_origin",
    latitude = "y (rounded for privacy reasons)",
    longitude = "x (rounded for privacy reasons)",
    publication
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

# map the EGDP super-populations on the the IGSR ones
super_map <- list(
  "Africa" = "African Ancestry",
  "Americas" = "American Ancestry",
  "Caucasus" = "West Eurasian Ancestry",
  "Central Asia" = "Central Asian and Siberian Ancestry",
  "East_Asia" = "East Asian Ancestry",
  "Europe" = "European Ancestry",
  "Sahul" = "Oceanian Ancestry",
  "Siberia" = "Central Asian and Siberian Ancestry",
  "South_Asia" = "Central and South Asian Ancestry",
  "Southeast Asia, Island" = "East Asian Ancestry",
  "Southeast Asia, Mainland" = "East Asian Ancestry",
  "West_Asia" = "West Eurasian Ancestry"
)

# remove trailing 's' and any non-alpha characters from the population name before abbreviating
meta <- meta %>%
  mutate(population_name = str_replace(population_name, "ies$", "y")) %>%
  mutate(population_name = str_replace(population_name, "s$", "")) %>%
  mutate(population_name = ifelse(population_name == "Abkhazian", "Abkhasian", population_name)) %>%
  mutate(population_name = ifelse(population_name == "Maasai", "Masai", population_name)) %>%
  mutate(population_name = str_replace_all(population_name, "[^A-Za-z]+", " "))

# create 3-letter population codes based on the population names
new_codes <- abbreviate(sort(unique(meta$population_name)), minlength = 3, strict = TRUE) %>% toupper()

# old_codes[old_codes == "ITL"]
# new_codes[new_codes == "ITL"]

# manual override of some duplicate codes
new_codes["Bakola Pygmy"] <- "BPY"
new_codes["Forest Nenet"] <- "FNN"
new_codes["Croat"] <- "CRO"
new_codes["Eskimo"] <- "EKM"
new_codes["Evenk"] <- "EVK"
new_codes["Italian"] <- "ITA"
new_codes["Kryashen Tatar"] <- "KTT"
new_codes["Mari"] <- "MRI"
new_codes["Shor"] <- "SHO"
new_codes["Sandawe"] <- "SDW"
new_codes["Turkmen"] <- "TKM"

# get the codes that have already been defined
old_codes <- igsr %>%
  filter(!is.na(population)) %>%
  select(population, population_name) %>%
  unique() %>%
  arrange(population_name)

old_codes <- setNames(as.list(old_codes$population), old_codes$population_name)

# keep the codes that have already been defined (e.g. YRI for Yoruba)
for (key in names(old_codes)) {
  new_codes[key] <- old_codes[key]
}

# codes must be unique
stopifnot(length(new_codes) == length(unlist(unique(new_codes))))

# show the duplicates
# new_codes[duplicated(new_codes)]

# apply the population codes
meta <- meta %>% mutate(population = as.character(new_codes[population_name]))

# check the population codes are unique
dup <- bind_rows(meta, igsr) %>%
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

# check the super-population pairings are consistent
dup <- meta %>%
  select(population_name, superpopulation_name) %>%
  unique() %>%
  group_by(population_name) %>%
  tally() %>%
  filter(n > 1)

stopifnot(nrow(dup) == 0)

# filter for this request publication only
meta <- meta %>%
  filter(publication == argv$publication) %>%
  select(-publication)

# save the metadata
write_tsv(meta, argv$output, na = "")

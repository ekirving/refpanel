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
  "Central Asian Ancestry" = "CAS",
  "East Asian Ancestry" = "EAS",
  "European Ancestry" = "EUR",
  "Middle Eastern Ancestry" = "MEA",
  "Oceanian Ancestry" = "OCE",
  "South Asian Ancestry" = "SAS",
  "West Asian Ancestry" = "WAS"
)

# map the HGDP and SGDP super-populations onto the 1000G names
super_map <- list(

  # HGDP
  "Africa (HGDP)" = "African Ancestry",
  "America (HGDP)" = "American Ancestry",
  "Central South Asia (HGDP)" = "South Asian Ancestry",
  "East Asia (HGDP)" = "East Asian Ancestry",
  "Europe (HGDP)" = "European Ancestry",
  "Middle East (HGDP)" = "Middle Eastern Ancestry",
  "Oceania (HGDP)" = "Oceanian Ancestry",

  # SGDP
  "Africa (SGDP)" = "African Ancestry",
  "America (SGDP)" = "American Ancestry",
  "Central Asia and Siberia (SGDP)" = "Central Asian Ancestry",
  "East Asia (SGDP)" = "East Asian Ancestry",
  "Oceania (SGDP)" = "Oceanian Ancestry",
  "South Asia (SGDP)" = "South Asian Ancestry",
  "West Eurasia (SGDP)" = "European Ancestry"
)


# manual override of some super population categories (mostly SGDP)
pop_map <- list(
  # Middle Eastern Ancestry
  "Druze" = "Middle Eastern Ancestry",
  "Mozabite" = "Middle Eastern Ancestry",
  "Palestinian" = "Middle Eastern Ancestry",
  "Iraqi Jew" = "Middle Eastern Ancestry",
  "Bedouin B" = "Middle Eastern Ancestry",
  "Druze" = "Middle Eastern Ancestry",
  "Iranian" = "Middle Eastern Ancestry",
  "Palestinian" = "Middle Eastern Ancestry",
  "Yemenite Jew" = "Middle Eastern Ancestry",
  "Jordanian" = "Middle Eastern Ancestry",
  "Samaritan" = "Middle Eastern Ancestry",

  # Central Asian Ancestry
  "Uygur" = "Central Asian Ancestry",
  "Yakut" = "Central Asian Ancestry",
  "Tajikistan" = "Central Asian Ancestry",

  # West Asian Ancestry
  "Armenian" = "West Asian Ancestry",
  "Turkish" = "West Asian Ancestry",
  "Abkhazia" = "West Asian Ancestry",
  "Georgia" = "West Asian Ancestry",

  # East Asian Ancestry
  "Igorot" = "East Asian Ancestry",
  "Dusun" = "East Asian Ancestry"
)

# create 3-letter population codes based on the population names
new_codes <- abbreviate(sort(unique(meta$population_name)), minlength = 3, strict = TRUE) %>% toupper()

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
  filter(n > 1) %>%
  pull(population)

stopifnot(length(dup) == 0)

# check the pairings are consistent
dup <- meta %>%
  select(population, population_name) %>%
  unique() %>%
  group_by(population_name) %>%
  tally() %>%
  filter(n > 1) %>%
  pull(population_name)

stopifnot(length(dup) == 0)

# manual override of some super population categories (mostly SGDP)
meta <- meta %>%
  mutate(superpopulation_name = ifelse(population_name %in% names(pop_map), as.character(pop_map[population_name]), superpopulation_name))

# apply the super-population name and code mapping
meta <- meta %>%
  mutate(superpopulation_name = ifelse(superpopulation_name %in% names(super_map), as.character(super_map[superpopulation_name]), superpopulation_name)) %>%
  mutate(superpopulation = ifelse(superpopulation_name %in% names(super_codes), as.character(super_codes[superpopulation_name]), superpopulation))

# check the pairings are consistent
dup <- meta %>%
  select(population_name, superpopulation_name) %>%
  unique() %>%
  group_by(population_name) %>%
  tally() %>%
  filter(n > 1) %>%
  pull(population_name)

stopifnot(length(dup) == 0)

# convert to spatial object
meta_sf <- meta %>%
  select(-population_desc) %>%
  unique() %>%
  st_as_sf(coords = c("longitude", "latitude"))

# get a world map
world <- map_data("world")

# colour mapping
cols <- c(
  "African Ancestry" = "#e41a1c",
  "American Ancestry" = "#377eb8",
  "Central Asian Ancestry" = "#4daf4a",
  "East Asian Ancestry" = "#984ea3",
  "European Ancestry" = "#ff7f00",
  "Middle Eastern Ancestry" = "#ffff33",
  "Oceanian Ancestry" = "#a65628",
  "South Asian Ancestry" = "#f781bf"
)

ggplot() +
  geom_map(data = world, map = world, mapping = aes(long, lat, map_id = region)) +
  geom_sf_label(data = meta_sf, mapping = aes(color = superpopulation_name, label = population_name), cex = 2) +
  scale_color_manual(values = cols)

# only keep one description per population
meta <- meta %>%
  group_by(population_name) %>%
  slice(1) %>%
  mutate(population_desc = str_replace(population_desc, " [(](SGDP|HGDP)[)]", ""))

# save the metadata
write_tsv(meta, argv$output, na = "")

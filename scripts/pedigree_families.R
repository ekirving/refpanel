#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2021, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

library(tidyverse)

# load the 1000G pedigree file
ped <- read_delim("data/source/1000g/20130606_g1k_3202_samples_ped_population.txt", delim = " ", na = "0")

# find individuals in more than one family
ped %>%
  select(FamilyID, SampleID, FatherID, MotherID) %>%
  pivot_longer(!FamilyID) %>%
  drop_na() %>%
  select(FamilyID, value) %>%
  unique() %>%
  group_by(value) %>%
  summarise(n = n_distinct(FamilyID)) %>%
  filter(n > 1)

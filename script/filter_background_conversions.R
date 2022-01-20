library(data.table)
library(dplyr)
library(magrittr)
library(tidyverse)

seq_error_rate <- 0.001
conversion_rate <- 1/40

### data
args <- commandArgs(trailingOnly = T)
input_data_per_umi <- args[1]
input_data_per_locus <- args[2]
output_data_per_locus <- args[3]

## merged table / data per umi
data_per_umi <-
  fread(file.path(sample_path, input_data_per_umi)) %>%
  as.data.frame
data_per_umi %<>%
  mutate(sample=sample_name) %>% 
  select(-a , -aG)

## tc table / data per locus
data_per_locus <- 
  fread(file.path(sample_path, input_data_per_locus)) %>%
  as.data.frame
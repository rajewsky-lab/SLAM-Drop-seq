library(data.table)
library(dplyr)
library(magrittr)
library(tidyverse)

### define sample path as the output file folder 
sample_path <- paste0(getwd(), "/output")
params <- read.table(file.path(sample_path, "conversion_error_rate_incorporation_rate.txt"), header = T)
# incorporation_rate <- 0.0477
# conversion_error_rate <- 1.96e-05

conversion_error_rate <- params$conversion_error_rate
incorporation_rate <- params$incorporation_rate

### data
args <- commandArgs(trailingOnly = T)
input_data_per_locus <- args[1]
input_data_per_umi <- args[2]

output_data_per_locus <- args[3]
output_data_per_umi <- args[4]


### get sample name 
sample_name <- unlist(strsplit(basename(input_data_per_umi), "_"))[1]

### read data 
data_per_umi <- 
  fread(input_data_per_umi) %>%
  mutate(sample = sample_name) %>%
  as.data.frame 

data_per_locus <-
  fread(input_data_per_locus) %>%
  as.data.frame

snp_file <- 
  readRDS(file.path(sample_path,
                    "snp_conversion_positions_in_control_samples.rds"))

#### filter background conversions identified from control samples
data_per_locus %<>%
    mutate(snp = ifelse(gene_pos %in% snp_file$gene_pos, 'SNP', 'non_SNP')) %>%
    filter(snp == 'non_SNP') %>%
    droplevels %>%
    dplyr::select(-snp)

### correct the number of TC conversion in each molecule after filtering background conversions 
data_per_umi %<>%
  left_join(data_per_locus %>%
              mutate(corrected_tc = ifelse(round(p_post)>0, 1, 0)) %>%
              group_by(cb_umi_gene) %>%
              summarise(corrected_tc = sum(corrected_tc)) %>%
              ungroup)
data_per_umi[is.na(data_per_umi)]=0

#### calculated the observed new molecule fractions
data_per_umi %<>%
  as.data.frame %>%
  separate(col=cb_umi_gene, into=c("cell", "umi", "gene"), sep="_") 

data_per_umi %<>%
  mutate(type_ori = ifelse(round(corrected_tc)>=1, "new", "old")) %>%
  group_by(gene, .drop=FALSE) %>%
  summarise(n_new = sum(type_ori=="new"),
            n=n()) %>%
  ungroup %>%
  mutate(new_frac=n_new/n) %>%
  select(gene, new_frac) %>%
  left_join(data_per_umi)


#### calculate the posterior probability for a molecule to be newly synthesized
data_per_umi %<>%
  select(sample, cell, umi, gene, new_frac, reads_count, length, ss, t, tC, corrected_tc) %>%
  rowwise %>%
  mutate(
    p_post = 1 / (1 + (1-new_frac)*dbinom(round(corrected_tc), t, conversion_error_rate) / (new_frac*dpois(round(corrected_tc), t*incorporation_rate)))) %>%
  ungroup %>%
  as.data.frame

fwrite(data_per_locus,
       file = output_data_per_locus,
       sep = "\t", col.names=TRUE, row.names=FALSE)

cat("Number of umis:", nrow(data_per_umi), "\n")
fwrite(data_per_umi, 
       file = output_data_per_umi, 
       sep = "\t", col.names=TRUE, row.names=FALSE)


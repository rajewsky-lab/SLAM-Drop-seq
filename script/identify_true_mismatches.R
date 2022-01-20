library(data.table)
library(dplyr)
library(magrittr)
library(tidyverse)

# wd <- getwd()
# sample_path <- paste0(wd, "/output/")
seq_error_rate <- 0.001
conversion_rate <- 1/40

### data
args <- commandArgs(trailingOnly = T)
input_data_per_locus <- args[1]
input_data_per_umi <- args[2]
output_data_per_locus <- args[3]

### get sample name 
sample_name <- unlist(strsplit(basename(input_data_per_locus), "_"))[1]
## tc table / data per locus
data_per_locus <- 
  #fread(file.path(sample_path, input_data_per_locus)) %>%
  fread(input_data_per_locus) %>%
  as.data.frame

data_per_locus %<>%
  mutate(cb_umi_gene=sub('_[^_]*$', "", id),
         gene_pos=sub("^[^_]*_[^_]*_", "", id),
         sample = sample_name)
## merged table / data per umi
data_per_umi <-
  #fread(file.path(sample_path, input_data_per_umi)) %>%
  fread(input_data_per_umi) %>%
  as.data.frame
data_per_umi %<>%
  mutate(sample=sample_name) %>% 
  select(-a , -aG)


### function 
### if conversion rate is fixed, use 1/40, if it is not fixed, compute it for each sample
##### The function of the bayes rule
p_post_real_conversion <-
  function(data_locus,
           data_umi,
           conversion_rate,
           error_rate,
           conversion_rate_fix=TRUE){
    if(!conversion_rate_fix){
      conversion_rate <-
        data_umi %>%
        select(cb_umi_gene, t) %>%
        left_join(data_locus %>%
                    rowwise %>%
                    mutate(ratio = tc_cover / ref_cover) %>%
                    filter(ratio > 0.5) %>%
                    group_by(cb_umi_gene) %>%
                    summarise(n=n()) %>%
                    ungroup) %>%
        as.data.frame
      
      conversion_rate[is.na(conversion_rate)] <- 0
      conversion_rate <- sum(conversion_rate$n) / sum(conversion_rate$t)
      cat(conversion_rate, "\n")
    }
    data_locus %<>%
      mutate(
        p_post=1/(1 + ((1-conversion_rate)*dbinom(tc_cover, ref_cover, error_rate))/(conversion_rate*dbinom(tc_cover, ref_cover, 1-error_rate))))
    
    return(data_locus)
  }


### identify true conversions from sequencing errors
data_per_locus <-
  p_post_real_conversion(
    data_locus = data_per_locus,
    data_umi = data_per_umi,
    conversion_rate = conversion_rate,
    error_rate = seq_error_rate,
    conversion_rate_fix = FALSE)

##To be more precise, calculate the posterior probability of each locus to be true conversion given that we see certain number of reads coverage and mismatches.
fwrite(data_per_locus,
       file = output_data_per_locus,
       sep = "\t", col.names=TRUE, row.names=FALSE)
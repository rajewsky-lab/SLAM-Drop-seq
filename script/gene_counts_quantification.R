library(data.table)
library(dplyr)
library(magrittr)
library(tidyverse)

#sample_path <- paste0(getwd(), "/output")
### data
args <- commandArgs(trailingOnly = T)
input_data_per_umi <- args[1]
output_gene_counts <- args[2]

### get sample name 
sample_name <- unlist(strsplit(basename(input_data_per_umi), "_"))[1]

### read data 
data_per_umi <- 
  fread(input_data_per_umi) %>%
  mutate(sample = sample_name) %>%
  as.data.frame 

################
# Data parsing #
################
startTime <- Sys.time()
cat(paste(Sys.time(), ': reading data: ', sep = ''))
raw.counts <-
  data_per_umi %>%
  mutate(cell_id = paste0(sample, "_", cell)) %>%
  mutate(label = ifelse(round(p_post)>0, "labeled", "unlabeled")) %>%
  mutate(type = paste(label, ss, sep = "_")) %>%
  dplyr::select(cell_id, umi, gene, type) %>%
  group_by(cell_id, gene, type) %>%
  summarise(n = n()) %>%
  ungroup %>%
  pivot_wider(names_from = type, values_from = n, values_fill = 0) %>%
  as.data.frame
  
raw.counts %<>%
  right_join(raw.counts %>% expand(cell_id, gene)) %>%
  arrange(cell_id, gene) %>%
  replace(is.na(.), 0) %>%
  as.data.frame
  
#### check the data columns
rna_type <- c("labeled_spliced", "unlabeled_spliced", "labeled_spliced", "unlabeled_spliced")
if(length(rna_type[!rna_type %in% colnames(raw.counts)]) != 0){
  raw.counts['new_col'] <- 0
  colnames(raw.counts)[which(colnames(raw.counts) == 'new_col')] <- rna_type[!rna_type %in% colnames(raw.counts)]
}
  
cat(paste0(round(Sys.time()-startTime, 2), attr(Sys.time()-startTime, 'units')), '\n')

fwrite(raw.counts, 
       file=output_gene_counts, 
       sep = "\t", col.names=TRUE, row.names=FALSE)



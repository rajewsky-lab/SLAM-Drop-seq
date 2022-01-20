##############################
### load libraries
##############################
library(data.table)
library(dplyr)
library(magrittr)
library(tidyverse)

## working directory
sample_path <- paste0(getwd(), "/output")

###############################################
### parameters
##############################################
snp_read_coverage <- 2  ### read coverage >= 2
snp_mismatch_rate <- 0.5 ### mismatche rate > 50% 

####################################################################
### read sample csv file and keep control and long-labeled samples
###################################################################
samples <- read.table("samples.tsv", header = T)
samples <- 
  samples %>%
  filter(labeling %in% c("control", "long-labeled"))
control_samples <- samples$sample[samples$labeling=="control"]
long_labeled_samples <- samples$sample[samples$labeling=="long-labeled"]
############################################################################
### read files:
####1. merged molecule data (merged.csv) 
####2. conversions at each locus (data_per_locus.csv)
###########################################################################
data_per_umi <- list() ### number of TC conversions and Ts in each UMI
data_per_locus <- list()  ### TC mismatch number and read depth at each locus
for (i in samples$sample){
  ## merged table / data per umi
  data_per_umi[[i]] <-
    fread(file.path(sample_path, paste0(i, "_gene_annotated_with_mismatches_SS_merged.csv"))) %>%
    mutate(sample = i) %>%
    as.data.frame
  data_per_umi[[i]] %<>% mutate(sample=i) %>% select(-a , -aG)

  ## tc table / data per locus
  data_per_locus[[i]] <-
    fread(file.path(sample_path, paste0(i, "_data_per_locus.csv"))) %>%
    as.data.frame
  data_per_locus[[i]] %<>%
    mutate(cb_umi_gene=sub('_[^_]*$', "", id),
           gene_pos=sub("^[^_]*_[^_]*_", "", id),
           sample = i)
}

#############################################
#### background conversions
#### generate background T->C conversions (SNP) file
#### if the coverage is >=2 and >50% of the reads are T->C mismatches, it is identified as conversion
#############################################
data_con <- as.data.frame(data.table::rbindlist(data_per_locus[control_samples]))
data_con %<>%
  mutate(corrected_tc_cover = ifelse(round(p_post)>0, tc_cover, 0)) %>%
  group_by(gene_pos) %>%
  summarise(ref_cover=sum(ref_cover),
            tc_cover=sum(tc_cover),
            corrected_tc_cover=sum(corrected_tc_cover)) %>%
  ungroup %>%
  mutate(ratio=corrected_tc_cover/ref_cover) %>%
  as.data.frame

data_con %<>%
  mutate(snp = ifelse(ratio>snp_mismatch_rate & ref_cover>=snp_read_coverage, 'SNP', 'non_SNP'))

cat(length(data_con$snp[data_con$snp=="SNP"]),
    "out of", nrow(data_con),
    "(", length(data_con$snp[data_con$snp=="SNP"])/nrow(data_con), ")",
    "positions with at least one mismatch are identified as SNPs in the control samples.\n")

data_con %<>% filter(snp=="SNP")

#############################################################
## save the background (SNP) conversion data frame as rds file
#############################################################
saveRDS(data_con, file = paste0(sample_path, "/snp_conversion_positions_in_control_samples.rds"))


########################################################
#### remove the background conversions from these samples
#########################################################
for (i in samples$samples){
  data_per_locus[[i]] %<>%
    mutate(snp = ifelse(gene_pos %in% data_con$gene_pos, 'SNP', 'non_SNP')) %>%
    filter(snp == 'non_SNP') %>%
    droplevels %>%
    dplyr::select(-snp)
}

for (i in samples$sample){
  cat(i)
  data_per_umi[[i]] %<>%
    left_join(data_per_locus[[i]] %>%
                mutate(corrected_tc = ifelse(round(p_post)>0, 1, 0)) %>%
                group_by(cb_umi_gene) %>%
                #summarise(corrected_tc = sum(p_post)) %>%
                summarise(corrected_tc = sum(corrected_tc)) %>%
                ungroup)
  data_per_umi[[i]][is.na(data_per_umi[[i]])]=0
}

###########################################################
#### conversion error rate
#### After correcting for background conversions (SNPs), the T->C conversion rate in control samples is calculated as the conversion error rate
#### calculate the ratio of sums
###########################################################
cat("ratio of sums: \n")
conversion_error_rate  <-
  data.table::rbindlist(data_per_umi[control_samples]) %>%
  filter(sample %in% control_samples, ss=="unspliced")

conversion_error_rate <- 
  sum(conversion_error_rate$corrected_tc)/sum(conversion_error_rate$t)
cat("The conversion error rate from control sample is: ", conversion_error_rate, "\n")

### incorporation rate 
cat("ratio of sums: \n")
incorporation_rate <-
  data.table::rbindlist(data_per_umi[long_labeled_samples]) %>%
  filter(sample %in% long_labeled_samples, ss=="unspliced") 

incorporation_rate <- 
  sum(incorporation_rate$corrected_tc) / sum(incorporation_rate$t)
cat("4sU incorporation rate in new RNAs calculated from 4sU 24h labeled sample is:", incorporation_rate, "\n")


###################################################################
## write the conversion error rate and 4sU incorpororation rate into a data frame 
##################################################################
write.table(data.frame(conversion_error_rate=conversion_error_rate,
                       incorporation_rate=incorporation_rate),
            file=paste0(sample_path, "/conversion_error_rate_incorporation_rate.txt"),
            sep = "\t", row.names=FALSE)


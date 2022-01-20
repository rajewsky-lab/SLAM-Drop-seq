library(data.table)
library(dplyr)
library(magrittr)
# agrments : 1. input--dge.summary_ori.txt, 2.output: best_cellbarcodes.txt
args <- commandArgs(trailingOnly = T)

h <- as.data.frame(fread(args[1], header=T, stringsAsFactors = F))
h %<>%  filter(NUM_GENES>=200)
write.table(h[,1], file=args[2], col.names = F, row.names = F, quote = F)
library(ggplot2)
library(data.table)
# agrments : 1. input--cell_read.counts.txt.gz, 2. output: png  3. output: cell_number.txt 4. output: best_cellbarcodes.txt
args <- commandArgs(trailingOnly = T)

h = as.data.frame(fread(args[1], header=F, stringsAsFactors = F))
y = cumsum(h$V1[1:10000])
y = y/max(y)
Y = length(y)
hdf <- data.frame("y" = y, "x" = 1:Y)

histplot <- (ggplot(data = hdf, aes(x = x, y = y)) 
             + geom_line(col = "steelblue", size = 0.3)
             + theme_minimal() + xlab("cell barcodes sorted by number of reads") 
             + scale_y_continuous(limits = c(0, 1), expand = c(0,0))
             + ylab("cumulative fraction of reads")
             + theme(text = element_text(family = "IPAPMincho", size = 5),
                     panel.grid.major.y = element_line(size = 0.3),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title.y = element_text(vjust = 1),
                     axis.title.x = element_text(vjust = 0)))
ggsave(histplot, file = args[2], width = 3, height = 2)

# Calculate the point where the elbow occurs and write it in a file
mindist = 1
dmax = dist(rbind(c(1/Y, y[1]), c(Y/Y, y[Y])), method="euclidean")
for (i in 1:Y) {
  d1 = dist(rbind(c(1/Y, y[1]), c(i/Y, y[i])), method="euclidean")
  d2 = dist(rbind(c(i/Y, y[i]), c(Y/Y, y[Y])), method="euclidean")
  if (d1 == 0 | d2 == 0) {next}
  min_temp = abs((d1^2 + d2^2 - dmax^2)/(2*d1*d2))
  if (min_temp < mindist) {
    mindist = min_temp
    whichmin = i
  }
}

# Write it to a file
write.table(whichmin, file=args[3], col.names = F, row.names = F)

# select vest cell barcode 
cb = data.frame(cb=h[1:whichmin,2])
write.table(cb, file=args[4], col.names = F, row.names = F, quote = F)

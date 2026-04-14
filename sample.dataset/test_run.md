```r

setwd('~/Desktop/scripts/EpiSignal/') 
source("EpiSignalEM.R")
data <- read.table(
  "~/Desktop/toy_H3K4me2_2.tsv",
  header = TRUE,
  sep = "\t"
)

output_EpiSignal <- EpiSignalEM(log(data[,1:2] + 1))

pdf("~/Desktop/scripts/EpiSignal/results/EpiSignal_plot.pdf")
plot.episignal_em(output_EpiSignal)
dev.off()

sink("~/Desktop/scripts/EpiSignal/results/EpiSignal_summary.txt")
summary.episignal_em(output_EpiSignal)
sink()

```

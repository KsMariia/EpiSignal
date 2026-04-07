# EpiSignal

EpiSignal is an implementation of an expectation maximisation (EM) algorithm for splitting genes into signal and noise based on the corresponding matrix obtained by 
RNA-seq, ChIP-seq or ATAC-seq. The model assumes a mixture of Normal distributions and estimates the probability of belonging to signal for every provided gene.

# Dependencies

install.packages(c(
  "BAMMtools",
  "modi",
  "bayestestR",
  "ggplot2",
  "gridExtra",
  "dplyr",
  "tidyr",
  "ggrepel",
  "RColorBrewer"
))

# Running the model

Input data.frame or matrix


output_EpiSig <- EpiSignalEM(log(rna+1))
summary(output_EpiSig) 
plot(output_EpiSig)

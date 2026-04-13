# EpiSignal

EpiSignal is an implementation of an expectation maximisation (EM) algorithm for classifying genes into signal and noise based on the RNA-seq, ChIP-seq or ATAC-seq matrix. The model assumes a mixture of two normal distributions corresponding to signal and noise (or 'on'/'off' states) and estimates the probability of belonging to signal for every provided gene.

# Dependencies

```r

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
```

# Running the model

```

EpiSignalEM(
  data,
  xl = NULL,
  mu_0m = NULL,
  sigma_0m = NULL,
  mu_1m = NULL,
  sigma_1 = NULL,
  conv_thresh = 1e-8
)

```

Input is a data.frame or a matrix of log-transformed counts with genes in rows and samples in columns. Each input parameter of the EpiSignalEM() function is described in the following table:

| Paremeter     | Description                                         |
| ------------- | --------------------------------------------------- |
| `data`        | data frame (genes × samples)                        |
| `xl`          | initial probabilities of belonging to signal        |
| `mu_0m`       | means of noise per sample                           |
| `sigma_0m`    | variances of noise per sample                       |
| `mu_1m`       | means of signal per sample                          |
| `sigma_1`     | variance of signal common for all samples           |
| `conv_thresh` | convergence threshold (default value is `1e-8`)     |

If the starting values are not provided they are initialised internally using Jenk's clustering and the converges of an algorithm is guaranteed to a local maximum. When the starting values are estimated the probabilities of belonging to signal are calculated on the E-step. On the M-step the parameters of the mixture distrubutions are renewed. Algorithm stops when the change of the log-likelihood between iterations is smaller than the conv_thresh parameter.

# Getting the output

```
output_EpiSigal <- EpiSignalEM(log(data+1))
summary.episignal_em(output_EpiSignal) 
plot.episignal_em(output_EpiSignal)
```

An episignal_em object is produced by the EpiSignalEM() function is a list that contains the log-likelihood value on the last iteration, the total number of iterations algorithm required to converge under a specified threshold, a vector of probabilities of belonging to signal for each gene, the final estimations of the model parameters, the filtered input dataset, the 'off' genes and the normalised signal. 

summary.episignal_em() provides a more detailed summary of the run, displaying the number of genes assigned to signal and noise, estimated parameters of signal and noise distributions, characteristics of the log-likelihood and input data characteristics such as KL-divergence and normality checks (Shapiro-Wilk).

plot.episignal_em() creates sets of plots to illustrate the input data structure, the algorithm convergense, the quality of the approximation and the classification into 'on'/'off' states. This includes multiple convergence plots of mean noise/signal parameters as well as the log-likelihood across the iterations. The quality of the approximation is assessed by the comparison of the empirical and fitted distributions for each sample. The normalised signal is represented by the densities, QQ-plots and plots comparing the estimated and actual means and standard deviations.   

# Sample dataset

In a 'sample.dataset' folder you can find a test count matrix for demonstrating purposes, which contains two wild type samples of CD8+ T cells H3K4me2 ChIP-seq in columns used in the original publication and downsampled to 3000 genes, with the corresponding R script for running EpiSignal algorithm in a 'test_run.R' file. The produced pictures and statistics are represented in an output folder.

# Citations

...





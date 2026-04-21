# EpiSignal

EpiSignal is a gene-specific model for signal detection and normalisation within a single cell-type that is applicable to, eg., DNA accessibility data (ATAC-seq), histone mark (ChIP-seq) data, and mRNA-seq gene expression data. Log-transformed counts of mapped reads in the broad promoter regions associated with epigenetic modifications and in the exonic regions of the transcriptome typically show bimodal occupancy patterns across genes. To account for the bimodal distribution, we use an expectation maximisation (EM) algorithm on a Gaussian mixture to optimally discriminate between genes in 'on' and 'off' (or 'signal' and 'noise’) states.
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

# Arguments

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

The obligate input argument of the EpiSignalEM() function is a data.frame or matrix of log-transformed counts with genes in rows and samples in columns. Further, optional arguments are described in the following table:

| Paremeter     | Description                                         |
| ------------- | --------------------------------------------------- |
| `data`        | data frame (genes × samples)                        |
| `xl`          | initial probabilities of belonging to signal        |
| `mu_0m`       | means of noise per sample                           |
| `sigma_0m`    | variances of noise per sample                       |
| `mu_1m`       | means of signal per sample                          |
| `sigma_1`     | variance of signal common for all samples           |
| `conv_thresh` | convergence threshold (default value is `1e-8`)     |

If the starting values are not provided assignments to 'on'/'off' states are initialised using Jenk's clustering. The Expectation-Maximization (EM) algorithm then iterates until the change in the Q-function is smaller than the  conv_thresh parameter, which corresponds to a local optimum.

# Values 

Running the model:

```
output_EpiSigal <- EpiSignalEM(log(data+1))
summary.episignal_em(output_EpiSignal) 
plot.episignal_em(output_EpiSignal)
```

EpiSignalEM() returns a list containing the following components: 

| Field | Description |
|------|-------------|
| `data` | the input data |
| `likelihood` | the final value of the Q-function at the last iteration |
| `k` | the total number of iterations algorithm until convergence under the threshold |
| `xl` | a vector of probabilities of belonging to signal for each gene |
| `estimates` | the final values of the parameter estimates |
| `samples` | all parameter estimates for all the samples |
| `filtered_out` | genes filtered out, because all samples had zero reads |
| `norm_hard_thresh` | normed values of genes assigned to signal for all samples |
| `hard_thresh_index` | assigment to signal |

summary.episignal_em() provides a more detailed summary of the run: 

| Field | Description |
|------|-------------|
| `filtered_out` | Number of Filtered Out Genes (no signal in any sample) |
| `noise_vs_signal` | Number of Noise and Signal Genes |
| `estimated_means_noise` | Estimated Means for Noise Genes per Sample |
| `estimated_sds_noise` | Estimated Sds for Noise Genes per Sample |
| `estimated_means_signal` | Estimated Means for Signal Genes per Sample |
| `estimated_sd_signal` | Estimated Sd for Signal Genes across Samples |
| `input_data_summary` | Input Data Characteristics per Sample: KL-divergences, overlap percentages between signal and noise, and p-values returned by the Shapiro-Wilks test |
| `output_data_summary` | Normalised Data Characteristics per Sample: KL-divergence, Shapiro-Wilks statistic and p-value |
| `log_likelihood | mean, sd and median of the likelihood |

plot.episignal_em() creates sets of plots: 

| Plot Group | Description |
|-----------|-------------|
| Convergence plots | Parameters and likelihood convergence |
| Sample histograms | Raw data distributions |
| Density fits | Data vs fitted mixture |
| QQ plots | Model fit diagnostics |
| Density | Raw data density distribution |
| Parameter evaluation plots | Mean and Sd agreement with estimates |
| Probability plots | Signal probability (xl) distribution |
| Gene summary plot | Expression vs signal probability bins |
| Normalised signal density | Signal vs reference normal |
| Normalised QQ plot | Normality check |
| Normalised stats plot | Mean and SD diagnostics |

# Example dataset

In a 'sample.dataset' folder you can find a test count matrix for demonstrating purposes, which contains two wild type samples of CD8+ T cells H3K4me2 ChIP-seq in columns used in the original publication and downsampled to 3000 genes, with the corresponding R script for running EpiSignal algorithm in a 'test_run.md' file. The produced statistics and pictures are represented in an EpiSignal_summary.txt and EpiSignal_plots.pdf files.

# Citations

...





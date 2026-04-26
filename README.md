# EpiSignal

EpiSignal is a gene-specific model for signal detection and normalisation within a single cell-type that is applicable to, eg., DNA accessibility data (ATAC-seq), histone mark (ChIP-seq) data, and mRNA-seq gene expression data. Log-transformed counts of mapped reads in the broad promoter regions associated with epigenetic modifications and in the exonic regions of the transcriptome typically show bimodal occupancy patterns across genes. To account for the bimodal distribution, we use an expectation maximisation (EM) algorithm on a Gaussian mixture to optimally discriminate between genes in 'on' and 'off' (or 'signal' and 'noise’) states.

<img width="1646" height="1158" alt="Overall_before_after_wLetters (1)-1" src="https://github.com/user-attachments/assets/3296f3ff-3c53-47f9-80a6-2c7394455910" />

On the picture above you can see the densities of log-transformed counts for the samples of mRNA-seq, ATAC-seq, and H3K4me2 ChIP-seq data in CD8+ T cells used in the original publication before (A, B, C, respectively) and after (D, E, F, respectively) EpiSignal filtering and normalisation. Bimodal distribution of the unfiltered data (A-C) corresponds to the two 'on' and 'off' categories, which transform into a single filtered 'signal' peak centered at 0 after filtering and normalisation (D-F).

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

| Paremeter | Description |
|------|-------------|
| `data` | the input data |
| `likelihood` | the final value of the Q-function, ie the expectation of the log-likelihood, at the last iteration |
| `k` | the total number of iterations algorithm until convergence under the threshold |
| `xl` | a vector of probabilities of belonging to signal for each gene |
| `estimates` | the final values of the parameter estimates |
| `samples` | all parameter estimates for all the samples |
| `filtered_out` | genes filtered out, because all samples had zero reads |
| `norm_hard_thresh` | normed values of genes assigned to signal for all samples |
| `hard_thresh_index` | assignment to signal |

summary.episignal_em() provides a more detailed summary of the run: 

| Paremeter | Description |
|------|-------------|
| `filtered_out` | number of filtered out genes (no signal in any sample) |
| `noise_vs_signal` | number of noise and signal genes |
| `estimated_means_noise` | estimated means for noise genes per sample |
| `estimated_sds_noise` | estimated sds for noise genes per sample |
| `estimated_means_signal` | estimated means for signal genes per sample |
| `estimated_sd_signal` | estimated sd for signal genes across samples |
| `input_data_summary` | input data characteristics per sample: KL-divergences, overlap percentages between signal and noise, and p-values returned by the Shapiro-Wilks test |
| `output_data_summary` | normalised data characteristics per sample: KL-divergence, Shapiro-Wilks statistic and p-value |
| `log_likelihood` | mean, sd and median of the likelihood |

plot.episignal_em() creates sets of plots: 

| Plot | Description |
|-----------|-------------|
| Convergence plots | parameters and likelihood convergence |
| Sample histograms | raw data distributions |
| Density fits | data vs fitted mixture |
| QQ plots | model fit diagnostics |
| Density | raw data density distribution |
| Parameter evaluation plots | mean and sd agreement with estimates |
| Probability plots | signal probability (xl) distribution |
| Gene summary plot | expression vs signal probability bins |
| Normalised signal density | signal vs reference normal |
| Normalised QQ plot | normality check |
| Normalised stats plot | mean and sd diagnostics |

# Example dataset

The 'sample.dataset' folder includes: ‘toy_H3K4me2.tsv’ a count matrix containing two wild type samples of CD8+ T cells H3K4me2 ChIP-seq in columns used in the original publication and downsampled to 3000 genes; 'test_run.md' the R script; ‘EpiSignal_summary.txt’ contains the output statistics and ‘EpiSignal_plots.pdf’ the output figures.

# Citations

...





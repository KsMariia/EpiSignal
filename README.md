# EpiSignal

EpiSignal is an implementation of an expectation maximisation (EM) algorithm for splitting genes into signal and noise based on the RNA-seq, ChIP-seq or ATAC-seq matrix. The model assumes a mixture of two Normal distributions corresponding to signal and noise (or 'on'/'off' states) and estimates the probability of belonging to signal for every provided gene.

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

Input is a data.frame or a matrix of log-transformed values with genes in rows and samples in columns.

| Параметр      | Описание                                            |
| ------------- | --------------------------------------------------- |
| `data`        | data frame (genes × samples)                        |
| `xl`          | initial probabilities of belonging to signal        |
| `mu_0m`       | means of noise per sample                           |
| `sigma_0m`    | variances of noise per sample                       |
| `mu_1m`       | means of signal per sample                          |
| `sigma_1`     | variance of signal common for all samples           |
| `conv_thresh` | convergence threshold (default value is `1e-8`)     |

If the starting values are not set they are initialised internally using Jenks clustering and the converges of an algorithm is guaranteed to a local maximum. When the starting values are estimated the probabilities of belonging to signal are calculated on the E-step. On the M-step the parameters of the distrubutions are renewed. Algorithm stops when the change of the log-likelihood is smaller than conv_thresh.

# Getting the output

```
output_EpiSig <- EpiSignalEM(log(data+1))
summary(output_EpiSig) 
plot(output_EpiSig)
```

The EpiSignalEM function provides the number of 'on'/'off' genes, estimated parameters, KL-divergence, overlap distributions, check of normality (Shapiro-Wilk) and log-likelihood. 


# Citations

@misc{episignalem,
  author = {Names},
  title = {Title},
  year = {2026},
  howpublished = {GitHub}
}





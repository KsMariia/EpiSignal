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

| Paremeter     | Description                                         |
| ------------- | --------------------------------------------------- |
| `data`        | data frame (genes × samples)                        |
| `xl`          | initial probabilities of belonging to signal        |
| `mu_0m`       | means of noise per sample                           |
| `sigma_0m`    | variances of noise per sample                       |
| `mu_1m`       | means of signal per sample                          |
| `sigma_1`     | variance of signal common for all samples           |
| `conv_thresh` | convergence threshold (default value is `1e-8`)     |

If the starting values are not provided they are initialised internally using Jenks clustering and the converges of an algorithm is guaranteed to a local maximum. When the starting values are estimated the probabilities of belonging to signal are calculated on the E-step. On the M-step the parameters of the mixture distrubutions are renewed. Algorithm stops when the change of the log-likelihood between iterations is smaller than the conv_thresh parameter.

# Getting the output

```
output_EpiSigal <- EpiSignalEM(log(data+1))
summary.episignal_em(output_EpiSignal) 
plot.episignal_em(output_EpiSignal)
```

The EpiSignalEM() returns an episignal_em object and provides the number of 'on'/'off' genes, estimated parameters of the distributions, KL-divergence, overlap distributions, check of normality (Shapiro-Wilk) and a log-likelihood. 

str(output)
Основные элементы:
Поле	Описание
xl	Вероятности принадлежности к сигналу
hard_thresh_index	Классификация (TRUE = signal)
estimates	Финальные параметры модели
samples	Траектории параметров (по итерациям)
likelihood	Log-likelihood
filtered_out	Гены без сигнала

summary.episignal_em() provides the summary of the run, displaying the number of genes assigned to signal or noise, estimated parameters of signal and noise distributions, input data characteristics such as KL-divergence and normality checks (Shapiro-Wilk) and characteristics of the log-likelihood.

plot.episignal_em() creates sets of plots to illustrate the input data structure, the algorithm convergense, the quality of the approximation and the classification into 'on'/'off' states. This includes multiple convergence plots of mean noise/signal parameters as well as the log-likelihood across the iterations. The quality of the approximation is assessed by the comparison of the empirical and fitted distributions for each sample. The normalised signal is represented by the densities, QQ-plots and comparisons of the estimated and actual means and standard deviations.   


# Citations

@misc{episignalem,
  author = {Your Name},
  title = {EpiSignalEM},
  year = {2026},
  howpublished = {GitHub}
}





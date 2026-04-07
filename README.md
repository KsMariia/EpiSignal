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

EpiSignalEM(
  data,
  xl = NULL,
  mu_0m = NULL,
  sigma_0m = NULL,
  mu_1m = NULL,
  sigma_1 = NULL,
  conv_thresh = 1e-8
)
Input data.frame or matrix with genes in rows and samples in columns.


output_EpiSig <- EpiSignalEM(log(rna+1))
summary(output_EpiSig) 
plot(output_EpiSig)

| Параметр      | Описание                                       |
| ------------- | ---------------------------------------------- |
| `data`        | Матрица экспрессии (гены × образцы)            |
| `xl`          | Начальные вероятности принадлежности к сигналу |
| `mu_0m`       | Средние значения noise по образцам             |
| `sigma_0m`    | SD noise                                       |
| `mu_1m`       | Средние значения signal                        |
| `sigma_1`     | SD signal (одно значение)                      |
| `conv_thresh` | Порог сходимости (по умолчанию `1e-8`)         |


Если параметры не заданы — они инициализируются автоматически (Jenks clustering).**








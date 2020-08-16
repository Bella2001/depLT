# depLT - Estimation of Survival under Dependent Left Truncation

The  `depLT` package is designed for estimation of a survival function from left-truncated data where truncation time and lifetime cannot be assumed independent. The approach assumes that the dependence between the truncation time and the studied lifetime is induced by covariates all of which are measured. The package implements two methods, a case-weights (CW) estimator and a time-varying-weights (TVW) estimator, both based on the  Inverse-Probability-of-Selection-Weighting as described in Vakulenko-Lagun, Qian, Chiou, Wang,  Betensky, *Estimation under covariate-induced dependent truncation using inverse probability weighting* (2020, submitted).

The  `depLT` package can be installed by
```{r}
devtools::install_github("Bella2001/depLT")
```

 
 

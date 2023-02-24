# Randomization inference for the 'Neyman' Null (zero average treatment effect)

Ding and Wu (2020) and Ding and Zhao (2021) regression adjustment
procedures for Randomization Inference on Neyman ('weak' =: Average
Treatment Effect = 0) as opposed to the traditional Fisher ('sharp' =:
Treatment effect = 0 for everyone). See relevant papers or [Ding's slides](https://www.dropbox.com/s/fj1q0fk57v8ctaz/uchicago.pdf?raw=1) for an introduction.

Single function `ri_weaknull` in `R/` directory. 

```r
set.seed(42)
source("R/ri_fn.R")
```
```r
# %% generate data with small effect
df = dgp(τf = \(x) 0.1, n = 1e4)
```

## difference in means

```r
# %% no covariates
lm_robust(y~w, df) %>% broom::tidy() %>% dplyr::filter(term == "w")
# term   estimate std.error statistic p.value
# w      0.1257   0.08216   1.53      0.126
```

## lin regression with robust SE

```r
# %% standard lin regression
lm_lin(y~w, ~ V1+V2+V3+V4+V5+V6, df) %>% broom::tidy() %>% dplyr::filter(term == "w")

# term   estimate std.error statistic p.value
# w       0.114   0.04428   2.573     0.01008
```

## RI for weak null without covariates

permute treatment label and report studentized test statistic (coef / robust SE on treatment indicator).

```r
# %% randomization inference for fisher null - no covariates
ri_weaknull(df, 'y', 'w', xn = NULL, K = 1e4)
# p-value 0.0687
```

## RI for weak null with covariates

```r
# %% randomization inference for fisher null - with covariates
ri_weaknull(df, 'y', 'w', xn = paste0("V", 1:6), K = 1e4, cores = 16)
# p-value 0.0052
```

Corrections / Suggestions welcome.

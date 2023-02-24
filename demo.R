# %% ####################################################
set.seed(42)
source("R/ri_fn.R")
# %% generate data with small effect
df = dgp(τf = \(x) 0.1, n = 1e4)
# %% no covariates
lm_robust(y~w, df) %>% broom::tidy() %>% dplyr::filter(term == "w")
# term	 estimate	std.error	statistic	p.value
# w	     0.1257	  0.08216	  1.53	    0.126
# %% standard lin regression
lm_lin(y~w, ~ V1+V2+V3+V4+V5+V6, df) %>% broom::tidy() %>% dplyr::filter(term == "w")
# term	 estimate	std.error	statistic	p.value
# w       0.114	  0.04428	  2.573	    0.01008
# %% randomization inference for fisher null - no covariates
ri_weaknull(df, 'y', 'w', xn = NULL, K = 1e4)
# p-value 0.0687
# %% randomization inference for fisher null - with covariates
ri_weaknull(df, 'y', 'w', xn = paste0("V", 1:6), K = 1e4, cores = 16)
# p-value 0.0052
# %%

# %%
#' Covariate-Adjusted Randomization inference for the Weak Null for the Average Treatment Effect
#'
#' @description Implements Ding and Zhao (2022) procedure for covariate-adjusted
#' randomization inference for the Neyman/'weak' null of 0 Average treatment
#' effect using studentized test statistic τ/robust_se from Lin regression. This
#' is typically more relevant for applied practice than the Fisher null of 0
#' treatment effect for all units. When no covariates are supplied, falls back
#' to Ding and Wu (2020) procedure for Randomization Inference for the Neyman
#' Null without covariates.
#' @param df data.table with outcome, treatment, and covariates
#' @param yn outcome name (string)
#' @param wn treatment name (string)
#' @param xn covariate names (string) or NULL (default)
#' @param K Number of permutations (integer) - defaults to 1000
#' @param cores Number of cores to parallelize on (if 0, uses replicate instead of mcReplicate)
#' @param fulldist Boolean for whether to return p-value or full distribution (TRUE returns list)
#' @return scalar p-value or list with elements 'nulldist' and 'real' depending on fulldist flag.

#' @references Wu and Ding (2020) Randomization tests for weak null hypotheses in randomized experiments. Journal of the American Statistical Association
#' @references Zhao and Ding (2020) Covariate-adjusted Fisher randomization tests for the average treatment effect. Journal of Econometrics
#' @export
#' @importFrom estimatr lm_robust lm_lin
ri_weaknull = function(df, yn, wn, xn = NULL, K = 1000L, cores = 8,
    fulldist = FALSE){
  setnames(df, c(yn, wn), c('y', 'w'))
  # different studentized stat depending on whether covariates are present or not
  if(is.null(xn)){
    # anonymous function for test statistic
    test_stat = \(df){
      m = estimatr::lm_robust(y ~ w, df)
      res = summary(m)$coefficients[2, 1:2]
      studentized_stat = res[1]/res[2]
    }
  } else {
    test_stat = \(df){
      m = estimatr::lm_lin(y ~ w,
                  as.formula(glue::glue("~{glue_collapse(xn, '+')}")), df)
      res = summary(m)$coefficients[2, 1:2]
      studentized_stat = res[1]/res[2]
    }
  }
  # store realised studentized stat
  realt = test_stat(df)
  # store n1
  n1 = sum(df[[wn]])
  ######################################################################
  # permute treatments, compute studentized statistic
  ######################################################################
  oneRep = \(dat, n1, test_stat){
    # treatment vector
    dat$w = 0
    # assign n1 of them to 1
    dat$w[sample(1:nrow(dat), n1)] = 1
    # compute test stat
    test_stat(dat)
  }
  ######################################################################
  # replicate
  df2 = copy(df)
  if(cores > 1) { # only tested on linux
    nullstats = mcreplicate::mc_replicate(n = K, oneRep(df2, n1, test_stat), mc.cores = cores)
  } else {
    nullstats = replicate(K, oneRep(df2, n1, test_stat))
  }
  # returns - either return full distribution for plotting
  if(fulldist) return(list(nulldist = nullstats, real = realt))
  # or just p-value
  return(mean(nullstats > realt))
}

# %% DGP to generate a synthetic dataset
#' @param n number of observations
#' @param p number of covariates
#' @param τf treatment effect function (may be constant, but still has to be a function)
#' @param y0f baseline PO effect function
#' @param πf propensity score function
#' @return data.table with outcome, treatment, and covariates
dgp = \(n = 1e3,
        p = 6,
        τf  = \(x) 1/3 ,
        y0f = \(x) 4*pmax(x[1] + x[2], 0) + sin(x[5]) * pmax(x[6], 0.5),
        πf  = \(x) 1/2
        ){
    X = matrix(runif(n*p, -2, 2), n, p)
    # generate treatment, heterogeneity, baseline
    W = rbinom(n, 1, plogis(apply(X, 1, πf)))
    τ = apply(X, 1, τf)
    Y0 = apply(X, 1, y0f)
    Y = Y0 + W * τ + rnorm(n)
    data.table(y = Y, w = W, X)
}

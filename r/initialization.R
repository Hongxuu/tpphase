
#' @description Random initialize parameters. Uniform(0, 1) for now to initialize etas, but for betas
#' use mnlogit to initialize since it is important to get roughly good betas to compute llk. According to 
#' Hongxu interaction term associated with position are not that significant compared with other terms, so 
#' exclude them for now, so we have 10 beta under each category
#' 
#' @param n_class number of class (A1 A2 D1 D2 here).
#' @param num_cat number of category (for multinomial logistic regression ATCG here)
#' @param num_beta number of betas
#' @param seed
#' @return parameters for EM

ini <- function(dat, n_observation, formula, n_class = 4, num_cat = 4, seed = 0) {
  par <- list()
  #set.seed(seed)
  #par$eta <- runif(n_class, 0, 1)
  #par$eta <- par$eta/sum(par$eta)
  par$eta <- rep(0.25, 4)
  par$wic <- 0
  par$ins_rate <- 1e-5
  par$del_rate <- 1e-5
  par$excluded_read <- rep(0, n_observation)
  Mpar <- Mstep(dat, par = par, formula = formula, weights = FALSE)
  par$beta <- Mpar$beta
  return(par)
}
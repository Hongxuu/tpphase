##### Try to use -1*(loglik + weights %*% log(baseProbVec)) 
dd <- fromat_data(dat_info = d, haplotype = haps[[9]])
dd$nuc <- as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G")[as.character(dd$nuc)])
dd$hap_nuc <- as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G")[as.character(dd$hap_nuc)])
dd <- dd[, !names(dd) %in% c("id")] # (modified_moligit exclude first column: id)

wic <- t(resu[[9]]$param$w_ic)
weights <- data.table::rbindlist(foreach(i = 1:ncol(wic)) %dopar% {
  data.frame(wic = rep(wic[, i], read_length[i]))
})
weights <- weights$wic

x <- resu[[10]]$param$beta %>% as.vector()

opt_beta <- function(x, dd, weights, change_var = -Inf, change_po) {
  choice.set <- c("A", "C", "G", "T")
  # choice.set <- unique(data[[choiceVar]])
  K <- length(choice.set)
  N <- nrow(data)/K
  p = 10
  choiceVar <- "nuc"
  ind <- sort_ind(dd[[choiceVar]], nrow(dd)) + 1
  dd <- dd[ind[, 1], ]
  response <- dd["mode"]
  response <- response[(N + 1):(K * N), ]
  X <- formDesignMat(dat = dd, N = N)
  if(change_var != -Inf) {
    x[change_po] <- change_var
  }
  probMat <- X %*% matrix(x, nrow = p, ncol = K-1, byrow=FALSE)
  
  loglik <- drop(as.vector(probMat) %*% (weights * response))
  # Convert utility to probabilities - use logit formula
  probMat <- exp(probMat)                           # exp(utility)
  baseProbVec <- 1/(1 + rowSums(probMat))           # P_i0
  probMat <- probMat * matrix(rep(baseProbVec, K-1),
                              nrow = N, ncol = K-1) # P_ik
  # Negative log-likelihood
  -1*(loglik + weights %*% log(baseProbVec))  
}


sapply(X = seq(-0.01, -4, by = 0.1), FUN = opt_beta, x= x, dd = dd, weights = weights,change_po = 1)
opt_beta(x, dd, weights, change_var, change_po)





####### Gradiant
responseMat <- matrix(response, nrow=N, ncol=(K-1))
baseResp <- rep(1, N) - rowSums(responseMat)
xgrad <- as.vector(crossprod(X, weights * (responseMat - probMat)))
-1 * c(xgrad)

library(nloptr)
nloptr(x0 = c(rep(0, 30)), eval_f =  function(par) opt_beta(par, dd, weights), eval_grad_f = 
         function(par) nl.grad(par, opt_beta, dd = dd, weights = weights), lb = rep(-20, 30), ub= rep(20, 30),
       opts = list("algorithm"="NLOPT_GD_STOGO", print_level = 3, maxeval = 10))








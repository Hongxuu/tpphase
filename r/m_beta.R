
#' @description fit mnlogit and update beta
#' @param dat expanded data
#' @param par parameters
#' @param ncores number of ncores to run mnlogit
#' @return betas and CE_llk

Mstep <- function(dat, read_length, par, weight_id, formula, num_cat = 4, ncores, weights = TRUE) {
  #give a initial value for optimization
  if (weights) {
    weights <- data.table::rbindlist(foreach(i = 1:ncol(par$wic)) %dopar% {
      data.frame(wic = rep.int(par$wic[, i], read_length[i]))
    })
    if(!is.null(weight_id))
      weights <- weights[-c(weight_id), ]
    weights <- weights$wic
    start <- par$beta %>% as.vector
    fit <- modified_mnlogit(formula = formula,
                   data = dat, weights = weights, choiceVar = "nuc", ncores = ncores, start = start)
  }
  else {
    fit <- modified_mnlogit(formula = formula,
                   data = dat, choiceVar = "nuc", ncores = ncores)
  }
  
  #A <- fit$coefficients[which(str_detect(attr(fit$coefficients, "names"), ":A") == 1)]
  C <- fit$coefficients[which(str_detect(attr(fit$coefficients, "names"), ":C") == 1)]
  G <- fit$coefficients[which(str_detect(attr(fit$coefficients, "names"), ":G") == 1)]
  T <- fit$coefficients[which(str_detect(attr(fit$coefficients, "names"), ":T") == 1)]
  
  mat <- matrix(cbind(C, G, T), ncol = 3)
  beta <- mat
  logLik <- fit$logLik
  
  res <- list()
  res$logLik <- logLik
  res$beta <- beta
  
  return(res)
}

#' @description prepare data and call mnlogit
#' @return logLik of mnlogit and betas

m_beta <- function(res, d, weight_id, data, id, formula, reads_lengths, ncores) {
  par <- list()
  
  if(length(res$excluded_id) != 0) {
    data_rm <- data %>% filter(mode == 1) 
    weight_id <- c(weight_id, which(data_rm$id %in% res$excluded_id))
    data <- data %>% filter(!(id %in% res$excluded_id))
  }
  
  par$wic <- t(res$param$w_ic)
  par$beta <- res$param$beta #beta from last step as starting value
  data <- data[, !names(data) %in% c("id")] # (modified_moligit exclude first column: id)
  
  Mpar <- Mstep(dat = data, read_length = reads_lengths, weight_id= weight_id, formula = formula, par = par, 
                ncores = ncores, weights = TRUE)
  par$beta <- Mpar$beta 
  par$wic <- res$param$w_ic ## For the use of update haplotype
  par$eta <- res$param$mixture_prop
  par$del_rate <- res$param$del_rate
  par$ins_rate <- res$param$ins_rate
  par$excluded_read <- res$param$excluded_read
  
  results <- list()
  results$par <- par
  results$CEllk <- Mpar$logLik
  
  return(results)
}

###################################################################################################

# read_hap <- function(path, n_class = 4) {
#   init_hap <- readFastq(path)
#   reads <- as(sread(init_hap), "matrix")[1:n_class, ]
#   ncol = ncol(reads)
#   reads_num <- as.numeric(c("A" = "0", "T" = "2", "C" = "1", "G" = "3")[t(reads)])
#   reads_num <- matrix(reads_num, ncol)
#   t(reads_num)
# }

#' @description Compute the likelihood for each position of a read and likelihood for that read
#' @param i ith observation
#' @param read_length a vector to store all of the observation length
#' @param haplotype one haplotype (out of 4)
#' @param dat_short not expanded data
#' @param par a list to store parameters
#' @param n_class number of class (A1 A2 D1 D2 here).
#' @param find_max indicate if using this function to infer TRUE haplotypes based on finding the maxi likelihood
#' @return list store likelihood for each position of a read and likelihood for that read

# llk <- function(i, read_length, haplotype, dat_short, par, n_class = 4, find_max = FALSE) {
#   res <- list()
#   site_llk <- rep(0, read_length[i])
#   for(j in 1:read_length[i]) {
#     ## predictors
#     if(find_max) {
#       if(haplotype == 'C') {
#         hap_nuc <- c(1, 0, 0)
#       } else if(haplotype == 'G') {
#         hap_nuc <- c(0, 1, 0)
#       } else if(haplotype == 'T') {
#         hap_nuc <- c(0, 0, 1)
#       } else if(haplotype == 'A') {
#         hap_nuc <- c(0, 0, 0)
#       }
#     } else {
#       if(haplotype[j] == 'C') {
#         hap_nuc <- c(1, 0, 0)
#       } else if(haplotype[j] == 'G') {
#         hap_nuc <- c(0, 1, 0)
#       } else if(haplotype[j] == 'T') {
#         hap_nuc <- c(0, 0, 1)
#       } else if(haplotype[j] == 'A') {
#         hap_nuc <- c(0, 0, 0)
#       }
#     }
#     subset <- dat_short %>% filter(idx == i)
#     qua <- subset$qua[1 + (j-1) * n_class]
#     read_pos <- subset$read_pos[1 + (j-1) * n_class]
#     ref_pos <- subset$ref_pos[1 + (j-1) * n_class]
#     
#     ##### notice the predictors should be chaged according to choice of number of betas
#     pred <- c(1, read_pos, ref_pos, qua, hap_nuc, qua*hap_nuc)
#     tail <- log(1/(1 + sum(exp(t(pred) %*% par$beta))))
#     
#     if(subset$read_nuc[1 + (j-1) * n_class] == 'C') {
#       xb <- t(pred) %*% par$beta[ ,1]
#     } else if(subset$read_nuc[1 + (j-1) * n_class] == 'G') {
#       xb <- t(pred) %*% par$beta[ ,2]
#     } else if(subset$read_nuc[1 + (j-1) * n_class] == 'T') {
#       xb <- t(pred) %*% par$beta[ ,3]
#     } else if(subset$read_nuc[1 + (j-1) * n_class] == 'A') {
#       xb <- 0
#     }
#     site_llk[j] <- exp(xb + tail)
#   }
#   res$site_llk <- site_llk
#   ### Notice here is not log likelihood since it's product 
#   res$read_llk <- prod(site_llk)
#   return(res)
# }

#' @description EM compute the mixture proportions and posterior probablities. Then update Hk.
#' @param n_observation number of observations
#' @param read_length a vector to store all of the observation length
#' @param par a list to store parameters
#' @param dat_short not expanded data
#' @param dat expanded data
#' @param n_class number of class (A1 A2 D1 D2 here).
#' @param num_cat number of category (for multinomial logistic regression ATCG here)
#' @return updated haplotypes, updated expanded dataframe, updated parameters(except for betas)

# Estep <- function(n_observation, read_length, par, dat_short, dat, n_class = 4, num_cat = 4) {
#   res <- list()
#   max_len <- max(read_length)
#   wic <- matrix(0, nrow = n_class, ncol = n_observation)
#   ## all_llk is the likelihood of each position of every read at each type of nucleotide (num_cat) and under each class
#   all_llk <- array(0, dim = c(n_class,  n_observation, num_cat, max_len))
#   ## read_class_llk is llk at each class for each read
#   read_class_llk <- matrix(0, ncol = n_class, nrow = n_observation)
#   
#   for (i in 1:n_observation) {
#     
#     subset <- dat_short %>% filter(idx == i)
#     ## eta * llk at each class for each read
#     weight_llk <- rep(0, n_class)
#     ## class_llk is likelihood of each position of each category under each class
#     class_llk <- array(0, dim = c(n_class, num_cat, max_len))
#     
#     for (k in 1:n_class) {
#       ### Extract haplotype
#       haplotype <- subset$hap_nuc[k + (0:(read_length[i] - 1)) * n_class]
#       llk_r <- llk(i = i, read_length = read_length, haplotype = haplotype, 
#                          dat = dat_short, par = par)
#       ###### compute the llk under ATCG (ONLY works for ATCG,  need to change to avoid hard coding)
#       llk_a <- llk(i = i, read_length = read_length, haplotype = 'A', 
#                    dat = dat_short, par = par, find_max = TRUE)
#       llk_t <- llk(i = i, read_length = read_length, haplotype = 'T', 
#                    dat = dat_short, par = par, find_max = TRUE)
#       llk_c <- llk(i = i, read_length = read_length, haplotype = 'C', 
#                    dat = dat_short, par = par, find_max = TRUE)
#       llk_g <- llk(i = i, read_length = read_length, haplotype = 'G', 
#                    dat = dat_short, par = par, find_max = TRUE)
#       
#       ## temp_llk store likelihood of each position and each nucleotide for each k and i 
#       temp_llk <- matrix(0, nrow = num_cat, ncol = max_len)
#       if(length(llk_r$site_llk) != max_len) {
#         temp_llk[1, ] <- c(llk_a$site_llk, rep(0, max_len - length(llk_a$site_llk)))
#         temp_llk[2, ] <- c(llk_c$site_llk, rep(0, max_len - length(llk_c$site_llk)))
#         temp_llk[3, ] <- c(llk_g$site_llk, rep(0, max_len - length(llk_g$site_llk)))
#         temp_llk[4, ] <- c(llk_t$site_llk, rep(0, max_len - length(llk_t$site_llk)))
#       } else {
#         temp_llk[1, ] <- llk_a$site_llk
#         temp_llk[2, ] <- llk_c$site_llk
#         temp_llk[3, ] <- llk_g$site_llk
#         temp_llk[4, ] <- llk_t$site_llk
#       }
#       class_llk[k, ,] <- temp_llk
#       read_class_llk[i, k] <- log(llk_r$read_llk) 
#       weight_llk[k] <- par$eta[k] * llk_r$read_llk # last_eta
#     }
#     all_llk[, i, ,] <- class_llk
#     wic[, i] <- (weight_llk)/sum(weight_llk)
#   }
#   ## Record the full likelihood for terminating EM
#   res$full_llk <- sum(wic * (t(read_class_llk) + log(par$eta)))
#   
#   par$eta <- rowMeans(wic)
#   par$wic <- wic
#   res$par <- par
#   
#   #### Find the most likely N for each j k (Should I chage to log likelihood??)
#   idx <- matrix(0, nrow = max_len, ncol = n_class)
#   for (j in 1:max_len) {
#     for (k in 1:n_class) {
#       lik <- rep(0, num_cat)
#       for (l in 1:num_cat) {
#         lik[l] <- sum(wic[k, ] %*% t(all_llk[k, , l, j]))
#       }
#       idx[j, k] <- which.max(lik) ## encode ACGT as 1234
#     }
#   }
#   ### Then update Hk, then length of each haplotype should be the same as the longest read? 
#   #### And when it matched back to dataset, should match according to the ref_pos?
#   hap <- matrix(0, ncol = n_class, nrow = max_len)
#   for (j in 1:max_len) {
#     for (k in 1:n_class) {
#       if(idx[j, k] == 1)
#         hap[j, k] <- "A"
#       else if(idx[j, k] == 2)
#         hap[j, k] <- "C"
#       else if(idx[j, k] == 3)
#         hap[j, k] <- "G"
#       else if(idx[j, k] == 4)
#         hap[j, k] <- "T"
#     }
#   }
#   
#   ##### Match new estimated haplotypes back to dataset
#   res$hap <- data.frame(hap) %>% `colnames<-`(paste0("hap", 1:n_class)) 
#   res$hap <- res$hap %>% 
#     add_column(ref_pos = 0:(nrow(res$hap) - 1)) %>% 
#     pivot_longer(cols = paste0("hap", 1:n_class), names_to = "hap_name", values_to = "hap_nuc")
#   res$dat <- dat %>% 
#     select(-c("hap_nuc")) %>% 
#     left_join(., res$hap, by = c("ref_pos", "hap_name"))
#   
#   return(res)
# }

################ TRY OTHER METHODS
## nnet got the same result  ## too slow
# fit_nnet <- nnet::multinom(formula(read_nuc ~ read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc),
#                            data = dat_long)
# fit_vglm <- vglm(read_nuc ~ read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc,
#                  multinomial, data = dat_long)
# TM <- mlogit.data(data, choice = "mode", shape = "long",
#                   alt.levels = c("A", "C", "G", "T"))
# fit_gmnl <- gmnl(formula(mode~1|read_pos + ref_pos + qua + hap_nuc + qua:hap_nuc|1),
#                  data = TM, model = "mnl")

  

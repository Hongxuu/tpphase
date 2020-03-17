fnlist <- function(x, fil){ 
  nams <- names(x) 
  for (i in seq_along(x)) {
    if(nams[i] == "haplotypes") {
      cat(nams[i], ":", "\n", file = fil, append = TRUE)
      write.table(x[[i]], file = fil, append = T, quote = FALSE, row.names=FALSE, col.names=FALSE)
    } else if(nams[i] == "snps") {
      cat(nams[i], ":", "\n", file = fil, append = TRUE)
      write.table(x[[i]], file = fil, append = T, quote = FALSE, row.names=FALSE, col.names=TRUE)
    } else {
      cat(nams[i], ":", "\n", x[[i]], "\n", file = fil, append = TRUE) 
    }
  }
}

dereplicate_res <- function(resu, haps, n_class) {
  final_res <- list()
  haplotypes <- matrix(to_char_r(haps$hap), nrow = n_class)
  idx <- duplicated(haplotypes)
  derepliacte_h <- haplotypes[!idx, ]
  #distinct(haplotypes)
  if(nrow(derepliacte_h) != nrow(haplotypes) & (nrow(derepliacte_h) != 1)) {
    flag <- which(idx == TRUE)
    ##TODO: hash talbe, this is only for the situation with four hap
    if (length(flag) == 1) {
      replicated <- c(flag - 1, flag)
      final_res$mixture_prop <- c(sum(resu$param$mixture_prop[replicated]), resu$param$mixture_prop[-replicated])
      weights <- cbind(rowSums(resu$param$w_ic[, replicated]), resu$param$w_ic[, -replicated])
    } else if ((length(flag) == 2)) {
      if((flag[2] - flag[1]) == 1) {
        replicated <- c(flag[1] - 1, flag)
        final_res$mixture_prop <- c(sum(resu$param$mixture_prop[replicated]), resu$param$mixture_prop[-replicated])
        weights <- cbind(rowSums(resu$param$w_ic[, replicated]), resu$param$w_ic[, -replicated])
      } else {
        final_res$mixture_prop <- c(sum(resu$param$mixture_prop[c(flag[1] - 1, flag[1])]), 
                                      sum(resu$param$mixture_prop[c(flag[2] - 1, flag[2])]))
        weights <- cbind(rowSums(resu$param$w_ic[, c(flag[1] - 1, flag[1])]),
                           rowSums(resu$param$w_ic[, c(flag[2] - 1, flag[2])]))
      }
    }
    final_res$assignments <- apply(weights, MARGIN = 1, FUN = which.max)
  } else {
    final_res$mixture_prop <- resu$param$mixture_prop
    final_res$assignments <- apply(resu$param$w_ic, MARGIN = 1, FUN = which.max)
  }
  cat(nrow(derepliacte_h), " haplotype(s) inferred\n")
  final_res$haplotypes <- derepliacte_h
  if (nrow(derepliacte_h) != 1) {
    snps_loci <- find_snp(hap = final_res$haplotypes) + 1
    final_res$snps <- final_res$haplotypes[, snps_loci] %>%  
      `colnames<-`(snps_loci)
    }
  final_res$logistic_coeff <- resu$param$beta
  final_res$full_llk <- resu$full_llk
  
  return(final_res)
}

to_char_r <- function(x) {
  as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G", "4" = "N")[as.character(x)])
}

to_xy_r <- function(x) {
  as.numeric(c("A" = "0", "T" = "2", "C" = "1", "G" = "3")[t(x)])
}

compare_par <- function(new, old, name) {
  mapply(FUN = function(A, B) {
    abs(A - B) < tol
  }, A = new[[name]], B = old[[name]]) %>% flatten_lgl()
} 

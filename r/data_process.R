fnlist <- function(x, fil){ 
  nams <- names(x) 
  for (i in seq_along(x)) {
    if(nams[i] == "haplotypes") {
      cat(nams[i], ":", "\n", file = fil, append = TRUE)
      write.table(x[[i]], file = fil, append = T, quote = FALSE, row.names=FALSE, col.names=FALSE)
    } else {
      cat(nams[i], ":", "\n", x[[i]], "\n", file = fil, append = TRUE) 
    }
  }
}

dereplicate_res <- function(resu, haps) {
  final_res <- list()
  haplotypes <- matrix(to_char_r(haps$hap), nrow = n_class)
  h <- data.table::data.table(haplotypes, key = paste0("V", 1:ncol(haplotypes)))
  derepliacte_h <- unique(h[duplicated(h)])
  
  if(nrow(derepliacte_h) != nrow(h)) {
    replicated <- lapply(1:nrow(derepliacte_h),function(i) h[derepliacte_h[i, ], which = T])
    final_res$mixture_prop <- sapply(1:nrow(derepliacte_h), function(i) sum(resu$param$mixture_prop[replicated[[i]]]))
    index <- sapply(replicated, "[[", 1)
    if (nrow(derepliacte_h) != 1)
      final_res$assignments <- apply(resu$param$w_ic[, index], MARGIN = 1, FUN = which.max)
  } else {
    final_res$mixture_prop <- resu$param$mixture_prop
    final_res$assignments <- apply(resu$param$w_ic, MARGIN = 1, FUN = which.max)
  }
  
  final_res$haplotypes <- derepliacte_h %>% as.matrix
  if (nrow(derepliacte_h) != 1)
    final_res$snps <- final_res$haplotypes[, find_snp(hap = final_res$haplotypes)]
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
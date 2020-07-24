library(prodlim)

read_res <- function(resultfileA = NULL, resultfileB = NULL, n_class = 4, select) {
  extract_hap <- function(resultfile, n_class, select = FALSE) {
    data <- read.delim(resultfile, header = FALSE, 
                       stringsAsFactors = FALSE) %>% `colnames<-`("name")
    mixture_prop <- data[which(data$name %in% "mixture_prop : ") + 1, ] %>% strsplit(" ") %>% unlist() %>% as.numeric() %>% na.omit()
    haps <- do.call(rbind, data[which(data$name %in% "haplotypes : ") + (1:n_class), ] %>% strsplit(" "))
    if(select == TRUE)
      haps[which((mixture_prop >= 0.1) == TRUE),]
    else
      haps
  }
  # if(!is.null(resultfileA)) {
  #   A <- extract_hap(resultfileA, n_class)
  #   return(A)
  # }
  # 
  # if(!is.null(resultfileB)) {
  #   B <- extract_hap(resultfileB, n_class)
  #   return(apply(A, 1, function(x) row.match(x, B)))
  # }
}

path = "../../data/tpphase_res_consensus/308TAN/308TAN_p46.txt"
hap46 <- extract_hap(resultfile = path, select = T, n_class = 4)

hap46[, find_snp(hap = hap38) + 1]
"../../data/tpphase/WGS/simu/L_SNP/sim0.fsa"
true_hap <- function(true_path, start_id, end_id) {
  t <- read_fasta(true_path)
  t_hap <- t$reads
  t_hap[, c(start_id:end_id)]
}
t_hap <- true_hap(true_path = "../../data/tpphase/WGS/simu/L_SNP/sim0.fsa", 
          start_id = dat_info$ref_start + 1, end_id = dat_info$ref_length_max + 1)


grasp <- function(id, start=0) {
  a <- c()
  for(i in 1:HMM$num_states[id]) {
    b <- hap_full_info$full_hap[[id]][[i]][, overlap_info$location[[id]] - HMM$time_pos[[id]] + 1+ start] 
    b <- as.vector(t(b))
    a <- rbind(a, b)
  }
  return(a)
}

grasp2 <- function(id, start=0) {
  for(i in 1:HMM$num_states[id]) {
    b <- hap_full_info$full_hap[[id]][[i]][, overlap_info$location[[id]] - HMM$time_pos[[id]] + 1+ start]
    print( b)
  }
}
id <- 18
grasp2(682, dat_info$ref_start)

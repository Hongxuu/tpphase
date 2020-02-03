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




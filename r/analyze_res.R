library(prodlim)

read_res <- function(resultfileA, resultfileB, n_class, select) {
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
  A <- extract_hap(resultfileA, 4)
  B <- extract_hap(resultfileB, n_class)
  apply(A, 1, function(x) row.match(x, B))
}

resultfileA = "../../../data/tpphase_res/308TAN_A_P40.txt"
resultfileB = "../../../data/tpphase_res/308TAN_B_P42.txt"
read_res(resultfileA, resultfileB, n_class = 4)





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
  if(!is.null(resultfileA)) {
    A <- extract_hap(resultfileA, n_class)
    return(A)
  }
  
  if(!is.null(resultfileB)) {
    B <- extract_hap(resultfileB, n_class)
    return(apply(A, 1, function(x) row.match(x, B)))
  }
}

path = "../../../data/tpphase_res_consensus"
data <- list()
dirs <- dir(path, full.names = TRUE)
filenames <- list.files(dirs, pattern = "*p*.txt", full.names=TRUE)[1:17]
haplotype <- list()
for(i in 1:length(filenames))
  haplotype[[i]] <- read_res(filenames[i], resultfileB = NULL, select = T)

haplotype[[7]][, 162]
haplotype[[7]][, 225]
haplotype[[7]][, 297]
haplotype[[7]][, 309]
haplotype[[7]][, 403]
haplotype[[7]] -> p3

haplotype[[8]][, 35]

for (i in 1:ncol(haplotype[[7]])) {
  if(length(unique(haplotype[[7]][, i])) != 1)
    cat(i, haplotype[[7]][, i], "\n")
}

54 T T G T 
87 C C A C 
162 C C A C 
225 A A T T 
237 G G C G 
297 G C G C 



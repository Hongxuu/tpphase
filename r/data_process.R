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

to_char_r <- function(x) {
  as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G", "4" = "N")[as.character(x)])
}

to_xy_r <- function(x) {
  as.numeric(c("A" = "0", "T" = "2", "C" = "1", "G" = "3")[t(x)])
}
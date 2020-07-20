######## Insertion: RefPos is -1; Deletion: Qua Nuc are -1 -

rea_dat <- function(filepath, filter = TRUE) {
  dat <- read.csv(filepath, sep = " ", header = F)
  names(dat) <- c("idx", "read_pos", "ref_pos", "qua", "read_nuc")
  #### Remove deletion and insertion:
  if(filter) {
    #dat <- dat %>% filter(!str_detect(qua, "-1"))
    dat <- dat %>% filter(!str_detect(ref_pos, "-1"))
  }
  dat %>% mutate_if(is.factor, as.character)
}

#### Initialize haplotype (ampliCI). Then expand the dataset.
# Mode <- function(x) {
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }
# dat <- dat %>% 
#   group_by(RefPos, refname) %>%
#   mutate(RefNuc = Mode(Nuc)) %>% 
#   ungroup

##### Find the number of reads and length in each read.

# n_observation <- function(dat) {
#   dat$idx <- 0
#   n <- 1
#   for (i in 1:500) {
#     if(dat$read_pos[i] > dat$read_pos[i+1]) {
#       n <- n + 1
#     }
#     dat$idx[i] <- n
#   }
# }

#n_observation(dat = dat)

###### Read fasta file: initial haplotype
##### All the haplotypes has the same length as the max length read:
###### M00259:19:000000000-AUAFV:1:2111:10994:25379	0	Adur105:427589_P3	1	60	251M250S	*	0	0

#reads <- readFastq("../data/test.fastq")
#init_hap@sread %>% width()
#reads@sread %>% width() %>% which.max()

clean_hap <- function(path, n_class = 4) {
  init_hap <- readFasta(path)
  haps <- init_hap@sread[1:n_class] %>% 
    str_split("") %>% 
    do.call(cbind, .) %>% 
    as.data.frame() %>% 
    `colnames<-`(paste0("hap", 1:n_class)) 
  haps <- haps %>% 
    add_column(ref_pos = 0:(nrow(haps)-1))
  haps
}

######### Merge two datasets and change the format:

format_data <- function(dat, haps, n_class = 4)  {
  dat1 <- left_join(dat, haps, by = "ref_pos")
  dat_long <- dat1 %>% 
    pivot_longer(cols = paste0("hap", 1:n_class), names_to = "hap_name", values_to = "hap_nuc") 
  dat_long %>% mutate_if(is.factor, as.character)
}


rep.data.frame <- function(x, times) {
  rnames <- attr(x, "row.names")
  x <- lapply(x, rep, each = times)
  class(x) <- "data.frame"
  if (!is.numeric(rnames))
    attr(x, "row.names") <- make.unique(rep.int(rnames, times))
  else
    attr(x, "row.names") <- .set_row_names(length(rnames) * times)
  x
}

lldat <- function(dat, dat_long, n_class = 4, num_cat = 4) {
  dat_rm <- dat_long %>% 
    select(-c("read_nuc")) %>% 
    rep.data.frame(num_cat)
  cat <- unique(dat_long$read_nuc)
  remain <- list()
  remain <- foreach (i = 1:length(dat$read_nuc)) %dopar% {
    left <- cat[-match(dat$read_nuc[i], cat)]
    order_nuc <- rep(c(dat$read_nuc[i], left), n_class)
    mode <- rep(rep(c(1, 0), c(1, (n_class - 1))), n_class)
    return(data.frame("read_nuc" = order_nuc, "mode" = mode))
  }
  remain_df <- as.data.frame(data.table::rbindlist(remain))
  bind_cols(dat_rm, remain_df)
}


### Main function

data_clean <- function(dat_path, hap_path, n_class = 4, num_cat = 4, filter = TRUE) {
  res <- list()
  data_r <- rea_dat(filepath = dat_path)
  haps <- clean_hap(path = hap_path)
  res$dat_long <- format_data(data_r, haps)
  res$lldat <- lldat(dat = data_r, dat_long = res$dat_long)
  return(res)
}

#lldat <- data_clean(dat_path = "test.txt", hap_path = "test_hap.fa")

write.table(dat_long, file = "long.txt", quote = FALSE, row.names=FALSE, col.names=FALSE)
data_r <- rea_dat(filepath = datafile)
datafile2 = "../../../../peanut_simu/homr0.005/cov3/out0a.txt"
data_r2 <- rea_dat(filepath = datafile2)
data_r %>% filter(!(idx %in% c(1:40) & ref_pos %in% c(150:200))) -> d
d <- d %>% group_by(idx) %>% mutate_at("read_pos", function(x) x - x[1])
data_r <- data_r %>% filter(idx <= 125)
write.table(data_r %>% filter(idx %in% c(1:10)), file = "test_30.txt", quote = FALSE, row.names=FALSE, col.names=FALSE)
##### add fake information for mnlogit format
data_rR <- data_r %>% filter(!(idx %in% c(1298, 1647, 1170, 1617, 1261, 1371, 220)))
write.table(d, file = "test_28.txt", quote = FALSE, row.names=FALSE, col.names=FALSE)
######### make fastq:
a <- vector()
for(i in 1:length(data_r$qua)) {
  a[i] <- intToUtf8(data_r$qua[i] + 33)
}
data_r$quality <- a
fil = "./test.fastq"
for(i in 1:length(unique(data_r$idx))) {
  subset <- data_r %>% filter(idx == i)
  cat(paste0("@", i), "\n", file = fil, append = TRUE)
  cat(subset$read_nuc, "\n", sep = "", file = fil, append = TRUE)
  cat("+", "\n", file = fil, append = TRUE)
  cat(subset$quality, "\n", sep = "", file = fil, append = TRUE)
}
######   idx   nuc   hap_nuc hap_name read_pos ref_pos qua mode
######    1     C       C       1        0       0     33   1
#####     1     C       C       2        0       0     33   1
######    1     C       A       3        0       0     33   1
#####     1     C       T       4        0       0     33   1
#####     1     A       T       4        0       0     33   0
####      1     A       C       2        0       0     33   0
######    1     A
#######   1     A
#######   1     G
#######   1     G
#######   1     G
#######   1     G
#######   1     T
#######   1     T
#######   1     T
#######   1     T

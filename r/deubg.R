# New one

beta1 <- resu[[7]]$param$beta
haplo1 <- haps[[5]] ##haps[[6]] == haps[[4]]
wic1 <- resu[[6]]$param$w_ic
data1 <- fromat_data(dat_info = d, haplotype = haplo1)
data1$nuc <- as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G")[as.character(data1$nuc)])
data1$hap_nuc <- as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G")[as.character(data1$hap_nuc)])
data1 <- data1[, !names(data1) %in% c("id")] # (modified_moligit exclude first column: id)
par1 <- list()
par1$wic <- t(wic1) #t(res$param$w_ic)
par1$beta <- resu[[6]]$param$beta #beta from last step as starting value

Mpar1 <- Mstep(dat = data1, read_length = read_length, par = par1, ncores = ncores)
###Old one

beta2 <- resu[[6]]$param$beta
haplo2 <- haps[[4]]
wic2 <- resu[[5]]$param$w_ic
data2 <- fromat_data(dat_info = d, haplotype = haplo2)
data2$nuc <- as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G")[as.character(data2$nuc)])
data2$hap_nuc <- as.character(c("0" = "A", "2" = "T", "1" = "C", "3" = "G")[as.character(data2$hap_nuc)])
data2 <- data2[, !names(data2) %in% c("id")] # (modified_moligit exclude first column: id)
par2 <- list()
par2$wic <- t(wic2)
par2$beta <- resu[[5]]$param$beta #beta from last step as starting value

Mpar2 <- Mstep(dat = data2, read_length = read_length, par = par2, ncores = ncores)
Mpar2$logLik


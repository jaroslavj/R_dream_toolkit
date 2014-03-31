setwd("/new_HD/jjelinek/dream112hc/mctables/sumtables")
getwd()
dir()
t1 <- read.table("t_dream11.1_hg19.txt", sep="\t", header=T)
t2 <- read.table("t_dream11.2_hg19.txt", sep="\t", header=T)
t12 <- t1[,4:10]+t2[,4:10]
head (t12)
mc91 <- read.table("mc9_dream11.1_hg19.txt", sep="\t", header=T)
mc92 <- read.table("mc9_dream11.2_hg19.txt", sep="\t", header=T)
mc912 <- {(t1[,4:10]*mc91[,4:10])+(t2[,4:10]*mc92[,4:10])}/t12
mc912 <- round(mc912, digits=2)
mc912f <-cbind(t1[,1:3],mc912)
write.table (mc912f, "mc9_HC_dream11_comb_NaN.txt", sep="\t",  row.names=F, quote=F)

mp1 <- read.table("mp_dream11.1_hg19.txt", sep="\t", header=T)
mp2 <- read.table("mp_dream11.2_hg19.txt", sep="\t", header=T)
mp12 <- {(t1[,4:10]*mp1[,4:10])+(t2[,4:10]*mp2[,4:10])}/t12

mp12 <- round(mp12, digits=2)
mp12f <-cbind(t1[,1:3],mp12)
write.table (mp12f, "mp_HC_dream11_comb_NaN.txt", sep="\t",  row.names=F, quote=F)

t12f <- cbind(t1[,1:3],t12)
write.table (t12f, "t_HC_dream11_comb_NaN.txt", sep="\t",  row.names=F, quote=F)


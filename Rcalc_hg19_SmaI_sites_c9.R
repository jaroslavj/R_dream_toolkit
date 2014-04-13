### Rcalc_hg19_SmaI_sites_c9.R ##
# ver 1.1
# JJ 04/12/2014
# fixed the wrong expected ln(m/u) values at 
# 	t268 to 4.7128
# 	SLa341 to -0.31955
# 	SLa374 to -0.0116
#
#
# Calculates 
#	
# 	uncorrected methylation 
#	c9 corrected methylation values
# Input file 1: _SmaI_sites.txt
# 	created as an output of "count_SmaI_CH3_c9.py" script
# removed # Input file 2: libraryName_spikes.txt
# removed #	created by an old script Rcalc_hum.R from smaple_meth.txt
# removed #	mapped by Shoudan's scripts
# removed #	takes c9 correction from the sample_spikes.txt file
# Output file: table with 
#	hg19_SmaI_ID, 
#	chromosome,
#	position, 
#	methylated reads, m
#	unmethylated reads, u
#	total reads, t
#	uncorrected methylation values, mp=100*(m+0.5)/(t+1)
#	c9-corrected meth values, mc9=100*c9*(m+0.5)/(c9*(m+0.5)+u+0.5)
#
# Methylation of spikes based on geometric means of mp values
# from dream8 to dream 13; n = 236 for c7 spikes, n=117 for c9 spikes
# La168=0.25; Eluc01=24.48 ; EGFP01=44.14; t245=t353=47.49; t268=t324=99.11; 
# not used t498=78.69; not used t883=44.63; not used t3307=53.08;
# SLa341=45.22, SLa341_137=39.00, SLa374_250=50.72, SLa374_130=48.70, Kan166=35.05


args  <- commandArgs(TRUE)

# read table __SmaI_sites.txt

met   <- read.table(args[1], sep="\t", header=T)
# create objects with methylated and unmethylated reads
m <- as.matrix(met[,4])
u <- as.matrix(met[,5])

# calculate total number of reads as t
t <- m+u

# calculate uncorrected methylation value
# avoid division by zero; if m=0 and u=0, mp=50
mp <- round(100*(m+0.5)/(t+1), digits=2)


# spike corrections calcs based on geometric means from >100 experiments
# see the above values
# calculate delta(expected_ln(m/u) - observed_ln((m+0.5)/(u+0.5))) for individual spikes
# calc average values from values at individual SmaI sites at spikes

# lambda corresponds to spike La168
lambda1 <- (-5.9978 -log((met[382290,4]+0.5)/(met[382290,5]+0.5)))
lambda2 <- (-5.9978 -log((met[382291,4]+0.5)/(met[382291,5]+0.5)))
lambda <- c(lambda1,lambda2)
avglambda <- mean(lambda)

# luc corresponds to spike 
luc1 <- (-1.1265 -log((met[382292,4]+0.5)/(met[382292,5]+0.5)))
luc2 <- (-1.1265 -log((met[382293,4]+0.5)/(met[382293,5]+0.5)))
luc <- c(luc1,luc2)
avgluc <- mean(luc)

# gfp corresponds to spike 
gfp1 <- (-0.2355 -log((met[382294,4]+0.5)/(met[382294,5]+0.5)))
gfp2 <- (-0.2355 -log((met[382295,4]+0.5)/(met[382295,5]+0.5)))
gfp <- c(gfp1, gfp2)
avggfp <- mean(gfp)

# t245 corresponds to spike T353
t245_1 <- (-0.1005 -log((met[382296,4]+0.5)/(met[382296,5]+0.5)))
t245_2 <- (-0.1005 -log((met[382297,4]+0.5)/(met[382297,5]+0.5)))
t245   <- c(t245_1, t245_2)
avgt245 <- mean(t245)

# t268 corresponds to spike T324
# 4.7128 is based on reference methylation 99.11%; ln(99.11/0.89)=4.712764
t268_1 <- (4.7128 -log((met[382298,4]+0.5)/(met[382298,5]+0.5)))
t268_2 <- (4.7128 -log((met[382299,4]+0.5)/(met[382299,5]+0.5)))
t268   <- c(t268_1, t268_2)
avgt268 <- mean(t268)

# Kan spike
kan1 <- (-0.6168 -log((met[382300,4]+0.5)/(met[382300,5]+0.5)))
kan2 <- (-0.6168 -log((met[382301,4]+0.5)/(met[382301,5]+0.5)))
kan   <- c(kan1,kan2)
avgkan <- mean(kan)

# SLa341 spike
# -0.31955 is the average of 45.22% (ln m/u= -0.19179) and 39% (ln m/u= -0.44731) 
SLa341_1 	<- (-0.31955 -log((met[382302,4]+0.5)/(met[382302,5]+0.5)))
SLa341_2 	<- (-0.31955 -log((met[382303,4]+0.5)/(met[382303,5]+0.5)))
SLa341_3 	<- (-0.31955 -log((met[382304,4]+0.5)/(met[382304,5]+0.5)))
sla341   	<- c(SLa341_1, SLa341_2, SLa341_3)
avgsla341	<- mean(sla341)

# SLa374 spike
# -0.0116 is the average of 50.72% (ln m/u= 0.028802) and 48.7% (ln m/u= -0.05201) 
SLa374_1 	<- (-0.0116 -log((met[382305,4]+0.5)/(met[382305,5]+0.5)))
SLa374_2 	<- (-0.0116 -log((met[382306,4]+0.5)/(met[382306,5]+0.5)))
SLa374_3 	<- (-0.0116 -log((met[382307,4]+0.5)/(met[382307,5]+0.5)))
sla374   	<- c(SLa374_1, SLa374_2, SLa374_3)
avgsla374	<- mean(sla374)

# average for differences from median ln(m/u) values for c9 spikes 
delta.c9spikes <- mean(c(avgluc, avggfp, avgt245, avgkan, avgsla341, avgsla374))
c9 <- exp(delta.c9spikes)

# average for differences from median ln(m/u) values for c10 spikes 
delta.c10spikes <- mean(c(avglambda, avgluc, avggfp, avgt245, avgt268, avgkan, avgsla341, avgsla374))
c10 <- exp(delta.c10spikes)

# average for differences from median ln(m/u) values for c3 spikes 
delta.c3spikes <- mean(c(avgluc, avggfp, avgt245))
c3 <- exp(delta.c3spikes)

# calculation of methylation values
mc3 <- round((100*c3*(m+0.5))/(c3*(m+0.5)+u+0.5), digits=2)
mc9 <- round((100*c9*(m+0.5))/(c9*(m+0.5)+u+0.5), digits=2)
mc10 <- round((100*c10*(m+0.5))/(c10*(m+0.5)+u+0.5), digits=2)

# cbind columns m, u, t, mp, mc9 and create column names
metc <- cbind(m, u, t, mp, mc9)
colnames (metc) <- c("m", "u", "t", "mp", "mc9")

# read table with SmaI IDs, chromosome, position 
#	and c9 correction in the last line
hg19_ID <- read.table("ID_hg19_SmaI_sites_c9.txt",sep="\t", header=T)

# 	make a matrix with c3,c9,c10 correction values
corrections <- round(c(c3,c9,c10), digits=3)
correctmatrix   <- matrix(corrections, nrow=3, ncol=5)

# bind all tables together and write an output table
metcorr <- rbind (metc,correctmatrix)
metcalc <- cbind(hg19_ID, metcorr)
outputfile <- paste(args[1], "_mc.txt", sep="")
write.table(metcalc, outputfile, sep="\t", row.names=F, quote=F)

# create a table of spikes and corrections
spikes <- (metcalc [382290:382307,])
spikesout <- paste(args[1], "_spikes.txt", sep="")
write.table (spikes, spikesout, sep="\t", row.names=F, quote=F)



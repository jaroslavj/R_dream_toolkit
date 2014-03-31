### mint10_ttable_mc9.R ###
# 03-31-2014 JJ #


# The script requires two tables with identical structure
# [1] table with the number of reads (t)
# [2] table with methylation values
# It uses the ifelse function to create a temporary table 'mint' with 0 for t<10 and 1 for t>=10.
# Next it divides the methylation values in 'mc9' table by the values in the 'mint' table.
# Methylation values divided by zero are converted to "Inf", methylation values divided by 1 do not change.
# The last step creates a new table where methylation values based on <10 reads are replaced by "Inf" . 


# loop for setting up minimum ten (10) reads
# t <- ttable[,4:dim(ttable)[2]]
# ttable <- read.table ("table with tvalues " sep="\t",header=T)
# mc9 <- read.table ("table with tvalues " sep="\t",header=T)

# begin script

args  <- commandArgs(TRUE)
ttable <- read.table (args[1], sep="\t",header=T)
mc9 <- read.table (args[2], sep="\t",header=T)


# starting the loop
for (i in 4:dim(ttable)[2]) 
{
if (i==4) {mint <-ifelse(ttable[,i]<10,0,ttable[,i])
table <-cbind(ttable[,1:3],mint)}

if (i>4) {mint <-ifelse(ttable[,i]<10,0,ttable[,i])
table <-cbind(table,mint)}
}
mint.final <- table[,4:dim(ttable)[2]]

mc9mint <- round(mc9[,4:dim(mc9)[2]]/mint.final, digits=2)
outputfile <- paste(args[2], "_mint10.txt", sep="")
# write.table(cbind(mc9[,1:3],mc9mint), "mc9mint10.txt",sep="\t",row.names=F)

write.table(cbind(mc9[,1:3],mc9mint), outputfile,sep="\t",row.names=F)

# end of script

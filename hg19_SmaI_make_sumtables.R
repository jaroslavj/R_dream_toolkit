### hg19_SmaI_make_sumtables.R ###

# modif JJ 03/23/2014

# The script creates merged t, mp, mc9 tables for a sampleSet from individual samples
# And deposits new tables into subdirectory tables/sumtables/
# Works for "_SmaI_sites.txt_mc9.txt" files created as an output of "Rcalc_hg19_SmaI_sites.R"


# Required input:  "_SmaI_sites.txt_mc9.txt" files in a subdirectory mctables/
# Makes output files:
#	sumtables/sampleSet_hg19_t.txt
#	sumtables/sampleSet_hg19_mp.txt
#	sumtables/sampleSet_hg19_mc9.txt

# Before running the script 
# Copy all ..._SmaI_sites.txt_mc9.txt files in a new directory mctables/

# Change the current directory to be above mctables/
# Open R and copy/paste the script below


setwd ("mctables")
# sets working directory to "mctables"

dir.create ("sumtables")
# creates a subdirectory "sumtables"

# this script makes a table with the total number of reads, t

file.name=dir()     #read all files name
length(file.name)   #check the dimension how many files you have
for (i in 1:length(file.name)){                                 #here start the loop
  table.file=read.table(file.name[i],sep="\t",header=T)                               #read the table number i
    if (i==1){                                                                        #first cycle
      final.table=as.data.frame(table.file[,c(1,2,3,6)])                                     #   create dataframe
      vec=as.character(c("hg19_ID", "chrom", "position", unlist(strsplit(file.name[i], "_SmaI"))[1]))         # and vector with column names
    }
  if (i>1){                                                                             #all the other cycles
    final.table=cbind(final.table,table.file[,6])                                         #bind columns together
    vec=c(vec,unlist(strsplit(file.name[i], "_SmaI"))[1])                                                 #bind colnames together
  }
}
colnames(final.table)=vec #assign colnames to the dataframe
write.table(final.table, "sumtables/sampleSet_hg19_t.txt",sep="\t",row.names=F)

### end of the script

# this script makes a table with uncorrected methylation values, mp

file.name=dir()     #read all files name
length(file.name)   #check the dimension how many files you have
for (i in 1:length(file.name)){                                 #here start the loop
  table.file=read.table(file.name[i],sep="\t",header=T)                               #read the table number i
    if (i==1){                                                                        #first cycle
      final.table=as.data.frame(table.file[,c(1,2,3,7)])                                     #   create dataframe
      vec=as.character(c("hg19_ID", "chrom", "position", unlist(strsplit(file.name[i], "_SmaI"))[1]))         # and vector with column names
    }
  if (i>1){                                                                             #all the other cycles
    final.table=cbind(final.table,table.file[,7])                                         #bind columns together
    vec=c(vec,unlist(strsplit(file.name[i], "_SmaI"))[1])                                                 #bind colnames together
  }
}
colnames(final.table)=vec #assign colnames to the dataframe
write.table(final.table, "sumtables/sampleSet_hg19_mp.txt",sep="\t",row.names=F)

### end of mp script


# this script makes a table with c9 corrected methylation values, mc9

file.name=dir()     #read all files name
length(file.name)   #check the dimension how many files you have
for (i in 1:length(file.name)){                                 #here start the loop
  table.file=read.table(file.name[i],sep="\t",header=T)                               #read the table number i
    if (i==1){                                                                        #first cycle
      final.table=as.data.frame(table.file[,c(1,2,3,8)])                                     #   create dataframe
      vec=as.character(c("hg19_ID", "chrom", "position", unlist(strsplit(file.name[i], "_SmaI"))[1]))         # and vector with column names
    }
  if (i>1){                                                                             #all the other cycles
    final.table=cbind(final.table,table.file[,8])                                         #bind columns together
    vec=c(vec,unlist(strsplit(file.name[i], "_SmaI"))[1])                                                 #bind colnames together
  }
}
colnames(final.table)=vec #assign colnames to the dataframe
write.table(final.table, "sumtables/sampleSet_hg19_mc9.txt",sep="\t",row.names=F)

### end of mc9 script
### summary tables are in subdirectory sumtables/



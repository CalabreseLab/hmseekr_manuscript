# organize multicov results for eCLIP data
setwd("XNM")

library(tidyverse)
library(ggplot2)

header<-read.table('multicov_files_list_XNM.txt',sep=' ')
# 140 files in total

header$V1<-sub("_.*", "", header$V1)


# normalize the reads with total counts in the sequencing fastq files
fastqcount<-read.table('readcounts_pre_filter_eCLIP_11_7_24.txt',sep=' ')
colnames(fastqcount)<-c('FileName','TotalReads')
fastqcount$FileName<-sub("_.*", "", fastqcount$FileName)


# organize fastqcount in the same order of cols in _mc
fastqorder<-header$V1
fastqcount.ordered<-vector(mode='numeric',length=length(fastqorder))

for (n in 1:length(fastqorder)) {
  fname<-fastqorder[n]
  fc<-fastqcount$TotalReads[which(fastqcount$FileName==fname)]
  fastqcount.ordered[n]<-fc
}


# remove extra info and create a header
headers<-header$V1


f.temp<-read.table('XNM_chunk_bedfile_multicov.out',sep='\t')

colnames(f.temp)<-c('chrom','chromstart','chromend','name','score','strand',headers)

# get rpm count
f_norm<- sweep(f.temp[,c(7:146)], MARGIN=2, STATS=fastqcount.ordered, FUN="/")
f_norm<- f_norm*1000000
f_norm<-cbind(f.temp[,c(1:6)],f_norm)


int9<-f_norm[which(f_norm$name=='XISTint9'),]
int9data<-int9[,c(7:145)]
int9colsum<-colSums(int9data)

f_norm$chromstart[14]<-min(int9$chromstart)
f_norm[14,c(7:145)]<-int9colsum

f_norm<-f_norm[-c(15:18),]

f_comb<-f_norm
# subtract ck from all samples
# if negative after subtraction, set to 0
f_comb[,7:146] <- apply(f_comb[,7:146], 2, function(x) x - f_comb$control)
f_comb[,7:146] <- as.data.frame(lapply(f_comb[,7:146], function(x) ifelse(x < 0, 0, x)))

f_comb$control<-NULL


########
# get the top 5 ranked RBP for each chunk

rbpdata<-f_comb

rbpdata$top5RBP<-''

rbpdata$top5rpm<-''

rbpdata$top5fc<-''

rbpdatafiltered<-rbpdata

for (n in 1:nrow(rbpdata)) {
  
  temp<-rbpdata[n,c(7:145)]
  
  temp.ori<-f_norm[n,c(7:146)]
  
  if (all(temp==0)) {
    
    rbpdata$top5RBP[n]<-''
    rbpdata$top5rpm[n]<-''
    rbpdata$top5fc[n]<-''
    
    rbpdatafiltered$top5RBP[n]<-''
    rbpdatafiltered$top5rpm[n]<-''
    rbpdatafiltered$top5fc[n]<-''
    
    next
  } 
  
  tempck<-temp.ori$control
  
  temp.ori$control<-NULL
  
  temp.ori<-temp.ori[order(as.numeric(temp),decreasing = T)]
  
  temp<-temp[order(as.numeric(temp),decreasing = T)]
  
  univalue <- unique(as.numeric(temp)) # Get unique values
  
  # Filter for positive values
  positive_values <- univalue[univalue > 0]
  
  # Determine the threshold
  if (length(positive_values) >= 5) {
    threshold <- positive_values[5] # Get the 5th ranked value
  } else {
    threshold <- min(positive_values) # Get the minimum positive value
  }
  

  
  temp.ori<-temp.ori[which(as.numeric(temp)>=threshold)]
  
  temp<-temp[which(as.numeric(temp)>=threshold)]
  
  t5rbpnames<-names(temp)
  
  rbpdata$top5RBP[n]<-paste(t5rbpnames,collapse=';')
  
  t5rpm<-round(as.numeric(temp),digits = 4)
  
  rbpdata$top5rpm[n]<-paste(t5rpm,collapse=';')
  
  t5fc<-round(as.numeric(temp.ori)/as.numeric(tempck),digits = 4)
  
  rbpdata$top5fc[n]<-paste(t5fc,collapse=';')
  
  # filter the recorded RBP based on rpm and fc
  
  flidx<-which(t5rpm>2 & t5fc>2)
  
  rbpdatafiltered$top5RBP[n]<-paste(t5rbpnames[flidx],collapse=';')
  rbpdatafiltered$top5rpm[n]<-paste(t5rpm[flidx],collapse=';')
  rbpdatafiltered$top5fc[n]<-paste(t5fc[flidx],collapse=';')
  
}

rbpdata<-rbpdata[,c(1:6,146:148,7:145)]

rbpdatafiltered<-rbpdatafiltered[,c(1:6,146:148,7:145)]

write.csv(rbpdatafiltered,'XNM_RBP_rpm_counts_CKnormed_chunks_filtered.csv',row.names = F)



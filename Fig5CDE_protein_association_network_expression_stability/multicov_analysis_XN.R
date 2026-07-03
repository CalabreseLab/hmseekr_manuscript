# organize multicov results and calculate pearson correlation between proteins
setwd("eCLIP_peakbams")
library(tidyverse)
library(ggplot2)

header<-read.table('multicov_files_list_XNpattern.txt',sep=' ')


header$V1<-gsub('_peaks.bam','',header$V1,fixed=T)

# normalize the reads with total counts in the sequencing fastq files
rawcount<-read.table('readcounts_pre_filter_eCLIP_11_7_24_all.txt',header=F,sep=':')
rawcount$V1<-gsub('_strand-specific.*','',rawcount$V1)
countdf<-read.csv('readcounts_eCLIP_110724.csv',header=T)

for (n in 1:nrow(countdf)){
  
  expname<-countdf$experiment[n]
  expname<-gsub('_strand-specific.*','',expname)
  
  ctrlname<-countdf$control[n]
  ctrlname<-gsub('_strand-specific.*','',ctrlname)
  
  expid<-which(rawcount$V1==expname)
  ctrlid<-which(rawcount$V1==ctrlname)
  
  if (length(expid)>1) {break}
  if (length(ctrlid)>1) {break}
  
  countdf$exp_count[n]<-as.numeric(rawcount$V2[expid])
  countdf$ctrl_count[n]<-as.numeric(rawcount$V2[ctrlid])
  
}

countdf$exp_name<-paste0(countdf$RBP,'_eCLIP')
countdf$ctrl_name<-paste0(countdf$RBP,'_control')

write.csv(countdf,'readcounts_eCLIP_110724.csv',row.names = F)

countdf1<-countdf[,c('exp_name','exp_count')]
countdf2<-countdf[,c('ctrl_name','ctrl_count')]
colnames(countdf1)<-c('name','count')
colnames(countdf2)<-c('name','count')

countdf<-rbind(countdf1,countdf2)
rm(countdf1)
rm(countdf2)

# organize fastqcount in the same order of cols in _mc
fastqorder<-header$V1
fastqcount.ordered<-vector(mode='numeric',length=length(fastqorder))


for (n in 1:length(fastqorder)) {
  fname<-fastqorder[n]
  fc<-countdf$count[which(countdf$name==fname)]
  fastqcount.ordered[n]<-fc
}


lowtriangle<-function(df) {
  
  df[upper.tri(df,diag=T)]<-NA
  
  df<-as.data.frame(df)
  
  # add rownames as one column
  df$RBP_1<-rownames(df)
  
  lt<-pivot_longer(df,cols=-RBP_1, names_to='RBP_2', values_to = 'cor_r', values_drop_na = T)
  
  lt<-as.data.frame(lt)
  
  return(lt)
}


mfiles<-dir('./multicov_out/',pattern='\\.out$')


# add counter to check if all multicov_out file has 70 columns (6 properties + 64 seq files)
mfile.err<-0

for (n in 1:length(mfiles)) {
  if (n %% 100==0){print(n)}
  
  fn<-mfiles[n]
  
  f.temp<-read.table(paste0('./multicov_out/',fn),sep='\t')
  

  # check property number
  if (ncol(f.temp)!=70) {mfile.err<-mfile.err+1}
  
  colnames(f.temp)<-c('chrom','chromstart','chromend','name','score','strand',header$V1)
  
  # get rpm count
  f_norm<- sweep(f.temp[,c(7:70)], MARGIN=2, STATS=fastqcount.ordered, FUN="/")
  f_norm<- f_norm*1000000
  f_norm<-cbind(f.temp[,c(1:6)],f_norm)
  
  allrbp <- header$V1[grepl('_eCLIP',header$V1)]
  allrbp<-gsub('_eCLIP','',allrbp)
  
  f_comb<-f_norm[,c(1:6)]
  
  for (rbp in allrbp){
    
    ename<-paste0(rbp,'_eCLIP')
    cname<-paste0(rbp,'_control')
    nname<-paste0(rbp,'_norm')
    
    f_comb[[nname]]<-f_norm[[ename]]-f_norm[[cname]]
    f_comb[[nname]][which(f_comb[[nname]]<0)]<-0
    
  }
  

  # select columns from 7 to 38
  f_data <- f_comb[,7:38]
  
  # compute the pairwise correlations using all complete pairs of observations
  # for each pair of columns, if there are NA entries the result will be NA
  # there is no NA in all datasets
  # NA is produced for RBP that has all 0 counts for a gene for all its cor with others
  f_cor <- cor(f_data, method = "pearson", use = "pairwise.complete.obs")

  
  ################# integrate all RBP profiles #########################
  # melt the correlation matrix to long format 
  # only keep the upper triangle and exclude the diagonal
  
  # Get indices of the upper triangle (excluding the diagonal)
  findices <- which(upper.tri(f_cor,diag = FALSE), arr.ind = TRUE)
  
  # Build data frame
  f_edges <- data.frame(
    RBP1 = rownames(f_cor)[findices[,1]],
    RBP2 = rownames(f_cor)[findices[,2]],
    Weight = f_cor[findices]
  )

  
  ltname<-gsub('_multicov.out','',fn)
  
  
  if (exists('comb_allgenes')) {
    if(all(comb_allgenes$RBP1==f_edges$RBP1) & all(comb_allgenes$RBP2==f_edges$RBP2)) {
      comb_allgenes[[ltname]]<-f_edges$Weight
    } else {
      print(ltname)
      print('RBP pairs do not match')
    }
    
  } else {
    colnames(f_edges)[3]<-ltname
    comb_allgenes<-f_edges
  }
  
}

write.csv(comb_allgenes,'./comb_allgenes_cor_XNpattern.csv',row.names = F)



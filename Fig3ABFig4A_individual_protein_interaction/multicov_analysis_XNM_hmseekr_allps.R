# process multicov results for XNM hmseekr hits vs control

setwd("XNM/seqstosummary")
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)


#########################################
# unspliced
##########################################

###################################
# normalize bams with control p 0.05
###################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

bamcounts<-read.table('readcounts_pre_filter_eCLIP_11_7_24.txt',sep = ' ',header=F)

bamcounts$V1<-sub("_.*", "", bamcounts$V1)

colnames(bamcounts)<-c('Filename','Total_Counts')


hitfiles<-list.files(path = './multicov/', pattern = "^tsc_condE_newpool_unspliced_.*_hits_multicov\\.out$", full.names = FALSE)

ckfiles<-list.files(path = './multicov/', pattern = "^tsc_condE_newpool_unspliced_.*_ck_multicov\\.out$", full.names = FALSE)



custom.filter<-function(x,y){
  temp<-as.data.frame(cbind(x,y))
  return(temp)
}


hmmstat<-as.data.frame(matrix(nrow=51,ncol=139))
xheaders <- c('XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
              'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
              'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14') 

nheaders<-paste0('NEAT1chk',c(1:22))
mheaders<-paste0('MALAT1chk',c(1:9))

xnm<-c(xheaders,nheaders,mheaders)

rownames(hmmstat)<-xnm


rbps<-read.table('multicov_files_list_XNM.txt',sep='\n')

rbps<-as.vector(rbps$V1)

rbps<-sub("_.*", "", rbps)

colnames(hmmstat)<-rbps[-12]

hmmcount<-hmmstat
hmmtcount<-hmmstat
hmmmedlen<-hmmstat

rip.vec<-rbps
ck.vec<-paste0(rbps,'_ck')

# remove control
rip.vec<-rip.vec[-12]
ck.vec<-ck.vec[-12]

# make sure the order in rbps are the same as in bamcounts
bamcounts <- bamcounts[match(rbps, bamcounts$Filename), ]

# counts file corresponds to multicov output cols
sum(rbps!=bamcounts$Filename) # 0

for (n in 1:length(hitfiles)) {
  
  print(n)
  
  hfile<-hitfiles[n]
  cfile<-ckfiles[n]
  
  # check if file is empty before processing
  # if empty store NA in stats and continue to the next n in loop
  filepath <- paste0('./multicov/', hfile)
  if (file.info(filepath)$size == 0) {
    print(paste0(hfile,' is empty'))
    
    chunk<-gsub('_4_viterbi_seekr_hits_multicov.out','',hfile)
    chunk<-gsub('tsc_condE_newpool_unspliced_','',chunk)
    # which row to store the stat
    rp<-which(xnm==chunk)
    
    hmmstat[rp,]<-NA
    hmmcount[rp,]<-0
    hmmtcount[rp,]<-0
    hmmmedlen[rp,]<-NA
    
    next
  }
  
  
  hitdata <- read.table(filepath, sep='\t')
  ckdata<-read.table(paste0('./multicov/',cfile),sep='\t')
  
  colnames(hitdata)<-c('chrom','chromStart','chromEnd','name','score','strand',rbps)
  colnames(ckdata)<-paste0(c('chrom','chromStart','chromEnd','name','score','strand',rbps),'_ck')
  
  hitdata$tlen<-hitdata$chromEnd-hitdata$chromStart
  ckdata$tlen_ck<-ckdata$chromEnd_ck-ckdata$chromStart_ck
  
  print(paste0('hitdata and ckdata has the same nrow_',nrow(hitdata)==nrow(ckdata)))
  
  print(paste0('hitdata and ckdata has the different length_',
               sum(hitdata$tlen!=ckdata$tlen_ck)))
  
  # get rpm count
  hitnorm<- sweep(hitdata[,c(7:146)], MARGIN=2, STATS=as.numeric(bamcounts$Total_Counts), FUN="/")
  hitnorm<- hitnorm*1000000
  hitnorm<-cbind(hitdata[,c(1:6)],hitnorm)
  hitnorm$tlen<-hitdata$tlen
  
  cknorm<- sweep(ckdata[,c(7:146)], MARGIN=2, STATS=as.numeric(bamcounts$Total_Counts), FUN="/")
  cknorm<- cknorm*1000000
  cknorm<-cbind(ckdata[,c(1:6)],cknorm)
  cknorm$tlen_ck<-ckdata$tlen_ck
  
  # subtract control from all samples
  # if negative after subtraction, set to 0
  hitnorm[,7:146] <- apply(hitnorm[,7:146], 2, function(x) x - hitnorm$control)
  hitnorm[,7:146] <- as.data.frame(lapply(hitnorm[,7:146], function(x) ifelse(x < 0, 0, x)))
  # drop the control col
  hitnorm$control<-NULL
  
  cknorm[,7:146] <- apply(cknorm[,7:146], 2, function(x) x - cknorm$control_ck)
  cknorm[,7:146] <- as.data.frame(lapply(cknorm[,7:146], function(x) ifelse(x < 0, 0, x)))
  
  # drop the control col
  cknorm$control_ck<-NULL
  
  
  comb<-cbind(hitnorm,cknorm)
  
  
  combname<-gsub('_4_viterbi_seekr_hits_multicov.out','',hfile)
  
  combname<-paste0('./combined/',combname,'_combined_CKnormed_p05.csv')
  
  write.csv(comb,combname,row.names = F)
  
  #########################################
  # filter the comb data
  # only filter on the score not score_ck
  # only filter on seekr p val score of hits not the control
  comb_fl<-comb[which(comb$score<0.05),]
  
  #######################################
  
  
  chunk<-gsub('_4_viterbi_seekr_hits_multicov.out','',hfile)
  chunk<-gsub('tsc_condE_newpool_unspliced_','',chunk)
  # which row to store the stat
  rp<-which(xnm==chunk)
  
  if (nrow(comb_fl)>0) {
    
    hmmmedlen[rp,]<-median(comb_fl$tlen,na.rm=T)
    hmmcount[rp,]<-nrow(comb_fl)
    
    thit<-sub("_[0-9]+$", "", comb_fl$name)
    hmmtcount[rp,]<-length(unique(thit))
    
    for (m in 1:length(rip.vec)) {
      
      rip<-rip.vec[m]
      ck<-ck.vec[m]
      
      temp<-custom.filter(comb_fl[rip],comb_fl[ck])
      
      hmmstat[rp,m]<-wilcox.test(temp[,1],temp[,2],alternative = 'greater',paired=T)$p.value
      
      
    }
  } else {
    hmmstat[rp,]<-NA
    hmmcount[rp,]<-0
    hmmtcount[rp,]<-0
    hmmmedlen[rp,]<-NA
  }
  
  
  
  
}

write.csv(hmmstat,'XNM_hmseekrhits_paired_unspliced_pval_CKnormed_p05.csv')
write.csv(hmmcount,'XNM_hmseekrhits_paired_unspliced_ncount_CKnormed_p05.csv')
write.csv(hmmtcount,'XNM_hmseekrhits_paired_unspliced_transcript_count_CKnormed_p05.csv')
write.csv(hmmmedlen,'XNM_hmseekrhits_paired_unspliced_medlen_CKnormed_p05.csv')



###################################
# normalize bams with control p 0.01
###################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

bamcounts<-read.table('readcounts_pre_filter_eCLIP_11_7_24.txt',sep = ' ',header=F)

bamcounts$V1<-sub("_.*", "", bamcounts$V1)

colnames(bamcounts)<-c('Filename','Total_Counts')


hitfiles<-list.files(path = './multicov/', pattern = "^tsc_condE_newpool_unspliced_.*_hits_multicov\\.out$", full.names = FALSE)

ckfiles<-list.files(path = './multicov/', pattern = "^tsc_condE_newpool_unspliced_.*_ck_multicov\\.out$", full.names = FALSE)



custom.filter<-function(x,y){
  
  temp<-as.data.frame(cbind(x,y))
  return(temp)
}


hmmstat<-as.data.frame(matrix(nrow=51,ncol=139))

xheaders <- c('XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
              'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
              'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14') 

nheaders<-paste0('NEAT1chk',c(1:22))
mheaders<-paste0('MALAT1chk',c(1:9))

xnm<-c(xheaders,nheaders,mheaders)

rownames(hmmstat)<-xnm


rbps<-read.table('multicov_files_list_XNM.txt',sep='\n')

rbps<-as.vector(rbps$V1)

rbps<-sub("_.*", "", rbps)

colnames(hmmstat)<-rbps[-12]

hmmcount<-hmmstat
hmmtcount<-hmmstat
hmmmedlen<-hmmstat

rip.vec<-rbps
ck.vec<-paste0(rbps,'_ck')

# remove control
rip.vec<-rip.vec[-12]
ck.vec<-ck.vec[-12]
  
# make sure the order in rbps are the same as in bamcounts
bamcounts <- bamcounts[match(rbps, bamcounts$Filename), ]

# counts file corresponds to multicov output cols
sum(rbps!=bamcounts$Filename) # 0


for (n in 1:length(hitfiles)) {
  
  print(n)
  
  hfile<-hitfiles[n]
  cfile<-ckfiles[n]
  
  # check if file is empty before processing
  # if empty store NA in stats and continue to the next n in loop
  filepath <- paste0('./multicov/', hfile)
  if (file.info(filepath)$size == 0) {
    print(paste0(hfile,' is empty'))
    
    chunk<-gsub('_4_viterbi_seekr_hits_multicov.out','',hfile)
    chunk<-gsub('tsc_condE_newpool_unspliced_','',chunk)
    # which row to store the stat
    rp<-which(xnm==chunk)
    
    hmmstat[rp,]<-NA
    hmmcount[rp,]<-0
    hmmtcount[rp,]<-0
    hmmmedlen[rp,]<-NA
    
    next
  }
  
  
  hitdata <- read.table(filepath, sep='\t')
  ckdata<-read.table(paste0('./multicov/',cfile),sep='\t')
  
  colnames(hitdata)<-c('chrom','chromStart','chromEnd','name','score','strand',rbps)
  colnames(ckdata)<-paste0(c('chrom','chromStart','chromEnd','name','score','strand',rbps),'_ck')
  
  hitdata$tlen<-hitdata$chromEnd-hitdata$chromStart
  ckdata$tlen_ck<-ckdata$chromEnd_ck-ckdata$chromStart_ck
  
  print(paste0('hitdata and ckdata has the same nrow_',nrow(hitdata)==nrow(ckdata)))
  
  print(paste0('hitdata and ckdata has the different length_',
               sum(hitdata$tlen!=ckdata$tlen_ck)))
  
  # get rpm count
  hitnorm<- sweep(hitdata[,c(7:146)], MARGIN=2, STATS=as.numeric(bamcounts$Total_Counts), FUN="/")
  hitnorm<- hitnorm*1000000
  hitnorm<-cbind(hitdata[,c(1:6)],hitnorm)
  hitnorm$tlen<-hitdata$tlen
  
  cknorm<- sweep(ckdata[,c(7:146)], MARGIN=2, STATS=as.numeric(bamcounts$Total_Counts), FUN="/")
  cknorm<- cknorm*1000000
  cknorm<-cbind(ckdata[,c(1:6)],cknorm)
  cknorm$tlen_ck<-ckdata$tlen_ck
  
  # subtract control from all samples
  # if negative after subtraction, set to 0
  hitnorm[,7:146] <- apply(hitnorm[,7:146], 2, function(x) x - hitnorm$control)
  hitnorm[,7:146] <- as.data.frame(lapply(hitnorm[,7:146], function(x) ifelse(x < 0, 0, x)))
  # drop the control col
  hitnorm$control<-NULL
  
  cknorm[,7:146] <- apply(cknorm[,7:146], 2, function(x) x - cknorm$control_ck)
  cknorm[,7:146] <- as.data.frame(lapply(cknorm[,7:146], function(x) ifelse(x < 0, 0, x)))
  
  # drop the control col
  cknorm$control_ck<-NULL
  
  
  comb<-cbind(hitnorm,cknorm)
  
  
  combname<-gsub('_4_viterbi_seekr_hits_multicov.out','',hfile)
  
  combname<-paste0('./combined/',combname,'_combined_CKnormed_p01.csv')
  
  write.csv(comb,combname,row.names = F)
  
  #########################################
  # filter the comb data
  # only filter on the score not score_ck
  # only filter on seekr p val score of hits not the control
  comb_fl<-comb[which(comb$score<0.01),]
  
  #######################################
  
  
  chunk<-gsub('_4_viterbi_seekr_hits_multicov.out','',hfile)
  chunk<-gsub('tsc_condE_newpool_unspliced_','',chunk)
  # which row to store the stat
  rp<-which(xnm==chunk)
  
  if (nrow(comb_fl)>0) {
    
    hmmmedlen[rp,]<-median(comb_fl$tlen,na.rm=T)
    hmmcount[rp,]<-nrow(comb_fl)
    
    thit<-sub("_[0-9]+$", "", comb_fl$name)
    hmmtcount[rp,]<-length(unique(thit))
    
    for (m in 1:length(rip.vec)) {
      
      rip<-rip.vec[m]
      ck<-ck.vec[m]
      
      temp<-custom.filter(comb_fl[rip],comb_fl[ck])
      
      hmmstat[rp,m]<-wilcox.test(temp[,1],temp[,2],alternative = 'greater',paired=T)$p.value
      
      
    }
  } else {
    hmmstat[rp,]<-NA
    hmmcount[rp,]<-0
    hmmtcount[rp,]<-0
    hmmmedlen[rp,]<-NA
  }
  
  
  
  
}

write.csv(hmmstat,'XNM_hmseekrhits_paired_unspliced_pval_CKnormed_p01.csv')
write.csv(hmmcount,'XNM_hmseekrhits_paired_unspliced_ncount_CKnormed_p01.csv')
write.csv(hmmtcount,'XNM_hmseekrhits_paired_unspliced_transcript_count_CKnormed_p01.csv')
write.csv(hmmmedlen,'XNM_hmseekrhits_paired_unspliced_medlen_CKnormed_p01.csv')



###################################
# normalize bams with control p 0.001
###################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

bamcounts<-read.table('readcounts_pre_filter_eCLIP_11_7_24.txt',sep = ' ',header=F)

bamcounts$V1<-sub("_.*", "", bamcounts$V1)

colnames(bamcounts)<-c('Filename','Total_Counts')


hitfiles<-list.files(path = './multicov/', pattern = "^tsc_condE_newpool_unspliced_.*_hits_multicov\\.out$", full.names = FALSE)

ckfiles<-list.files(path = './multicov/', pattern = "^tsc_condE_newpool_unspliced_.*_ck_multicov\\.out$", full.names = FALSE)



custom.filter<-function(x,y){
  
  temp<-as.data.frame(cbind(x,y))
  return(temp)
}


hmmstat<-as.data.frame(matrix(nrow=51,ncol=139))

xheaders <- c('XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
              'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
              'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14') 


nheaders<-paste0('NEAT1chk',c(1:22))
mheaders<-paste0('MALAT1chk',c(1:9))

xnm<-c(xheaders,nheaders,mheaders)

rownames(hmmstat)<-xnm


rbps<-read.table('multicov_files_list_XNM.txt',sep='\n')

rbps<-as.vector(rbps$V1)

rbps<-sub("_.*", "", rbps)

colnames(hmmstat)<-rbps[-12]

hmmcount<-hmmstat
hmmtcount<-hmmstat
hmmmedlen<-hmmstat

rip.vec<-rbps
ck.vec<-paste0(rbps,'_ck')

# remove control
rip.vec<-rip.vec[-12]
ck.vec<-ck.vec[-12]

# make sure the order in rbps are the same as in bamcounts
bamcounts <- bamcounts[match(rbps, bamcounts$Filename), ]

# counts file corresponds to multicov output cols
sum(rbps!=bamcounts$Filename) # 0

for (n in 1:length(hitfiles)) {
  
  print(n)
  
  hfile<-hitfiles[n]
  cfile<-ckfiles[n]
  
  # check if file is empty before processing
  # if empty store NA in stats and continue to the next n in loop
  filepath <- paste0('./multicov/', hfile)
  if (file.info(filepath)$size == 0) {
    print(paste0(hfile,' is empty'))
    
    chunk<-gsub('_4_viterbi_seekr_hits_multicov.out','',hfile)
    chunk<-gsub('tsc_condE_newpool_unspliced_','',chunk)
    # which row to store the stat
    rp<-which(xnm==chunk)
    
    hmmstat[rp,]<-NA
    hmmcount[rp,]<-0
    hmmtcount[rp,]<-0
    hmmmedlen[rp,]<-NA
    
    next
  }
  
  
  hitdata <- read.table(filepath, sep='\t')
  ckdata<-read.table(paste0('./multicov/',cfile),sep='\t')
  
  colnames(hitdata)<-c('chrom','chromStart','chromEnd','name','score','strand',rbps)
  colnames(ckdata)<-paste0(c('chrom','chromStart','chromEnd','name','score','strand',rbps),'_ck')
  
  hitdata$tlen<-hitdata$chromEnd-hitdata$chromStart
  ckdata$tlen_ck<-ckdata$chromEnd_ck-ckdata$chromStart_ck
  
  print(paste0('hitdata and ckdata has the same nrow_',nrow(hitdata)==nrow(ckdata)))
  
  print(paste0('hitdata and ckdata has the different length_',
               sum(hitdata$tlen!=ckdata$tlen_ck)))
  
  # get rpm count
  hitnorm<- sweep(hitdata[,c(7:146)], MARGIN=2, STATS=as.numeric(bamcounts$Total_Counts), FUN="/")
  hitnorm<- hitnorm*1000000
  hitnorm<-cbind(hitdata[,c(1:6)],hitnorm)
  hitnorm$tlen<-hitdata$tlen
  
  cknorm<- sweep(ckdata[,c(7:146)], MARGIN=2, STATS=as.numeric(bamcounts$Total_Counts), FUN="/")
  cknorm<- cknorm*1000000
  cknorm<-cbind(ckdata[,c(1:6)],cknorm)
  cknorm$tlen_ck<-ckdata$tlen_ck
  
  # subtract control from all samples
  # if negative after subtraction, set to 0
  hitnorm[,7:146] <- apply(hitnorm[,7:146], 2, function(x) x - hitnorm$control)
  hitnorm[,7:146] <- as.data.frame(lapply(hitnorm[,7:146], function(x) ifelse(x < 0, 0, x)))
  # drop the control col
  hitnorm$control<-NULL
  
  cknorm[,7:146] <- apply(cknorm[,7:146], 2, function(x) x - cknorm$control_ck)
  cknorm[,7:146] <- as.data.frame(lapply(cknorm[,7:146], function(x) ifelse(x < 0, 0, x)))
  
  # drop the control col
  cknorm$control_ck<-NULL
  
  
  comb<-cbind(hitnorm,cknorm)
  
  
  combname<-gsub('_4_viterbi_seekr_hits_multicov.out','',hfile)
  
  combname<-paste0('./combined/',combname,'_combined_CKnormed_p001.csv')
  
  write.csv(comb,combname,row.names = F)
  
  #########################################
  # filter the comb data
  # only filter on the score not score_ck
  # only filter on seekr p val score of hits not the control
  comb_fl<-comb[which(comb$score<0.001),]
  
  #######################################
  
  
  chunk<-gsub('_4_viterbi_seekr_hits_multicov.out','',hfile)
  chunk<-gsub('tsc_condE_newpool_unspliced_','',chunk)
  # which row to store the stat
  rp<-which(xnm==chunk)
  
  if (nrow(comb_fl)>0) {
    
    hmmmedlen[rp,]<-median(comb_fl$tlen,na.rm=T)
    hmmcount[rp,]<-nrow(comb_fl)
    
    thit<-sub("_[0-9]+$", "", comb_fl$name)
    hmmtcount[rp,]<-length(unique(thit))
    
    for (m in 1:length(rip.vec)) {
      
      rip<-rip.vec[m]
      ck<-ck.vec[m]
      
      temp<-custom.filter(comb_fl[rip],comb_fl[ck])
      
      hmmstat[rp,m]<-wilcox.test(temp[,1],temp[,2],alternative = 'greater',paired=T)$p.value

      
    }
  } else {
    hmmstat[rp,]<-NA
    hmmcount[rp,]<-0
    hmmtcount[rp,]<-0
    hmmmedlen[rp,]<-NA
  }
  
  
}

write.csv(hmmstat,'XNM_hmseekrhits_paired_unspliced_pval_CKnormed_p001.csv')
write.csv(hmmcount,'XNM_hmseekrhits_paired_unspliced_ncount_CKnormed_p001.csv')
write.csv(hmmtcount,'XNM_hmseekrhits_paired_unspliced_transcript_count_CKnormed_p001.csv')
write.csv(hmmmedlen,'XNM_hmseekrhits_paired_unspliced_medlen_CKnormed_p001.csv')



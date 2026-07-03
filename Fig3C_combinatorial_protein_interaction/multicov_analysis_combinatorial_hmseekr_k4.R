# check if two or more RBP tends to bind together in hits rather than control

#########################################
# unspliced
##########################################

###################################
# XNM normalize bams with control
###################################
setwd("XNM/seqstosummary")
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

bamcounts<-read.table('readcounts_pre_filter_eCLIP_11_7_24.txt',sep = ' ',header=F)

bamcounts$V1<-sub("_.*", "", bamcounts$V1)

colnames(bamcounts)<-c('Filename','Total_Counts')


hitfiles<-list.files(path = './multicov/', pattern = "^tsc_condE_newpool_unspliced_.*_hits_multicov\\.out$", full.names = FALSE)

ckfiles<-list.files(path = './multicov/', pattern = "^tsc_condE_newpool_unspliced_.*_ck_multicov\\.out$", full.names = FALSE)

ck2files<-list.files(path = './multicov/', pattern = "^tsc_condE_newpool_unspliced_.*_ck2_multicov\\.out$", full.names = FALSE)

xheaders <- c('XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
              'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
              'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14') 

nheaders<-paste0('NEAT1chk',c(1:22))
mheaders<-paste0('MALAT1chk',c(1:9))

xnm<-c(xheaders,nheaders,mheaders)


rbps<-read.table('multicov_files_list_XNM.txt',sep='\n')

rbps<-as.vector(rbps$V1)

rbps<-sub("_.*", "", rbps)

rip.vec<-rbps
ck.vec<-paste0(rbps,'_ck')
ck2.vec<-paste0(rbps,'_ck2')

# remove control
rip.vec<-rip.vec[-12]
ck.vec<-ck.vec[-12]
ck2.vec<-ck2.vec[-12] 

# make sure the order in rbps are the same as in bamcounts
bamcounts <- bamcounts[match(rbps, bamcounts$Filename), ]

# counts file corresponds to multicov output cols
sum(rbps!=bamcounts$Filename) # 0

########################
# set parameters
# list all query files
hitfiles
# choose the target query 
# XrA XrE M1 M3 M6
# use XrE as exmaple
n<-6
(hfile<-hitfiles[n])

# select the RBPs that are expected to bind altogether to the query
# if only testing two RBPs, comment out m3 and all corresponding codes
m1<-which(rip.vec=='U2AF1')
m2<-which(rip.vec=='U2AF2')
# m3<-which(rip.vec=='TIA1')
#################################

# check the control files are correct
(cfile<-ckfiles[n])
(c2file<-ck2files[n])


filepath <- paste0('./multicov/', hfile)
hitdata <- read.table(filepath, sep='\t')
ckdata<-read.table(paste0('./multicov/',cfile),sep='\t')
ck2data<-read.table(paste0('./multicov/',c2file),sep='\t')

colnames(hitdata)<-c('chrom','chromStart','chromEnd','name','score','strand',rbps)
colnames(ckdata)<-paste0(c('chrom','chromStart','chromEnd','name','score','strand',rbps),'_ck')
colnames(ck2data)<-paste0(c('chrom','chromStart','chromEnd','name','score','strand',rbps),'_ck2')


hitdata$tlen<-hitdata$chromEnd-hitdata$chromStart
ckdata$tlen_ck<-ckdata$chromEnd_ck-ckdata$chromStart_ck
ck2data$tlen_ck<-ck2data$chromEnd_ck-ck2data$chromStart_ck

print(paste0('hitdata and ck2data has the same nrow_',nrow(hitdata)==nrow(ck2data)))

print(paste0('hitdata and ck2data has the different length_',
             sum(hitdata$tlen!=ck2data$tlen_ck)))

# get rpm count
hitnorm<- sweep(hitdata[,c(7:146)], MARGIN=2, STATS=as.numeric(bamcounts$Total_Counts), FUN="/")
hitnorm<- hitnorm*1000000
hitnorm<-cbind(hitdata[,c(1:6)],hitnorm)
hitnorm$tlen<-hitdata$tlen

cknorm<- sweep(ckdata[,c(7:146)], MARGIN=2, STATS=as.numeric(bamcounts$Total_Counts), FUN="/")
cknorm<- cknorm*1000000
cknorm<-cbind(ckdata[,c(1:6)],cknorm)
cknorm$tlen_ck<-ckdata$tlen_ck

ck2norm<- sweep(ck2data[,c(7:146)], MARGIN=2, STATS=as.numeric(bamcounts$Total_Counts), FUN="/")
ck2norm<- ck2norm*1000000
ck2norm<-cbind(ck2data[,c(1:6)],ck2norm)
ck2norm$tlen_ck<-ck2data$tlen_ck

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


ck2norm[,7:146] <- apply(ck2norm[,7:146], 2, function(x) x - ck2norm$control_ck2)
ck2norm[,7:146] <- as.data.frame(lapply(ck2norm[,7:146], function(x) ifelse(x < 0, 0, x)))
# drop the control col
ck2norm$control_ck2<-NULL



comb<-cbind(hitnorm,cknorm,ck2norm)


#########################################
# filter the comb data
# only filter on the score not score_ck
# only filter on seekr p val score of hits not the control

# set the p val to intended thresholds 0.05, 0.01 or 0.001
comb_fl<-comb[which(comb$score<0.05),]

#######################################

(nrow(comb_fl)>0) 
  

(rip_1<-rip.vec[m1])
(rip_2<-rip.vec[m2])
# comment this out if only has two RBP
(rip_3<-rip.vec[m3])

(ck_1<-ck.vec[m1])
(ck_2<-ck.vec[m2])
# comment this out if only has two RBP
(ck_3<-ck.vec[m3])

(ck2_1<-ck2.vec[m1])
(ck2_2<-ck2.vec[m2])
# comment this out if only has two RBP
(ck2_3<-ck2.vec[m3])

rip_comp_1<-comb_fl[rip_1]>comb_fl[ck_1]
ck_comp_1<-comb_fl[ck2_1]>comb_fl[ck_1]

rip_comp_2<-comb_fl[rip_2]>comb_fl[ck_2]
ck_comp_2<-comb_fl[ck2_2]>comb_fl[ck_2]

# comment this out if only has two RBP
rip_comp_3<-comb_fl[rip_3]>comb_fl[ck_3]
ck_comp_3<-comb_fl[ck2_3]>comb_fl[ck_3]

# comment this out if only has two RBP
rip_comp<-rip_comp_1 & rip_comp_2 & rip_comp_3
ck_comp<-ck_comp_1 & ck_comp_2 & ck_comp_3

# use this block of code for counting if only have two RBP
rip_comp<-rip_comp_1 & rip_comp_2
ck_comp<-ck_comp_1 & ck_comp_2


#c(exp_true, exp_false,ck_true, ck_false)
mat <- matrix(c(sum(rip_comp), sum(!rip_comp),
                sum(ck_comp), sum(!ck_comp)),
              nrow = 2,
              byrow = TRUE,
              dimnames = list(
                group = c("rip", "ck"),
                outcome = c("true", "false")
              ))

mat

(res<-fisher.test(mat, alternative = "greater"))

res$p.value


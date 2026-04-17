# intersect spliced + unspliced v47 expressed transcripts
# with semi-extractable list
setwd("XNM/semi_extractable")

library(Biostrings)

semi<-read.csv('semi_extractables.csv',header=F)

# remove the _RI by the end of some gene names
semi$V1<-gsub('_RI','',semi$V1)
# 1074

pool<-readDNAStringSet('v47_filtered_ot1fixed_combined.fa')
v47names<-names(pool)
# 15293

v47list<-strsplit(v47names,'_',fixed=T)

v47genes<-sapply(v47list,'[[',1)

length(intersect(v47genes,semi$V1))
# 254 unique genes
sum(v47genes %in% semi$V1)
# 305 spliced + unspliced

v47semi<-v47names[v47genes %in% semi$V1]
write.csv(v47semi,'v47_semi-extractable_genes.csv',row.names = F)

v47semiunsplice<-v47semi[grepl('unspliced',v47semi)]

# check spliced and unspliced
v47unsplice<-v47names[grepl('unspliced',v47names)]
# 6663

v47unsplice<-strsplit(v47unsplice,'_',fixed=T)

v47unsplicenames<-sapply(v47unsplice,'[[',1)

sum(v47unsplicenames %in% semi$V1)
# 69


###########
# read in the NEAT1 similar genes 0.05
neat1<-read.csv('../pattern/NEAT1_pattern_matched_p05.csv',header=T)

unsplicesemicount<-sum(v47unsplicenames %in% semi$V1)
unsplicerest<-length(v47unsplicenames)-unsplicesemicount

neat1semicount<-sum(neat1$seqName %in% v47semi)
neat1rest<-nrow(neat1)-neat1semicount


######### unspliced
neat1uns<-neat1$seqName[grepl('unspliced',neat1$seqName)]

(neat1uns_semicount<-sum(neat1uns %in% v47semiunsplice))
(neat1uns_rest<-length(neat1uns)-neat1uns_semicount)

mat<-matrix(c(neat1uns_semicount,neat1uns_rest,unsplicesemicount,unsplicerest),nrow=2,ncol=2,byrow = T)
mat

(ftest<-fisher.test(mat,alternative = 'greater'))
ftest$p.value


###########
# read in the NEAT1 similar genes 0.01
neat1<-read.csv('../pattern/NEAT1_pattern_matched_p01.csv',header=T)

unsplicesemicount<-sum(v47unsplicenames %in% semi$V1)
unsplicerest<-length(v47unsplicenames)-unsplicesemicount

neat1semicount<-sum(neat1$seqName %in% v47semi)
neat1rest<-nrow(neat1)-neat1semicount


######### unspliced 
neat1uns<-neat1$seqName[grepl('unspliced',neat1$seqName)]

(neat1uns_semicount<-sum(neat1uns %in% v47semiunsplice))
(neat1uns_rest<-length(neat1uns)-neat1uns_semicount)

mat<-matrix(c(neat1uns_semicount,neat1uns_rest,unsplicesemicount,unsplicerest),nrow=2,ncol=2,byrow = T)
mat

(ftest<-fisher.test(mat,alternative = 'greater'))
ftest$p.value


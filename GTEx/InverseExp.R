# analyze XN like genes' inverse expression
# nearby genes would inversly correlate with X/N-like RNA
setwd("/GTEx")

library(data.table)
library(Biostrings) # needed to be loaded first, before rtracklayer
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)

# merge DEGs from all tissue types
getdeg<-function(filename) {
  temp<-read.csv(filename,header=T)
  temp<-temp[which(temp$padj<0.05 & abs(temp$log2FoldChange)>=1),]
  if (nrow(temp)>0) {
    return(temp$X)
  }
}

filelist<-dir(pattern='*.genelist_SB.csv')
fulldeg<-vector(mode='character',length = 0L)

for (file in filelist) {
  temp<-getdeg(file)
  fulldeg<-unique(c(fulldeg,temp))
  rm(temp)
}


fullexp<-fread('DESeqNormCounts_sb.csv')

degexp<-fullexp[(fullexp$V1 %in% fulldeg),]
rm(fullexp)

geneID<-degexp$V1
degexp<-as.data.frame(degexp)
degexp<-degexp[,c(2:ncol(degexp))]
rownames(degexp)<-geneID

# get gene properties
gtf<-import('comprehensive_annotation_v47_unspliced_new_kcnq1ot1_1_22_25.gtf')

# get the transcript id, gene id, strand, start and end coordinates for each gene
genes <- gtf[ mcols(gtf)$type == "gene" ]
genes <- as.data.table(genes)
genes <- genes[,c(1:5,10:12),with=F]

deggenes<-genes[gene_id %in% geneID]
deggenes<-as.data.frame(deggenes)

# read in XN like genes
xlike<-read.csv('k6XIST_pattern_matched_p05.csv',header=T)
nlike<-read.csv('k6NEAT1_pattern_matched_p01.csv',header=T)

xtemp<-strsplit(xlike$seqName,'_',fixed=T)
ntemp<-strsplit(nlike$seqName,'_',fixed=T)

xlike$geneID<-sapply(xtemp,'[[',6)
nlike$geneID<-sapply(ntemp,'[[',6)

xgeneID<-unique(xlike$geneID)
ngeneID<-unique(nlike$geneID)

length(intersect(xgeneID,geneID)) 
length(intersect(ngeneID,geneID)) 

xgeneID<-xgeneID[(xgeneID %in% geneID)]
ngeneID<-ngeneID[(ngeneID %in% geneID)]


nearbycor<-function(geneIDvec,deggenes,degexp){
  
  comb<-data.frame(gene1=as.character(),
                   gene2=as.character(),
                   cor_r=as.numeric())
  
  for (n in 1:length(geneIDvec)) {
    
    if (n %% 10 ==0) {print(n)}
    
    g<-geneIDvec[n]
    
    # get chr
    chrtemp<-as.character(deggenes$seqnames[which(deggenes$gene_id==g)])
    
    geneIDtemp<-deggenes$gene_id[deggenes$seqnames==chrtemp]
    
    degexptemp<-degexp[(rownames(degexp) %in% geneIDtemp),]
    
    targetgene<-degexp[which(rownames(degexp)==g),]
    
    pcor<-cor(t(targetgene),t(degexptemp),method = 'pearson',use='pairwise.complete.obs')
    
    n <- length(targetgene)
    
    t_stat <- pcor * sqrt((n - 2) / (1 - pcor^2))
    pvals <- 2 * pt(-abs(t_stat), df = n - 2)
    adjp<-p.adjust(pvals,method='BH')
    
    keep<-(abs(pcor)>0.25 & adjp<0.05)
    
    pcorval<-pcor[keep]
    pcorname<-colnames(pcor)[keep]
    
    res<-data.frame(gene1=rep(g,times=length(pcorval)),
                    gene2=pcorname,
                    cor_r=as.numeric(pcorval))
    
    comb<-rbind(comb,res)
  }
  
  return(comb)
  
}


xlikecor<-nearbycor(xgeneID,deggenes,degexp)

g1prop<-left_join(xlikecor[,1,drop=F],deggenes,by=join_by(gene1==gene_id))
colnames(g1prop)<-c('g1_geneID','g1_chr','g1_start','g1_end','g1_width','g1_strand','g1_geneType','g1_geneName')

g2prop<-left_join(xlikecor[,2,drop=F],deggenes,by=join_by(gene2==gene_id))
colnames(g2prop)<-c('g2_geneID','g2_chr','g2_start','g2_end','g2_width','g2_strand','g2_geneType','g2_geneName')

# distance of gene on + strand is neg if g2 is on the left side of g1, pos otherwise
# gene on - strand, distance is neg if g2 is on the right side of g1, pos otherwise
cal_dist<-function(df) {
  for (n in 1:nrow(df)) {
    
    if(n %% 10000 ==0) {print(n)}
    
    if (df$g1_strand[n]=='+') {
      df$distance[n]<-df$g2_start[n]-df$g1_start[n]
    } else if (df$g1_strand[n]=='-'){
      df$distance[n]<-df$g1_end[n]-df$g2_end[n]
    } else {
      print('strand not standard')
    }
    
  }
  return(df)
}


if (sum(g1prop$g1_geneID != xlikecor$gene1)==0 & sum(g2prop$g2_geneID != xlikecor$gene2)==0) {
  
  xres<-cbind(g1prop,xlikecor[,3,drop=F],g2prop)
  
  xres$distance<-NA
  
  xres<-cal_dist(xres)
  # only keep distance <1Mb
  xres<-xres[which(abs(xres$distance)<=500000),]
  
} else {
  print('error')
}

write.csv(xres,'Xlike_highcor_0.25_500kb_sb.csv',row.names = F)


nlikecor<-nearbycor(ngeneID,deggenes,degexp)

g1prop<-left_join(nlikecor[,1,drop=F],deggenes,by=join_by(gene1==gene_id))
colnames(g1prop)<-c('g1_geneID','g1_chr','g1_start','g1_end','g1_width','g1_strand','g1_geneType','g1_geneName')

g2prop<-left_join(nlikecor[,2,drop=F],deggenes,by=join_by(gene2==gene_id))
colnames(g2prop)<-c('g2_geneID','g2_chr','g2_start','g2_end','g2_width','g2_strand','g2_geneType','g2_geneName')

if (sum(g1prop$g1_geneID != nlikecor$gene1)==0 & sum(g2prop$g2_geneID != nlikecor$gene2)==0) {
  
  nres<-cbind(g1prop,nlikecor[,3,drop=F],g2prop)
  
  nres$distance<-NA
  
  nres<-cal_dist(nres)
  # only keep distance <1Mb
  nres<-nres[which(abs(nres$distance)<=500000),]
  
} else {
  print('error')
}

write.csv(nres,'Nlike_highcor_0.25_500kb_sb.csv',row.names = F)



#########################################################
# follow up stats
library(data.table)
library(dplyr)
library(tidyr)


xlike<-read.csv('Xlike_highcor_0.25_500kb_sb.csv',header=T)
nlike<-read.csv('Nlike_highcor_0.25_500kb_sb.csv',header=T)

# merge with network r values
networkr<-read.csv('XN_rnetwork_stab_exp_pattern_semiext_k6.csv',header=T)

temp<-strsplit(networkr$transcript,'_',fixed=T)
networkr$geneID<-sapply(temp,'[[',6)

xnet<-networkr[which(networkr$Xsimi_p0.05==TRUE),]
# remove spliced XIST
xnet<-xnet[which(xnet$transcript!='XIST_chrX_73820649_73852723_-_ENSG00000229807.14_ENST00000429829.6'),]

sum(!xlike$g1_geneID %in% xnet$geneID)

xlike<- xlike %>%
  left_join(xnet[,c(1:3,15)],by=join_by(g1_geneID==geneID))

write.csv(xlike,'Xlike_highcor_0.25_500kb_sb_networkr.csv',row.names = F)

nnet<-networkr[which(networkr$Nsimi_p0.01==TRUE),]

sum(!nlike$g1_geneID %in% nnet$geneID)

nlike<- nlike %>%
  left_join(nnet[,c(1,4,5,15)],by=join_by(g1_geneID==geneID))

write.csv(nlike,'Nlike_highcor_0.25_500kb_sb_networkr.csv',row.names = F)


# for each g1 gene, calculate how many genes are pos r and neg r
xsum<-xlike %>%
  group_by(g1_geneID) %>%
  summarise(
    g1_geneName=first(g1_geneName),
    g1_geneType=first(g1_geneType),
    XIST_r=first(XIST_r),
    XIST_perct=first(XIST_perct),
    ttl_pos_r=sum(cor_r>0,na.rm=T)-1,
    # -1 as we have the gene itself with r=1 in there
    ttl_neg_r=sum(cor_r<0,na.rm=T)
  )

xsum<-as.data.frame(xsum)
xsum$cor_bias<-xsum$ttl_pos_r-xsum$ttl_neg_r

sum(xsum$cor_bias>=0)
sum(xsum$cor_bias<0)

# for each g1 gene, calculate how many genes are pos r and neg r
nsum<-nlike %>%
  group_by(g1_geneID) %>%
  summarise(
    g1_geneName=first(g1_geneName),
    g1_geneType=first(g1_geneType),
    NEAT1_r=first(NEAT1_r),
    NEAT1_perct=first(NEAT1_perct),
    ttl_pos_r=sum(cor_r>0,na.rm=T)-1,
    # -1 as we have the gene itself with r=1 in there
    ttl_neg_r=sum(cor_r<0,na.rm=T)
  )

nsum<-as.data.frame(nsum)
nsum$cor_bias<-nsum$ttl_pos_r-nsum$ttl_neg_r

write.csv(xsum,'Xlike_highcor_SB_genestats.csv',row.names = F)
write.csv(nsum,'Nlike_highcor_SB_genestats.csv',row.names = F)

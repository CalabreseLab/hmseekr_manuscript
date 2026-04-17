# quantify 32 protein RBP network similarities between all unspliced transcripts
setwd("eCLIP_peakbams")
#################################################
# Make histograms and csv table of all 19295 transcripts network Pearson r values relative to XN
# replace NAs in the edge weights with 0
#################################################
library(dplyr)
library(data.table)

dt<-fread('comb_allgenes_cor_XNpattern.csv',header=T)

gene_cols <- setdiff(names(dt), c('RBP1','RBP2'))

# convert non-id columns to double (in place)
dt[, (gene_cols) := lapply(.SD, as.double), .SDcols = gene_cols]



long_dt <- melt(
  dt,
  id.vars = c('RBP1','RBP2'),
  variable.name = "transcript",
  value.name = "Weight",
  variable.factor = FALSE
)


long_dt$Weight[is.na(long_dt$Weight)]<-0

rm(dt)
dt<-long_dt
rm(long_dt)

dt$RBP1<-gsub('_norm','',dt$RBP1)
dt$RBP2<-gsub('_norm','',dt$RBP2)

######## subset X and N RBP list for its own pearson
xlist<-read.csv('top5sigRBP_Xmodule.csv',header=T)
# 19
nlist<-read.csv('top5sigRBP_Nmodule.csv',header=T)
# 18

xdt<-dt[which(dt$RBP1 %in% xlist$x & dt$RBP2 %in% xlist$x),]
ndt<-dt[which(dt$RBP1 %in% nlist$x & dt$RBP2 %in% nlist$x),]

xist<-'XIST_chrX_73820649_73852723_-_ENSG00000229807.14_ENSG00000229807.14.unspliced'
neat1<-'NEAT1_chr11_65422278_65445540_+_ENSG00000245532.12_ENST00000501122.3.monoexonic.unspliced'


# loop thru all genes against XN
xn <- data.table(
  transcript = character(),
  xist_r = numeric(),
  neat1_r = numeric()
)

genelist<-unique(dt$transcript)

# Filter rows for each gene using data.table's syntax
txdt <- xdt[transcript == xist]
# 171
tndt <- ndt[transcript == neat1]
# 153

for (n in 1:length(genelist)) {
  
  if (n %% 100 ==0) {print(n)}
  
  g<-genelist[n]
  
  txdt2 <- xdt[transcript == g]
  tndt2 <- ndt[transcript == g]
  
  # Check if the RBP1 and RBP2 are the same
  if (all(txdt$RBP1 == txdt2$RBP1 & txdt$RBP2 == txdt2$RBP2) &
          all(tndt$RBP1 == tndt2$RBP1 & tndt$RBP2 == tndt2$RBP2)) {
    xcor <- cor(txdt$Weight, txdt2$Weight, method = "pearson")
    ncor <- cor(tndt$Weight, tndt2$Weight, method = "pearson")
  } else {
    print('RBP1 and RBP2 matching problem')
    break
  }
  
  
  
  temp <- data.table(
    transcript = g,
    xist_r = xcor,
    neat1_r = ncor
  )
  
  xn<-rbindlist(list(xn, temp))
  
}

# still NA exists, as there are genes with all NAs as the edge weights before
# list all NA genes
xn$transcript[is.na(xn$xist_r)]
# 518
xn$transcript[is.na(xn$neat1_r)]
# 563

xn$xist_perct<-NA
xn$neat1_perct<-NA

for (n in 1:nrow(xn)) {
  
  if (!is.na(xn$xist_r[n])){
    xn$xist_perct[n]<-sum(xn$xist_r[n]>=xn$xist_r,na.rm=T)*100/sum(!is.na(xn$xist_r))
  }
  
  if (!is.na(xn$neat1_r[n])){
    xn$neat1_perct[n]<-sum(xn$neat1_r[n]>=xn$neat1_r,na.rm=T)*100/sum(!is.na(xn$neat1_r))
  }
  
}

colnames(xn)<-c('transcript','XIST_r','NEAT1_r','XIST_perct','NEAT1_perct')
xn<-xn[,c(1,2,4,3,5)]
write.csv(xn,'XN_vs_all_pearson_r_network_NA0.csv',row.names = F)


###################################
# wilcox test for XN similar ones against all unspliced
setwd("sigRBP_network")

networkr<-read.csv('XN_vs_all_pearson_r_network_NA0.csv',header=T)

nlike<-read.csv('../pattern/NEAT1_pattern_matched_p01.csv',header=T)
xlike<-read.csv('../pattern/XIST_pattern_matched_p05.csv',header=T)



ngene<-nlike$seqName[grepl('unspliced',nlike$seqName)]
xgene<-xlike$seqName[grepl('unspliced',xlike$seqName)]
# removed XIST spliced and 
# NOPCHAP1_chr12_104986316_105074197_+_ENSG00000151131.11_ENST00000552951.7

# remove itself
ngene<-ngene[!grepl('NEAT1',ngene)]
xgene<-xgene[!grepl('XIST',xgene)]

ndist<-networkr$NEAT1_r[networkr$transcript %in% ngene]
xdist<-networkr$XIST_r[networkr$transcript %in% xgene]

wilcox.test(xdist,networkr$XIST_r,alternative='greater')$p.value
# 2.354049e-15
wilcox.test(ndist,networkr$NEAT1_r,alternative='greater')$p.value
# 1.080915e-48

xplot<-data.frame(rval=c(xdist,networkr$XIST_r),
                  grp=c(rep('XISTlike',times=length(xdist)),
                        rep('All',times=nrow(networkr))))

nplot<-data.frame(rval=c(ndist,networkr$NEAT1_r),
                  grp=c(rep('NEAT1like',times=length(ndist)),
                        rep('All',times=nrow(networkr))))


library(ggplot2) 

p<-ggplot() +
  geom_density(data = xplot, aes(rval,color=grp),size=2) +
  labs(x = "Network Pearson's r values",y = "Density") +
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  coord_cartesian(xlim=c(-0.25,0.65))+
  scale_x_continuous(breaks=seq(from=-0.25, to=0.65, by=0.25),labels = scales::number_format(accuracy = 0.01))+
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color="grey", linewidth = 0.5),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='top',
        legend.title = element_blank(),
        legend.text = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.x=element_text(size=26, angle = 45, vjust=0.5),
        axis.text.y=element_text(size=26))


ggsave("Xlikevsall_pearson_network_density.pdf", plot = p, width = 5.5, height = 5, units = "in")



p<-ggplot() +
  geom_density(data = nplot, aes(rval,color=grp),size=2) +
  labs(x = "Network Pearson's r values",y = "Density") +
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  coord_cartesian(xlim=c(-0.25,0.7))+
  scale_x_continuous(breaks=seq(from=-0.25, to=0.75, by=0.25),labels = scales::number_format(accuracy = 0.01))+
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color="grey", linewidth = 0.5),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='top',
        legend.title = element_blank(),
        legend.text = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.x=element_text(size=26, angle = 45, vjust=0.5),
        axis.text.y=element_text(size=26))


ggsave("Nlikevsall_pearson_network_density.pdf", plot = p, width = 5.5, height = 5, units = "in")




################################

setwd("/GTEx")

library(tximport)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(data.table)

#####################################
cts<-fread('GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.gct',skip = 2)
headers<-as.data.frame(colnames(cts))
anno<-read.csv('sampleAnnotation.csv',header=T)

colnames(headers)[1]<-'colName'
headers$Index<-c(1:nrow(headers))

headers<-headers[3:nrow(headers),]

newheaders<-headers %>%
  left_join(anno[,c(1,6,7,11)],by=join_by(colName==SAMPID))

sum(!newheaders$colName==headers$colName)

headers<-newheaders
rm(newheaders)

sbidx<-headers$Index[which(headers$SMTSD %in% c("Cells - Cultured fibroblasts",
                                                "Whole Blood", "Spleen",
                                                "Cells - EBV-transformed lymphocytes",
                                                "Lung","Liver"))]
sbidx<-c(1,2,sbidx)

sbidx<-sort(sbidx,decreasing = F)

sbheaders<-headers[(headers$Index %in% sbidx),]

sum(!sbheaders$Index==sbidx[3:2927])

write.csv(sbheaders,'smallbatch_metadata.csv',row.names = F)

cts<-as.data.frame(cts)

sb<-cts[,sbidx]
write.csv(sb,'smallbatch.csv',row.names = F)


###################### 
# import data

cts <- fread('smallbatch.csv')

cts<-as.data.frame(cts)
# save gene feature as anther df
genefeature<-cts[,c(1:2)]
cts<-cts[,c(3:ncol(cts))]

cts <-as.matrix(cts)
row.names(cts)<-genefeature$Name


# prepare meta data for DESeq2 object construction
# It is absolutely critical that the columns of the count matrix and 
# the rows of the sampleTable (information about samples) 
# are in the same order.
sampleTable<-read.csv('smallbatch_metadata.csv',header=T)

unique(sampleTable$SMTSD)
# check that everything matched
all(sampleTable$colName %in% colnames(cts))
all(sampleTable$colName==colnames(cts))

#remove weird characters
sampleTable$SMTSD<-gsub(' ','',sampleTable$SMTSD,fixed=T)
sampleTable$SMTSD<-gsub('-','_',sampleTable$SMTSD,fixed=T)
unique(sampleTable$SMTSD)

# set reference level
sampleTable$SMTSD<-factor(sampleTable$SMTSD)
sampleTable$SMTSD <- relevel(sampleTable$SMTSD, ref = "Cells_Culturedfibroblasts")



# create Deseq2 obj from matrix
# the log2 fold change and Wald test p value will be for the last variable in the design formula, 
# and if this is a factor, the comparison will be the last level of this variable over the reference level

sirna <- DESeqDataSetFromMatrix(countData = cts,
                                colData = sampleTable,
                                design = ~ SMTSD)

sirna <- estimateSizeFactors(sirna)

#####################################
# explore and visualization
#####################################

############## pre-filter##################

# we can remove the rows that have no or nearly no information about the amount of gene expression. 
# minimal-removing rows of the DESeqDataSet that have no counts, or only a single count across all samples
# or at least 1500 samples (smallest group size) with a count of x or higher
nrow(sirna) # 74628
keep <- rowSums(counts(sirna) >= 1) >= 1500
sum(keep) #30838
sirna <- sirna[keep,]
nrow(sirna) 

##################### vst and rlog#########################
# remove heterskedasticity
# remove the dependence of variance on the mean 

sirna.vst <- vst(sirna, blind = FALSE)

########################## PCA ######################

p<-plotPCA(sirna.vst, intgroup = c("SMTSD"))

pdf('PCA_SB.pdf',width=18,height=15)
print(p,heatmap_legend_side = "right",padding = unit(c(30, 30, 30, 30), "mm"))
dev.off()

rm(p)

#######################################
# differential expression analysis
#######################################
# use DESeq to calculate differential expression

sirna.test<-DESeq(sirna)
resultsNames(sirna.test)

print('finished preanalysis')
############### generate results for a specific genotype ############################

tps<-as.character(unique(sampleTable$SMTSD))
tps<-setdiff(tps,'Cells_Culturedfibroblasts')
for (tp in tps) {
  sirna.res<- results(sirna.test,contrast=c('SMTSD',tp,'Cells_Culturedfibroblasts'),alpha=0.05)
  write.csv(sirna.res,paste0(tp,'_0.05_allgenelist_SB.csv'))
}

print('start normcount')
# get normalized counts
normalcounts<-counts(sirna.test, normalized=TRUE)
write.csv(normalcounts,'DESeqNormCounts_sb.csv')


print('start saving')
save.image('SBbatch.RData')

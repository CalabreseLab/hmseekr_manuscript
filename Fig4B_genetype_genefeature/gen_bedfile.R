# generate bedfiles for protein coding, lncRNA and others from gtf

# setwd set working directory properly
setwd("XNM")
library(Biostrings) # needed to be loaded first, before rtracklayer
library(rtracklayer)
library(data.table) # better data handling

options(scipen=999)

gtf_file<-'comprehensive_annotation_v47_unspliced_new_kcnq1ot1_1_22_25.gtf'

# extract info from gtf
gtf<-import(gtf_file)

# get the transcript id, gene id, strand, start and end coordinates for each exon
genes <- gtf[ mcols(gtf)$type == "gene" ]

# exons<-as.data.frame(exons)
genes <- as.data.table(genes)

genes$new_gene_type<-ifelse(genes$gene_type %in% c("lncRNA", "protein_coding"), genes$gene_type ,"other")

cds<-genes[genes$new_gene_type=='protein_coding']
lnc<-genes[genes$new_gene_type=='lncRNA']
other<-genes[genes$new_gene_type=='other']

chromlist<-c(paste0('chr',c(1:22)),'chrX','chrY')

organize_bed<-function(dt) {
  
  dt<-dt[,c(1,2,3,4,5,10,12)]
  
  dt$name<-paste0(dt$gene_id,'_',dt$gene_name)
  
  dt<-dt[,c(1,2,3,8,4,5)]
  
  colnames(dt)<-c('chrom','chromStart','chromEnd','name','score','strand')
  
  dt<-dt[dt$chrom %in% chromlist]
  
  return(dt)
  
}

cds_bed<-organize_bed(cds)
lnc_bed<-organize_bed(lnc)
other_bed<-organize_bed(other)

fwrite(cds_bed, 'XNM/v47_gene_types/protein_coding.bed', 
       sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

fwrite(lnc_bed, 'XNM/v47_gene_types/lncRNA.bed', 
       sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

fwrite(other_bed, 'XNM/v47_gene_types/other.bed', 
       sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)


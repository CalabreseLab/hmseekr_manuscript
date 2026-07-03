# analyze gene feature results with Fisher exact test

setwd("XNM")

library(gtools)
#################################################
### p 0.05 ####
# get the list of files

filelist<-list.files(path = './k6realdata05', pattern = '.*\\.csv$', full.names = F)

# loop through each file
# get the median of the three shuffle data
# then do a Fisher exact test for each query check

getchunk<-function(df) {
  
  df$qchunk<-gsub('tsc_condE_newpool_unspliced_','',df$file)
  df$qchunk<-gsub('_6_viterbi_seekr.*','',df$qchunk)
  return(df)
}

if (exists('comb_long')) {rm(comb_long)}
if (exists('comb_wide')) {rm(comb_wide)}

for (file in filelist) {
  
  rd<-read.csv(paste0('./k6realdata05/',file),header=T)
  
  s1<-read.csv(paste0('./k6shuffledata1_05/',file),header=T)
  s2<-read.csv(paste0('./k6shuffledata2_05/',file),header=T)
  s3<-read.csv(paste0('./k6shuffledata3_05/',file),header=T)
  
  # order which chk1 chk2 ... then chk10
  rd<-rd[mixedorder(rd$file), ]
  s1<-s1[mixedorder(s1$file), ]
  s2<-s2[mixedorder(s2$file), ]
  s3<-s3[mixedorder(s3$file), ]
  
  rd<-getchunk(rd)
  s1<-getchunk(s1)
  s2<-getchunk(s2)
  s3<-getchunk(s3)
  # check names are aligned
  print(paste0('namecheck_',file,'_',sum(rd$qchunk!=s1$qchunk),sum(rd$qchunk!=s2$qchunk),sum(rd$qchunk!=s3$qchunk)))
  # check if counts are the same
  rdttl<-rd$overlap+rd$nonoverlap
  s1ttl<-s1$overlap+s1$nonoverlap
  s2ttl<-s2$overlap+s2$nonoverlap
  s3ttl<-s3$overlap+s3$nonoverlap
  print(paste0('sumcheck_',file,'_',sum(rdttl!=s1ttl),sum(rdttl!=s2ttl),sum(rdttl!=s3ttl)))
  
  scomb<-s1
  soverlap<-cbind(s1$overlap,s2$overlap,s3$overlap)
  scomb$overlap<-apply(soverlap,1,median)
  snonover<-cbind(s1$nonoverlap,s2$nonoverlap,s3$nonoverlap)
  scomb$nonoverlap<-apply(snonover,1,median)
  
  pcomb<-data.frame(qchunk=rd$qchunk,pval=NA)
  
  # test whether rd has a significantly higher ratio of overlap events over nonoverlap
  # compared to shuffle
  pcomb$pval <- sapply(seq_len(nrow(rd)), function(i) {
      # Create a 2x2 matrix with the counts from df1 (first row) and df2 (second row)
      tab <- matrix(c(rd$overlap[i], rd$nonoverlap[i],
                      scomb$overlap[i], scomb$nonoverlap[i]), 
                    nrow = 2, byrow = TRUE)
      # Perform Fisher's Exact Test
      fisher.test(tab, alternative = "greater")$p.value
    })
  
  
  gf<-gsub('_stats.csv','',file)
  pcomb$genefeature<-gf
  
  
  if (exists('comb_long')) {
    comb_long<-rbind(comb_long,pcomb)
  } else {
    comb_long<-pcomb
  }
  
  
  pcomb$genefeature<-NULL
  colnames(pcomb)[2]<-gf

  if (exists('comb_wide')) {
    comb_wide<-cbind(comb_wide,pcomb[2])
  } else {
    comb_wide<-pcomb
  }
  
}

write.csv(comb_long,'k6genefeature_pval0.05_long.csv',row.names = F)
write.csv(comb_wide,'k6genefeature_pval0.05_wide.csv',row.names = F)

# adj pval for each gene type across all queries

comb_wide_adj<-data.frame(qchunk=comb_wide$qchunk,
                          unsplice_3_UTR_adjp=p.adjust(comb_wide$unsplice_3_UTR,method = 'BH'),
                          unsplice_5_UTR_adjp=p.adjust(comb_wide$unsplice_5_UTR,method = 'BH'),
                          unsplice_CDS_regions_adjp=p.adjust(comb_wide$unsplice_CDS_regions,method = 'BH'),
                          unsplice_distal_introns_adjp=p.adjust(comb_wide$unsplice_distal_introns,method = 'BH'),
                          unsplice_exons_adjp=p.adjust(comb_wide$unsplice_exons,method = 'BH'),
                          unsplice_proximal_introns_adjp=p.adjust(comb_wide$unsplice_proximal_introns,method = 'BH'),
                          unsplice_splice_sites_3_adjp=p.adjust(comb_wide$unsplice_splice_sites_3,method = 'BH'),
                          unsplice_splice_sites_5_adjp=p.adjust(comb_wide$unsplice_splice_sites_5,method = 'BH'),
                          unsplice_TES_regions_adjp=p.adjust(comb_wide$unsplice_TES_regions,method = 'BH'),
                          unsplice_TSS_regions_adjp=p.adjust(comb_wide$unsplice_TSS_regions,method = 'BH'))

write.csv(comb_wide_adj,'k6genefeature_adjpval0.05_wide.csv',row.names = F)


##################################################plot heatmap#####
# only unsplice ones

setwd("XNM")
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)


comb_wide<-read.csv('k6genefeature_adjpval0.05_wide.csv',header=T)

features<-colnames(comb_wide)[-1]

colnames(comb_wide)<-gsub('unsplice_','',colnames(comb_wide))

colnames(comb_wide)<-gsub('_adjp','',colnames(comb_wide))

comb_wide<-comb_wide[,c("qchunk","TSS_regions","5_UTR","CDS_regions","exons",
                        "splice_sites_5","proximal_introns","distal_introns",
                        "splice_sites_3","3_UTR","TES_regions")]

xheaders <- c('XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
              'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
              'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14') 

nheaders<-paste0('NEAT1chk',c(1:22))
mheaders<-paste0('MALAT1chk',c(1:9))

xnm<-c(xheaders,nheaders,mheaders)

xnm_plot<- comb_wide %>% 
  mutate(qchunk=factor(qchunk,levels=xnm)) %>%
  arrange(qchunk)


xnm_plot$group<-gsub('chk.*|int.*|r.*','',xnm_plot$qchunk)
xnm_plot$group<-factor(xnm_plot$group,levels=c('XIST','NEAT1','MALAT1'))
gcat<-xnm_plot$group
xnm_plot$group<-NULL

xnm_plot$rname<-gsub('ISTint|IST|EAT1chk|ALAT1chk','',xnm_plot$qchunk)
rownames(xnm_plot)<-xnm_plot$rname

xnm_plot$qchunk<-NULL
xnm_plot$rname<-NULL
xnm_plot<-as.matrix(xnm_plot)



# Define color function for the heatmap
col_fun <- colorRamp2(
  breaks = c(0, 0.05, 1),
  colors = c("#dd1c77", "#FCFC03", "#77dd1c")
)


# Create heatmap
ht<-Heatmap(
  xnm_plot,
  name = "p-val\n",
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,  
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  column_names_side = "bottom",
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 30,fontface = 3),
  row_names_gp = gpar(fontsize = 30,fontface = 3),
  row_split = gcat,
  row_gap = unit(15, "mm"),
  rect_gp = gpar(col = "lightgrey", lwd = 1),
  row_title = NULL,
  heatmap_legend_param = list(
    at = c(0, 0.05, 1),              
    labels = c("0", "0.05", "1"),
    labels_gp = gpar(fontsize = 30),
    title_gp = gpar(fontsize = 30),
    legend_direction = "vertical",
    legend_height = unit(5, "cm"),
    grid_width = unit(1, "cm")
  )
)

pdf('k6genefeature_adjpval0.05_heatmap_unspliceXNM.pdf',width=10,height=22)
draw(ht,heatmap_legend_side = "right",padding = unit(c(30, 30, 30, 30), "mm"))
dev.off()

#################################################
### p 0.01 ####
# get the list of files

filelist<-list.files(path = './k6realdata01', pattern = '.*\\.csv$', full.names = F)

# loop through each file
# get the median of the three shuffle data
# then do a Fisher exact test for each query check

getchunk<-function(df) {
  
  df$qchunk<-gsub('tsc_condE_newpool_unspliced_','',df$file)
  df$qchunk<-gsub('_6_viterbi_seekr.*','',df$qchunk)
  return(df)
}

if (exists('comb_long')) {rm(comb_long)}
if (exists('comb_wide')) {rm(comb_wide)}

for (file in filelist) {
  
  rd<-read.csv(paste0('./k6realdata01/',file),header=T)
  
  s1<-read.csv(paste0('./k6shuffledata1_01/',file),header=T)
  s2<-read.csv(paste0('./k6shuffledata2_01/',file),header=T)
  s3<-read.csv(paste0('./k6shuffledata3_01/',file),header=T)
  
  # order which chk1 chk2 ... then chk10
  rd<-rd[mixedorder(rd$file), ]
  s1<-s1[mixedorder(s1$file), ]
  s2<-s2[mixedorder(s2$file), ]
  s3<-s3[mixedorder(s3$file), ]
  
  rd<-getchunk(rd)
  s1<-getchunk(s1)
  s2<-getchunk(s2)
  s3<-getchunk(s3)
  # check names are aligned
  print(paste0('namecheck_',file,'_',sum(rd$qchunk!=s1$qchunk),sum(rd$qchunk!=s2$qchunk),sum(rd$qchunk!=s3$qchunk)))
  # check if counts are the same
  rdttl<-rd$overlap+rd$nonoverlap
  s1ttl<-s1$overlap+s1$nonoverlap
  s2ttl<-s2$overlap+s2$nonoverlap
  s3ttl<-s3$overlap+s3$nonoverlap
  print(paste0('sumcheck_',file,'_',sum(rdttl!=s1ttl),sum(rdttl!=s2ttl),sum(rdttl!=s3ttl)))
  
  scomb<-s1
  soverlap<-cbind(s1$overlap,s2$overlap,s3$overlap)
  scomb$overlap<-apply(soverlap,1,median)
  snonover<-cbind(s1$nonoverlap,s2$nonoverlap,s3$nonoverlap)
  scomb$nonoverlap<-apply(snonover,1,median)
  
  pcomb<-data.frame(qchunk=rd$qchunk,pval=NA)
  
  # test whether rd has a significantly higher ratio of overlap events over nonoverlap
  # compared to shuffle
  pcomb$pval <- sapply(seq_len(nrow(rd)), function(i) {
    # Create a 2x2 matrix with the counts from df1 (first row) and df2 (second row)
    tab <- matrix(c(rd$overlap[i], rd$nonoverlap[i],
                    scomb$overlap[i], scomb$nonoverlap[i]), 
                  nrow = 2, byrow = TRUE)
    # Perform Fisher's Exact Test
    fisher.test(tab, alternative = "greater")$p.value
  })
  
  
  gf<-gsub('_stats.csv','',file)
  pcomb$genefeature<-gf
  
  
  if (exists('comb_long')) {
    comb_long<-rbind(comb_long,pcomb)
  } else {
    comb_long<-pcomb
  }
  
  
  pcomb$genefeature<-NULL
  colnames(pcomb)[2]<-gf
  
  if (exists('comb_wide')) {
    comb_wide<-cbind(comb_wide,pcomb[2])
  } else {
    comb_wide<-pcomb
  }
  
}

write.csv(comb_long,'k6genefeature_pval0.01_long.csv',row.names = F)
write.csv(comb_wide,'k6genefeature_pval0.01_wide.csv',row.names = F)

comb_wide_adj<-data.frame(qchunk=comb_wide$qchunk,
                          unsplice_3_UTR_adjp=p.adjust(comb_wide$unsplice_3_UTR,method = 'BH'),
                          unsplice_5_UTR_adjp=p.adjust(comb_wide$unsplice_5_UTR,method = 'BH'),
                          unsplice_CDS_regions_adjp=p.adjust(comb_wide$unsplice_CDS_regions,method = 'BH'),
                          unsplice_distal_introns_adjp=p.adjust(comb_wide$unsplice_distal_introns,method = 'BH'),
                          unsplice_exons_adjp=p.adjust(comb_wide$unsplice_exons,method = 'BH'),
                          unsplice_proximal_introns_adjp=p.adjust(comb_wide$unsplice_proximal_introns,method = 'BH'),
                          unsplice_splice_sites_3_adjp=p.adjust(comb_wide$unsplice_splice_sites_3,method = 'BH'),
                          unsplice_splice_sites_5_adjp=p.adjust(comb_wide$unsplice_splice_sites_5,method = 'BH'),
                          unsplice_TES_regions_adjp=p.adjust(comb_wide$unsplice_TES_regions,method = 'BH'),
                          unsplice_TSS_regions_adjp=p.adjust(comb_wide$unsplice_TSS_regions,method = 'BH'))

write.csv(comb_wide_adj,'k6genefeature_adjpval0.01_wide.csv',row.names = F)

##################################################plot heatmap#####
# only unsplice ones

setwd("XNM")
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)


comb_wide<-read.csv('k6genefeature_adjpval0.01_wide.csv',header=T)

features<-colnames(comb_wide)[-1]

colnames(comb_wide)<-gsub('unsplice_','',colnames(comb_wide))

colnames(comb_wide)<-gsub('_adjp','',colnames(comb_wide))

comb_wide<-comb_wide[,c("qchunk","TSS_regions","5_UTR","CDS_regions","exons",
                        "splice_sites_5","proximal_introns","distal_introns",
                        "splice_sites_3","3_UTR","TES_regions")]

xheaders <- c('XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
              'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
              'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14') 

nheaders<-paste0('NEAT1chk',c(1:22))
mheaders<-paste0('MALAT1chk',c(1:9))

xnm<-c(xheaders,nheaders,mheaders)

xnm_plot<- comb_wide %>% 
  mutate(qchunk=factor(qchunk,levels=xnm)) %>%
  arrange(qchunk)


xnm_plot$group<-gsub('chk.*|int.*|r.*','',xnm_plot$qchunk)
xnm_plot$group<-factor(xnm_plot$group,levels=c('XIST','NEAT1','MALAT1'))
gcat<-xnm_plot$group
xnm_plot$group<-NULL

xnm_plot$rname<-gsub('ISTint|IST|EAT1chk|ALAT1chk','',xnm_plot$qchunk)
rownames(xnm_plot)<-xnm_plot$rname

xnm_plot$qchunk<-NULL
xnm_plot$rname<-NULL
xnm_plot<-as.matrix(xnm_plot)



# Define color function for the heatmap
col_fun <- colorRamp2(
  breaks = c(0, 0.05, 1),
  colors = c("#dd1c77", "#FCFC03", "#77dd1c")
)


# Create heatmap
ht<-Heatmap(
  xnm_plot,
  name = "p-val\n",
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,  
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  column_names_side = "bottom",
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 30,fontface = 3),
  row_names_gp = gpar(fontsize = 30,fontface = 3),
  row_split = gcat,
  row_gap = unit(15, "mm"),
  rect_gp = gpar(col = "lightgrey", lwd = 1),
  row_title = NULL,
  heatmap_legend_param = list(
    at = c(0, 0.05, 1),              
    labels = c("0", "0.05", "1"),
    labels_gp = gpar(fontsize = 30),
    title_gp = gpar(fontsize = 30),
    legend_direction = "vertical",
    legend_height = unit(5, "cm"),
    grid_width = unit(1, "cm")
  )
)

pdf('k6genefeature_adjpval0.01_heatmap_unspliceXNM.pdf',width=10,height=22)
draw(ht,heatmap_legend_side = "right",padding = unit(c(30, 30, 30, 30), "mm"))
dev.off()



#################################################
### p 0.001 ####
# get the list of files

filelist<-list.files(path = './k6realdata001', pattern = '.*\\.csv$', full.names = F)

# loop through each file
# get the median of the three shuffle data
# then do a Fisher exact test for each query check

getchunk<-function(df) {
  
  df$qchunk<-gsub('tsc_condE_newpool_unspliced_','',df$file)
  df$qchunk<-gsub('_6_viterbi_seekr.*','',df$qchunk)
  return(df)
}

if (exists('comb_long')) {rm(comb_long)}
if (exists('comb_wide')) {rm(comb_wide)}

for (file in filelist) {
  
  rd<-read.csv(paste0('./k6realdata001/',file),header=T)
  
  s1<-read.csv(paste0('./k6shuffledata1_001/',file),header=T)
  s2<-read.csv(paste0('./k6shuffledata2_001/',file),header=T)
  s3<-read.csv(paste0('./k6shuffledata3_001/',file),header=T)
  
  # order which chk1 chk2 ... then chk10
  rd<-rd[mixedorder(rd$file), ]
  s1<-s1[mixedorder(s1$file), ]
  s2<-s2[mixedorder(s2$file), ]
  s3<-s3[mixedorder(s3$file), ]
  
  rd<-getchunk(rd)
  s1<-getchunk(s1)
  s2<-getchunk(s2)
  s3<-getchunk(s3)
  # check names are aligned
  print(paste0('namecheck_',file,'_',sum(rd$qchunk!=s1$qchunk),sum(rd$qchunk!=s2$qchunk),sum(rd$qchunk!=s3$qchunk)))
  # check if counts are the same
  rdttl<-rd$overlap+rd$nonoverlap
  s1ttl<-s1$overlap+s1$nonoverlap
  s2ttl<-s2$overlap+s2$nonoverlap
  s3ttl<-s3$overlap+s3$nonoverlap
  print(paste0('sumcheck_',file,'_',sum(rdttl!=s1ttl),sum(rdttl!=s2ttl),sum(rdttl!=s3ttl)))
  
  scomb<-s1
  soverlap<-cbind(s1$overlap,s2$overlap,s3$overlap)
  scomb$overlap<-apply(soverlap,1,median)
  snonover<-cbind(s1$nonoverlap,s2$nonoverlap,s3$nonoverlap)
  scomb$nonoverlap<-apply(snonover,1,median)
  
  pcomb<-data.frame(qchunk=rd$qchunk,pval=NA)
  
  # test whether rd has a significantly higher ratio of overlap events over nonoverlap
  # compared to shuffle
  pcomb$pval <- sapply(seq_len(nrow(rd)), function(i) {
    # Create a 2x2 matrix with the counts from df1 (first row) and df2 (second row)
    tab <- matrix(c(rd$overlap[i], rd$nonoverlap[i],
                    scomb$overlap[i], scomb$nonoverlap[i]), 
                  nrow = 2, byrow = TRUE)
    # Perform Fisher's Exact Test
    fisher.test(tab, alternative = "greater")$p.value
  })
  
  
  gf<-gsub('_stats.csv','',file)
  pcomb$genefeature<-gf
  
  
  if (exists('comb_long')) {
    comb_long<-rbind(comb_long,pcomb)
  } else {
    comb_long<-pcomb
  }
  
  
  pcomb$genefeature<-NULL
  colnames(pcomb)[2]<-gf
  
  if (exists('comb_wide')) {
    comb_wide<-cbind(comb_wide,pcomb[2])
  } else {
    comb_wide<-pcomb
  }
  
}

write.csv(comb_long,'k6genefeature_pval0.001_long.csv',row.names = F)
write.csv(comb_wide,'k6genefeature_pval0.001_wide.csv',row.names = F)

comb_wide_adj<-data.frame(qchunk=comb_wide$qchunk,
                          unsplice_3_UTR_adjp=p.adjust(comb_wide$unsplice_3_UTR,method = 'BH'),
                          unsplice_5_UTR_adjp=p.adjust(comb_wide$unsplice_5_UTR,method = 'BH'),
                          unsplice_CDS_regions_adjp=p.adjust(comb_wide$unsplice_CDS_regions,method = 'BH'),
                          unsplice_distal_introns_adjp=p.adjust(comb_wide$unsplice_distal_introns,method = 'BH'),
                          unsplice_exons_adjp=p.adjust(comb_wide$unsplice_exons,method = 'BH'),
                          unsplice_proximal_introns_adjp=p.adjust(comb_wide$unsplice_proximal_introns,method = 'BH'),
                          unsplice_splice_sites_3_adjp=p.adjust(comb_wide$unsplice_splice_sites_3,method = 'BH'),
                          unsplice_splice_sites_5_adjp=p.adjust(comb_wide$unsplice_splice_sites_5,method = 'BH'),
                          unsplice_TES_regions_adjp=p.adjust(comb_wide$unsplice_TES_regions,method = 'BH'),
                          unsplice_TSS_regions_adjp=p.adjust(comb_wide$unsplice_TSS_regions,method = 'BH'))

write.csv(comb_wide_adj,'k6genefeature_adjpval0.001_wide.csv',row.names = F)

##################################################plot heatmap#####
# only unsplice ones

setwd("XNM")
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)


comb_wide<-read.csv('k6genefeature_adjpval0.001_wide.csv',header=T)

features<-colnames(comb_wide)[-1]

colnames(comb_wide)<-gsub('unsplice_','',colnames(comb_wide))

colnames(comb_wide)<-gsub('_adjp','',colnames(comb_wide))

comb_wide<-comb_wide[,c("qchunk","TSS_regions","5_UTR","CDS_regions","exons",
                        "splice_sites_5","proximal_introns","distal_introns",
                        "splice_sites_3","3_UTR","TES_regions")]

xheaders <- c('XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
              'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
              'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14') 

nheaders<-paste0('NEAT1chk',c(1:22))
mheaders<-paste0('MALAT1chk',c(1:9))

xnm<-c(xheaders,nheaders,mheaders)


xnm_plot<- comb_wide %>% 
  mutate(qchunk=factor(qchunk,levels=xnm)) %>%
  arrange(qchunk)


xnm_plot$group<-gsub('chk.*|int.*|r.*','',xnm_plot$qchunk)
xnm_plot$group<-factor(xnm_plot$group,levels=c('XIST','NEAT1','MALAT1'))
gcat<-xnm_plot$group
xnm_plot$group<-NULL

xnm_plot$rname<-gsub('ISTint|IST|EAT1chk|ALAT1chk','',xnm_plot$qchunk)
rownames(xnm_plot)<-xnm_plot$rname

xnm_plot$qchunk<-NULL
xnm_plot$rname<-NULL
xnm_plot<-as.matrix(xnm_plot)



# Define color function for the heatmap
col_fun <- colorRamp2(
  breaks = c(0, 0.05, 1),
  colors = c("#dd1c77", "#FCFC03", "#77dd1c")
)


# Create heatmap
ht<-Heatmap(
  xnm_plot,
  name = "p-val\n",
  col = col_fun,
  na_col = "white",
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,  
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  column_names_side = "bottom",
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 30,fontface = 3),
  row_names_gp = gpar(fontsize = 30,fontface = 3),
  row_split = gcat,
  row_gap = unit(15, "mm"),
  rect_gp = gpar(col = "lightgrey", lwd = 1),
  row_title = NULL,
  heatmap_legend_param = list(
    at = c(0, 0.05, 1),              
    labels = c("0", "0.05", "1"),
    labels_gp = gpar(fontsize = 30),
    title_gp = gpar(fontsize = 30),
    legend_direction = "vertical",
    legend_height = unit(5, "cm"),
    grid_width = unit(1, "cm")
  )
)

pdf('genefeature_adjpval0.001_heatmap_unspliceXNM.pdf',width=10,height=22)
draw(ht,heatmap_legend_side = "right",padding = unit(c(30, 30, 30, 30), "mm"))
dev.off()



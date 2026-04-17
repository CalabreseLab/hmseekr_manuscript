setwd("regulatory")

library(Biostrings)
library(dplyr)
library(tidyr)

# calculate total length for activator and repressor for spliced and unspliced

fa <- readDNAStringSet("../activator_repressor_seqs.fa")

# 2) Extract headers
hdr <- names(fa)

# 3) Split on `_` or `(` or `)`
#    Note: () are regex metacharacters, so we escape them
parts <- strsplit(hdr, "[_()]+")

category <- sapply(parts, '[[', 1)
splice <- sapply(parts, '[[', 7)

# Sequence lengths
seqlen <- width(fa)

# Build a data.frame
df <- data.frame(
  category = category,
  splice   = splice,
  seqlen   = seqlen,
  stringsAsFactors = FALSE
)

# Sum sequence lengths for all category–splice combinations
summary_len <- aggregate(seqlen ~ category + splice, data = df, sum)

write.csv(summary_len,'total_seqlen.csv',row.names = F)

###################### process single files 0.05

filelist05<-dir(path='./p0.05/',pattern='*._wide.csv')

dnames<-gsub('_wide.csv','',filelist05)
dnames<-dnames[c(1:17,20,21,18,19,22,33,37:43,23:32,34:36,44:49,64:69,50,56:63,51:55)]

fstats<-as.data.frame(matrix(ncol=length(dnames),nrow=6))
colnames(fstats)<-dnames
rownames(fstats)<-c('p0.05_unspliced','p0.05_spliced',
                    'p0.01_unspliced','p0.01_spliced',
                    'p0.001_unspliced','p0.001_spliced')

for (n in 1:length(filelist05)) {
  
  temp<-read.csv(paste0('./p0.05/',filelist05[n]),header=T)
  
  colidx<-which(dnames==gsub('_wide.csv','',filelist05[n]))
  
  prop<-strsplit(temp$seqName,'_|\\(|\\)')
  
  print(sum(lengths(prop)!=8))
  
  temp$category<-sapply(prop,'[[',1)
  
  temp$splice<-sapply(prop,'[[',7)
  

  # make sure all combination of cate and splice exist in res even if some does not exist
  # fill in 0 for the nonexist ones
  tempstats<- temp %>%
    group_by(category, splice) %>%
    summarize(total_unique_coverage = sum(total_unique_coverage), .groups = "drop") %>%
    complete(category=c('activator','repressor'), 
             splice=c('spliced','unspliced'), 
             fill = list(total_unique_coverage = 0))
  
  tempstats<-as.data.frame(tempstats)
  
  tempmerge<-merge(summary_len,tempstats,by=c('category','splice'))
  
  tempmerge$nonsig_lensum<-tempmerge$seqlen-tempmerge$total_unique_coverage
  
  tempmerge$seqlen<-NULL
  
  colnames(tempmerge)[3]<-'hmseekrunique_lensum'
  
  tempmerge$hmseekrunique_lensum<-round(tempmerge$hmseekrunique_lensum/100,digits = 0)
  tempmerge$nonsig_lensum<-round(tempmerge$nonsig_lensum/100,digits = 0)
  # fisher for splice and unsplice 
  
  fsp<-fisher.test(tempmerge[which(tempmerge$splice=='spliced'),c(3,4)],alternative = "two.sided")
  
  funsp<-fisher.test(tempmerge[which(tempmerge$splice=='unspliced'),c(3,4)],alternative = "two.sided")
  
  fstats['p0.05_unspliced',colidx]<-funsp$p.value
  
  fstats['p0.05_spliced',colidx]<-fsp$p.value
  
  tempmerge$funcdomain<-dnames[colidx]
  
  tempmerge$pthreshold<-'p0.05'
  
  tempmerge<-tempmerge[order(tempmerge$splice,decreasing = T),]
  
  tempmerge$fisherp<-NA
  
  tempmerge$fisherp[which(tempmerge$splice=='unspliced')]<-funsp$p.value
  tempmerge$fisherp[which(tempmerge$splice=='spliced')]<-fsp$p.value
  
  tempmerge$fisherp_bias<-''
  
  if (funsp$p.value < 0.05) {
    if (funsp$estimate > 1) {
      tempmerge$fisherp_bias[which(tempmerge$splice=='unspliced')]<-'activator'
    } else {
      tempmerge$fisherp_bias[which(tempmerge$splice=='unspliced')]<-'repressor'
    }
  }
  
  if (fsp$p.value < 0.05) {
    if (fsp$estimate > 1) {
      tempmerge$fisherp_bias[which(tempmerge$splice=='spliced')]<-'activator'
    } else {
      tempmerge$fisherp_bias[which(tempmerge$splice=='spliced')]<-'repressor'
    }
  }
  
  # calculate sig perct and perct diff
  tempmerge$sigperct<-tempmerge$hmseekrunique_lensum*100/(tempmerge$hmseekrunique_lensum+tempmerge$nonsig_lensum)
  
  tempmerge$perct_diff[which(tempmerge$splice=='unspliced')]<-tempmerge$sigperct[which(tempmerge$splice=='unspliced' & tempmerge$category=='activator')] -
    tempmerge$sigperct[which(tempmerge$splice=='unspliced' & tempmerge$category=='repressor')]
  
  tempmerge$perct_diff[which(tempmerge$splice=='spliced')]<-tempmerge$sigperct[which(tempmerge$splice=='spliced' & tempmerge$category=='activator')] -
    tempmerge$sigperct[which(tempmerge$splice=='spliced' & tempmerge$category=='repressor')]
  
  if (exists('combraw')) {
    combraw<-rbind(combraw,tempmerge)
  } else {
    combraw<-tempmerge
  }
}



############ 0.01

filelist01<-dir(path='./p0.01/',pattern='*._wide.csv')

for (n in 1:length(filelist01)) {
  
  temp<-read.csv(paste0('./p0.01/',filelist01[n]),header=T)
  
  colidx<-which(dnames==gsub('_wide.csv','',filelist01[n]))
  
  prop<-strsplit(temp$seqName,'_|\\(|\\)')
  
  print(sum((lengths(prop)!=8)))
  
  temp$category<-sapply(prop,'[[',1)
  
  temp$splice<-sapply(prop,'[[',7)
  
  
  # make sure all combination of cate and splice exist in res even if some does not exist
  # fill in 0 for the nonexist ones
  tempstats<- temp %>%
    group_by(category, splice) %>%
    summarize(total_unique_coverage = sum(total_unique_coverage), .groups = "drop") %>%
    complete(category=c('activator','repressor'), 
             splice=c('spliced','unspliced'), 
             fill = list(total_unique_coverage = 0))
  
  tempstats<-as.data.frame(tempstats)
  
  tempmerge<-merge(summary_len,tempstats,by=c('category','splice'))
  
  tempmerge$nonsig_lensum<-tempmerge$seqlen-tempmerge$total_unique_coverage
  
  tempmerge$seqlen<-NULL
  
  colnames(tempmerge)[3]<-'hmseekrunique_lensum'
  # fisher for splice and unsplice 
  
  tempmerge$hmseekrunique_lensum<-round(tempmerge$hmseekrunique_lensum/100,digits = 0)
  tempmerge$nonsig_lensum<-round(tempmerge$nonsig_lensum/100,digits = 0)
  
  fsp<-fisher.test(tempmerge[which(tempmerge$splice=='spliced'),c(3,4)],alternative = "two.sided")
  
  funsp<-fisher.test(tempmerge[which(tempmerge$splice=='unspliced'),c(3,4)],alternative = "two.sided")
  
  fstats['p0.01_unspliced',colidx]<-funsp$p.value
  
  fstats['p0.01_spliced',colidx]<-fsp$p.value
  
  tempmerge$funcdomain<-dnames[colidx]
  
  tempmerge$pthreshold<-'p0.01'
  
  tempmerge<-tempmerge[order(tempmerge$splice,decreasing = T),]
  
  tempmerge$fisherp<-NA
  
  tempmerge$fisherp[which(tempmerge$splice=='unspliced')]<-funsp$p.value
  tempmerge$fisherp[which(tempmerge$splice=='spliced')]<-fsp$p.value
  
  tempmerge$fisherp_bias<-''
  
  if (funsp$p.value < 0.05) {
    if (funsp$estimate > 1) {
      tempmerge$fisherp_bias[which(tempmerge$splice=='unspliced')]<-'activator'
    } else {
      tempmerge$fisherp_bias[which(tempmerge$splice=='unspliced')]<-'repressor'
    }
  }
  
  if (fsp$p.value < 0.05) {
    if (fsp$estimate > 1) {
      tempmerge$fisherp_bias[which(tempmerge$splice=='spliced')]<-'activator'
    } else {
      tempmerge$fisherp_bias[which(tempmerge$splice=='spliced')]<-'repressor'
    }
  }
  
  # calculate sig perct and perct diff
  tempmerge$sigperct<-tempmerge$hmseekrunique_lensum*100/(tempmerge$hmseekrunique_lensum+tempmerge$nonsig_lensum)
  
  tempmerge$perct_diff[which(tempmerge$splice=='unspliced')]<-tempmerge$sigperct[which(tempmerge$splice=='unspliced' & tempmerge$category=='activator')] -
    tempmerge$sigperct[which(tempmerge$splice=='unspliced' & tempmerge$category=='repressor')]
  
  tempmerge$perct_diff[which(tempmerge$splice=='spliced')]<-tempmerge$sigperct[which(tempmerge$splice=='spliced' & tempmerge$category=='activator')] -
    tempmerge$sigperct[which(tempmerge$splice=='spliced' & tempmerge$category=='repressor')]
  
  
  
  combraw<-rbind(combraw,tempmerge)
  
}



############ 0.001

filelist001<-dir(path='./p0.001/',pattern='*._wide.csv')

for (n in 1:length(filelist001)) {
  
  temp<-read.csv(paste0('./p0.001/',filelist001[n]),header=T)
  
  colidx<-which(dnames==gsub('_wide.csv','',filelist001[n]))
  
  
  prop<-strsplit(temp$seqName,'_|\\(|\\)')
  
  print(sum((lengths(prop)!=8)))
  
  temp$category<-sapply(prop,'[[',1)
  
  temp$splice<-sapply(prop,'[[',7)
  
  
  # make sure all combination of cate and splice exist in res even if some does not exist
  # fill in 0 for the nonexist ones
  tempstats<- temp %>%
    group_by(category, splice) %>%
    summarize(total_unique_coverage = sum(total_unique_coverage), .groups = "drop") %>%
    complete(category=c('activator','repressor'), 
             splice=c('spliced','unspliced'), 
             fill = list(total_unique_coverage = 0))
  
  tempstats<-as.data.frame(tempstats)
  
  tempmerge<-merge(summary_len,tempstats,by=c('category','splice'))
  
  tempmerge$nonsig_lensum<-tempmerge$seqlen-tempmerge$total_unique_coverage
  
  tempmerge$seqlen<-NULL
  
  colnames(tempmerge)[3]<-'hmseekrunique_lensum'
  # fisher for splice and unsplice 
  
  tempmerge$hmseekrunique_lensum<-round(tempmerge$hmseekrunique_lensum/100,digits = 0)
  tempmerge$nonsig_lensum<-round(tempmerge$nonsig_lensum/100,digits = 0)
  
  fsp<-fisher.test(tempmerge[which(tempmerge$splice=='spliced'),c(3,4)],alternative = "two.sided")
  
  funsp<-fisher.test(tempmerge[which(tempmerge$splice=='unspliced'),c(3,4)],alternative = "two.sided")
  
  fstats['p0.001_unspliced',colidx]<-funsp$p.value
  
  fstats['p0.001_spliced',colidx]<-fsp$p.value
  
  tempmerge$funcdomain<-dnames[colidx]
  
  tempmerge$pthreshold<-'p0.001'
  
  tempmerge<-tempmerge[order(tempmerge$splice,decreasing = T),]
  
  tempmerge$fisherp<-NA
  
  tempmerge$fisherp[which(tempmerge$splice=='unspliced')]<-funsp$p.value
  tempmerge$fisherp[which(tempmerge$splice=='spliced')]<-fsp$p.value
  
  tempmerge$fisherp_bias<-''
  
  if (funsp$p.value < 0.05) {
    if (funsp$estimate > 1) {
      tempmerge$fisherp_bias[which(tempmerge$splice=='unspliced')]<-'activator'
    } else {
      tempmerge$fisherp_bias[which(tempmerge$splice=='unspliced')]<-'repressor'
    }
  }
  
  if (fsp$p.value < 0.05) {
    if (fsp$estimate > 1) {
      tempmerge$fisherp_bias[which(tempmerge$splice=='spliced')]<-'activator'
    } else {
      tempmerge$fisherp_bias[which(tempmerge$splice=='spliced')]<-'repressor'
    }
  }
  
  # calculate sig perct and perct diff
  tempmerge$sigperct<-tempmerge$hmseekrunique_lensum*100/(tempmerge$hmseekrunique_lensum+tempmerge$nonsig_lensum)
  
  tempmerge$perct_diff[which(tempmerge$splice=='unspliced')]<-tempmerge$sigperct[which(tempmerge$splice=='unspliced' & tempmerge$category=='activator')] -
    tempmerge$sigperct[which(tempmerge$splice=='unspliced' & tempmerge$category=='repressor')]
  
  tempmerge$perct_diff[which(tempmerge$splice=='spliced')]<-tempmerge$sigperct[which(tempmerge$splice=='spliced' & tempmerge$category=='activator')] -
    tempmerge$sigperct[which(tempmerge$splice=='spliced' & tempmerge$category=='repressor')]
  
  
  combraw<-rbind(combraw,tempmerge)
  
  
}


write.csv(fstats,'single_functionaldomain_fisher_allp_hb.csv')
write.csv(combraw,'single_functionaldomain_rawdata_allp_hb.csv',row.names = F)

#### organize to have the final table
# only keep the activator row (same as repressor for pval)
combraw<-combraw[which(combraw$category=='activator'),]
combraw$hmseekrunique_lensum<-NULL
combraw$nonsig_lensum<-NULL
combraw$sigperct<-NULL

comb_wide <- combraw %>%
  pivot_wider(
    id_cols = c(funcdomain, splice),
    names_from = pthreshold,
    values_from = c(fisherp, fisherp_bias, perct_diff),
    names_glue = "{.value}_{pthreshold}",
    values_fill = list(
      fisherp = NA,
      fisherp_bias = "",
      perct_diff = NA
      
    )
  )


write.csv(comb_wide,'single_functionaldomain_orgdata_allp_hb.csv',row.names = F)

###################################
# plot order perct diff and plot barplot
# only plot top10 and bottom 10

library(ggplot2)
library(dplyr)
library(tidytext)
library(stringr)


# uncomment the lines of p0.01 and p0.001 for other plots
plotdf<-comb_wide[,c("funcdomain","splice","perct_diff_p0.05")]
# plotdf<-comb_wide[,c("funcdomain","splice","perct_diff_p0.01")]
# plotdf<-comb_wide[,c("funcdomain","splice","perct_diff_p0.001")]

colnames(plotdf)[3]<-'perct_diff'
plotdf$perct_diff[is.na(plotdf$perct_diff)]<-0

n_keep <- 10

plotdf <- plotdf %>%
  group_by(splice) %>%
  # take bottom 10 and top 10 within each splice
  slice_min(perct_diff, n = n_keep, with_ties = FALSE) %>%
  bind_rows(
    plotdf %>%
      group_by(splice) %>%
      slice_max(perct_diff, n = n_keep, with_ties = FALSE)
  ) %>%
  ungroup() %>%
  mutate(
    splice = factor(splice, levels = c("unspliced","spliced")),
    name_ord = reorder_within(funcdomain, perct_diff, splice),
    sign = if_else(perct_diff >= 0, "pos", "neg")
  )

range(plotdf$perct_diff)


p<-ggplot(plotdf,aes(x = perct_diff, y = name_ord, fill = sign)) +
  geom_col(width = 0.8) +
  geom_vline(xintercept = 0) +
  facet_wrap(.~ splice, scales = "free_y") +
  scale_y_reordered() +
  scale_x_continuous(limits = c(-21, 10), breaks = seq(-20, 10, 5)) +
  scale_fill_manual(values = c(pos = "red", neg = "blue"), guide = "none") +
  labs(x = "p 0.05 pert diff (act - rep)", y = NULL) +
  theme(
    panel.background=element_rect(fill='white'),
    panel.grid.major=element_blank(),
    axis.line.x = element_line(color="black", linewidth = 0.5),
    axis.line.y = element_line(color="black", linewidth = 0.5),
    legend.position='none',
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x=element_text(size=20, angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_text(hjust = 1,size = 20)
  )

ggsave("perct_diff_p0.05_top10.pdf", plot = p, width = 9, height = 8, units = "in")

###################################
# plot order GC content for all top10 btm10 perct diff
setwd("regulatory")

library(ggplot2)
library(dplyr)
library(tidytext)
library(stringr)
library(Biostrings)

comb_wide<-read.csv('single_functionaldomain_orgdata_allp_hb.csv',header=T)

seqs<-DNAStringSet(comb_wide$Sequence)
comb_wide$GCperct<-rowSums(letterFrequency(seqs, letters = c("G", "C"), as.prob = TRUE))*100
write.csv(comb_wide,'single_functionaldomain_orgdata_allp_hb.csv',row.names = F)

# uncomment the lines of p0.01 and p0.001 for other plots
plotdf<-comb_wide[,c("funcdomain","splice","perct_diff_p0.05","GCperct")]
# plotdf<-comb_wide[,c("funcdomain","splice","perct_diff_p0.01","GCperct")]
# plotdf<-comb_wide[,c("funcdomain","splice","perct_diff_p0.001","GCperct")]

colnames(plotdf)[3]<-'perct_diff'
plotdf$perct_diff[is.na(plotdf$perct_diff)]<-0

n_keep <- 10

plotdf <- plotdf %>%
  group_by(splice) %>%
  # take bottom 10 and top 10 within each splice
  slice_min(perct_diff, n = n_keep, with_ties = FALSE) %>%
  bind_rows(
    plotdf %>%
      group_by(splice) %>%
      slice_max(perct_diff, n = n_keep, with_ties = FALSE)
  ) %>%
  ungroup() %>%
  mutate(
    splice = factor(splice, levels = c("unspliced","spliced")),
    name_ord = reorder_within(funcdomain, perct_diff, splice),
    sign = if_else(perct_diff >= 0, "pos", "neg")
  )

range(plotdf$GCperct)

# calculate GC perct and perct_diff Pearson's correlation
temp<-plotdf[which(plotdf$splice=='spliced'),]
#temp<-plotdf[which(plotdf$splice=='unspliced'),]

(tcor<-cor.test(temp$perct_diff,temp$GCperct,alternative='two.sided',method='pearson'))
tcor$estimate
tcor$p.value

p<-ggplot(plotdf,aes(x = GCperct, y = name_ord)) +
  geom_col(width = 0.8,fill='#4daf4a') +
  geom_vline(xintercept = 0) +
  facet_wrap(.~ splice, scales = "free_y") +
  scale_y_reordered() +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_fill_manual(values = '#4daf4a', guide = "none") +
  labs(x = "GC percentage", y = NULL) +
  theme(
    panel.background=element_rect(fill='white'),
    panel.grid.major=element_blank(),
    axis.line.x = element_line(color="black", linewidth = 0.5),
    axis.line.y = element_line(color="black", linewidth = 0.5),
    legend.position='none',
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x=element_text(size=20, angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_text(hjust = 1,size = 20)
  )

ggsave("GCperct_p0.05_top10.pdf", plot = p, width = 9, height = 8, units = "in")



######################################
# get spreadsheet of rows as all transcripts in the pool
# column as each functional domain + total seqlen
# fill in the unique coverage for each search
# create one table for all three p vals

setwd("regulatory")

library(Biostrings)
library(dplyr)
library(tidyr)

# calculate total length for activator and repressor for spliced and unspliced

fa <- readDNAStringSet("../activator_repressor_seqs.fa")

# 2) Extract headers
seqNames <- names(fa)

# Sequence lengths
total_seqlen <- width(fa)

# Build a data.frame
comb<-as.data.frame(cbind(seqNames,total_seqlen))


filelist05<-dir(path='./p0.05/',pattern='*._wide.csv')

dnames<-gsub('_wide.csv','',filelist05)
dnames<-dnames[c(1:17,20,21,44:48,18,19,22,33,37:43,23:32,34:36,49,64:69,50,56:63,51:55)]

for (n in 1:length(dnames)) {
  
  temp<-read.csv(paste0('./p0.05/',dnames[n],'_wide.csv'),header=T)
  
  temp<-temp[,c(1,6)]
  colnames(temp)[2]<-paste0(dnames[n],'_tucp0.05')
  
  comb<-merge(comb,temp,by.x='seqNames',by.y='seqName',all.x=T)
  
}


filelist01<-dir(path='./p0.01/',pattern='*._wide.csv')

dnames<-gsub('_wide.csv','',filelist01)
dnames<-dnames[c(1:15,18,19,41:44,16:17,20,30,34:40,21:29,31:33,45,60:65,46,52:59,47:51)]

for (n in 1:length(dnames)) {
  
  temp<-read.csv(paste0('./p0.01/',dnames[n],'_wide.csv'),header=T)
  
  temp<-temp[,c(1,6)]
  colnames(temp)[2]<-paste0(dnames[n],'_tucp0.01')
  
  comb<-merge(comb,temp,by.x='seqNames',by.y='seqName',all.x=T)
  
}


filelist001<-dir(path='./p0.001/',pattern='*._wide.csv')

dnames<-gsub('_wide.csv','',filelist001)
dnames<-dnames[c(1:9,20,10,16:19,11:15,21,35:40,22,28:34,23:27)]

for (n in 1:length(dnames)) {
  
  temp<-read.csv(paste0('./p0.001/',dnames[n],'_wide.csv'),header=T)
  
  temp<-temp[,c(1,6)]
  colnames(temp)[2]<-paste0(dnames[n],'_tucp0.001')
  
  comb<-merge(comb,temp,by.x='seqNames',by.y='seqName',all.x=T)
  
}


write.csv(comb,'total_unique_coverage_table.csv',row.names = F)

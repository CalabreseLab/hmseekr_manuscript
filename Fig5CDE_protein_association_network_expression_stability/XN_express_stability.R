#############################################
# check relative stability of X and N like transcripts against all unspliced
setwd("sigRBP_network")

library(dplyr)

unspliced<-read.csv('XN_vs_all_pearson_r_network_NA0.csv',header=T)

nlike<-read.csv('../pattern/NEAT1_pattern_matched_p01.csv',header=T)
xlike<-read.csv('../pattern/XIST_pattern_matched_p05.csv',header=T)

stab<-read.table('combined_tpm_table_v47_biotype_length_stability_fraction_chromatin_association_01_22_25.tsv',header=T,sep='\t')

# use all unspliced as the background + XIST spliced 
bkglist<-c(unspliced$transcript,
           "XIST_chrX_73820649_73852723_-_ENSG00000229807.14_ENST00000429829.6")
# 6664
# XIST spliced 

# only keep the bkglist of transcript
stab<-stab[(stab$transcript_name %in% bkglist),]
# 6664
sum(stab$K_bru0_1+stab$K_bru0_2<=0)
# 62

stab$stability<-(stab$K_bru6_1+stab$K_bru6_2)/(stab$K_bru0_1+stab$K_bru0_2+0.125)
stab$expression<-(stab$K_tot_1+stab$K_tot_2)/2
# 6664 genes

stab<-stab[,c('transcript_name','stability','expression')]


ngene<-data.frame(transcript_name=nlike$seqName[grepl('unspliced',nlike$seqName)]) # 917
xgene<-data.frame(transcript_name=xlike$seqName[grepl('unspliced|XIST',xlike$seqName)]) # 87


ngenestab<-merge(ngene,stab,by.x = 'transcript_name',by.y = 'transcript_name',all.x=T)
xgenestab<-merge(xgene,stab,by.x = 'transcript_name',by.y = 'transcript_name',all.x=T)

ngenestab$grp='NEAT1like'
xgenestab$grp='XISTlike'
stab$grp='All'

sum(stab$stability==0) # 68
sum(xgenestab$stability==0) # 0 
sum(ngenestab$stability==0) # 0


stabplot<-rbind(xgenestab,ngenestab,stab)
# 7668

xnlines<-stab[grepl('NEAT1|XIST',stab$transcript_name),]
xnlines<-xnlines[c(1,3),]
xnlines$geneName<-c('NEAT1','XIST')

range(stabplot$stability)


median(ngenestab$stability)
median(xgenestab$stability)
median(stab$stability)

library(ggplot2) 
library(ggrepel)

p<-ggplot() +
  geom_density(data = stabplot, aes(stability,color=grp),
               linewidth=1.5,n=10000) +
  labs(x = "Bruseq stability",y = "Density") +
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  coord_cartesian(xlim=c(0,2))+
  scale_x_continuous(breaks=seq(from=0, to=2, by=0.5))+
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
        axis.text.y=element_text(size=26))+
  geom_vline(data = xnlines, aes(xintercept = stability), color = "black", linetype = "solid",linewidth=1.0)+
  geom_text_repel(data = xnlines, aes(x = stability, y = c(2.2,1.5), label = geneName), 
                  size = 8, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -0.3,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=3)


ggsave("XNlikevsunspliced_Bruseqstab.pdf", plot = p, width = 5.5, height = 5, units = "in")


range(stabplot$expression)

median(ngenestab$expression)
median(xgenestab$expression)
median(stab$expression)

p<-ggplot() +
  geom_density(data = stabplot, aes(expression,color=grp),
               linewidth=1.5,n=50000) +
  labs(x = "Total Expression",y = "Density") +
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  coord_cartesian(xlim=c(0,2.1))+
  scale_x_continuous(breaks=seq(from=0, to=2, by=0.25))+
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
        axis.text.y=element_text(size=26))+
  geom_vline(data = xnlines, aes(xintercept = expression), color = "black", linetype = "solid",linewidth=1.0)+
  geom_text_repel(data = xnlines, aes(x = expression, y = c(2,4), label = geneName), 
                  size = 8, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -0.3,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=3)


ggsave("XNlikevsunspliced_expression.pdf", plot = p, width = 5.5, height = 5, units = "in")

# stats of stability
wilcox.test(xgenestab$stability,stab$stability,alternative='greater')$p.value
wilcox.test(ngenestab$stability,stab$stability,alternative='greater')$p.value
wilcox.test(xgenestab$stability,ngenestab$stability,alternative='greater')$p.value

# stats of expression
wilcox.test(xgenestab$expression,stab$expression,alternative='greater')$p.value
wilcox.test(ngenestab$expression,stab$expression,alternative='greater')$p.value
wilcox.test(xgenestab$expression,ngenestab$expression,alternative='greater')$p.value


wilcox.test(xgenestab$expression,stab$expression,alternative='less')$p.value
wilcox.test(ngenestab$expression,stab$expression,alternative='less')$p.value
wilcox.test(xgenestab$expression,ngenestab$expression,alternative='less')$p.value

#############################################
# check relative stability of X and N like transcripts against all unspliced + spliced
setwd("sigRBP_network")

library(dplyr)
library(Biostrings)

bkg<-readDNAStringSet('v47_filtered_ot1fixed_combined.fa')
bkglist<-names(bkg)
# 15293

nlike<-read.csv('../pattern/NEAT1_pattern_matched_p01.csv',header=T)
xlike<-read.csv('../pattern/XIST_pattern_matched_p05.csv',header=T)

stab<-read.table('combined_tpm_table_v47_biotype_length_stability_fraction_chromatin_association_01_22_25.tsv',header=T,sep='\t')


# only keep the bkglist of transcript
stab<-stab[(stab$transcript_name %in% bkglist),]
# 15293
sum(stab$K_bru0_1+stab$K_bru0_2<=0)
# 1647

stab$stability<-(stab$K_bru6_1+stab$K_bru6_2)/(stab$K_bru0_1+stab$K_bru0_2+0.125)
stab$expression<-(stab$K_tot_1+stab$K_tot_2)/2


stab<-stab[,c('transcript_name','stability','expression')]


ngene<-data.frame(transcript_name=nlike$seqName) # 918
xgene<-data.frame(transcript_name=xlike$seqName) # 88


ngenestab<-merge(ngene,stab,by.x = 'transcript_name',by.y = 'transcript_name',all.x=T)
xgenestab<-merge(xgene,stab,by.x = 'transcript_name',by.y = 'transcript_name',all.x=T)

ngenestab$grp='NEAT1like'
xgenestab$grp='XISTlike'
stab$grp='All'

sum(stab$stability==0) # 1093
sum(xgenestab$stability==0) # 0 
sum(ngenestab$stability==0) # 0


stabplot<-rbind(xgenestab,ngenestab,stab)
# 16299

xnlines<-stab[grepl('NEAT1|XIST',stab$transcript_name),]
xnlines<-xnlines[c(1,3),]
xnlines$geneName<-c('NEAT1','XIST')

range(stabplot$stability)


median(ngenestab$stability)
median(xgenestab$stability)
median(stab$stability)

library(ggplot2) 
library(ggrepel)

p<-ggplot() +
  geom_density(data = stabplot, aes(stability,color=grp),
               linewidth=1.5,n=10000) +
  labs(x = "Bruseq stability",y = "Density") +
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  coord_cartesian(xlim=c(0,5))+
  scale_x_continuous(breaks=seq(from=0, to=5, by=1))+
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
        axis.text.y=element_text(size=26))+
  geom_vline(data = xnlines, aes(xintercept = stability), color = "black", linetype = "solid",linewidth=1.0)+
  geom_text_repel(data = xnlines, aes(x = stability, y = c(2.2,1.5), label = geneName), 
                  size = 8, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -0.3,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=3)


ggsave("XNlikevsALL_Bruseqstab.pdf", plot = p, width = 5.5, height = 5, units = "in")


range(stabplot$expression)

median(ngenestab$expression)
median(xgenestab$expression)
median(stab$expression)

p<-ggplot() +
  geom_density(data = stabplot, aes(expression,color=grp),
               linewidth=1.5,n=50000) +
  labs(x = "Total Expression",y = "Density") +
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  coord_cartesian(xlim=c(0,3))+
  scale_x_continuous(breaks=seq(from=0, to=3, by=1))+
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
        axis.text.y=element_text(size=26))+
  geom_vline(data = xnlines, aes(xintercept = expression), color = "black", linetype = "solid",linewidth=1.0)+
  geom_text_repel(data = xnlines, aes(x = expression, y = c(2,4), label = geneName), 
                  size = 8, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -0.3,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=3)


ggsave("XNlikevsALL_expression.pdf", plot = p, width = 5.5, height = 5, units = "in")

# stats of stability
wilcox.test(xgenestab$stability,stab$stability,alternative='greater')$p.value
wilcox.test(ngenestab$stability,stab$stability,alternative='greater')$p.value
wilcox.test(xgenestab$stability,ngenestab$stability,alternative='greater')$p.value

wilcox.test(xgenestab$stability,stab$stability,alternative='less')$p.value
wilcox.test(ngenestab$stability,stab$stability,alternative='less')$p.value
wilcox.test(xgenestab$stability,ngenestab$stability,alternative='less')$p.value


# stats of expression
wilcox.test(xgenestab$expression,stab$expression,alternative='greater')$p.value
wilcox.test(ngenestab$expression,stab$expression,alternative='greater')$p.value
wilcox.test(xgenestab$expression,ngenestab$expression,alternative='greater')$p.value


wilcox.test(xgenestab$expression,stab$expression,alternative='less')$p.value
wilcox.test(ngenestab$expression,stab$expression,alternative='less')$p.value
wilcox.test(xgenestab$expression,ngenestab$expression,alternative='less')$p.value




#############################
# append XN like at 3 pvals and semi-extractable 

setwd("sigRBP_network")

comb<-read.csv('XN_vs_all_pearson_r_network_NA0.csv',header=T)

# add spliced XN like back in
spliced<-data.frame(transcript=c("XIST_chrX_73820649_73852723_-_ENSG00000229807.14_ENST00000429829.6",
                                 "NOPCHAP1_chr12_104986316_105074197_+_ENSG00000151131.11_ENST00000552951.7",
                                 "ATP11A_chr13_112690038_112887168_+_ENSG00000068650.19_ENST00000375645.8"),
                    XIST_r=c(NA,NA,NA),
                    NEAT1_r=c(NA,NA,NA))

comb<-rbind(comb,spliced)

stab$grp<-NULL

colnames(stab)<-c("transcript_name","Bruseq_stability","total_expression")

comb<-merge(comb,stab,by.x='transcript',by.y='transcript_name',all.x=T)

comb$Xsimi_p0.05<-'FALSE'
comb$Xsimi_p0.01<-'FALSE'
comb$Xsimi_p0.001<-'FALSE'

comb$Nsimi_p0.05<-'FALSE'
comb$Nsimi_p0.01<-'FALSE'
comb$Nsimi_p0.001<-'FALSE'

comb$semi_extract<-'FALSE'

xsimi<-read.csv('../pattern/XIST_pattern_matched_p05.csv',header=T)
comb$Xsimi_p0.05[comb$transcript %in% xsimi$seqName]<-'TRUE'

xsimi<-read.csv('../pattern/XIST_pattern_matched_p01.csv',header=T)
comb$Xsimi_p0.01[comb$transcript %in% xsimi$seqName]<-'TRUE'

xsimi<-read.csv('../pattern/XIST_pattern_matched_p001.csv',header=T)
comb$Xsimi_p0.001[comb$transcript %in% xsimi$seqName]<-'TRUE'

nsimi<-read.csv('../pattern/NEAT1_pattern_matched_p05.csv',header=T)
comb$Nsimi_p0.05[comb$transcript %in% nsimi$seqName]<-'TRUE'

nsimi<-read.csv('../pattern/NEAT1_pattern_matched_p01.csv',header=T)
comb$Nsimi_p0.01[comb$transcript %in% nsimi$seqName]<-'TRUE'

nsimi<-read.csv('../pattern/NEAT1_pattern_matched_p001.csv',header=T)
comb$Nsimi_p0.001[comb$transcript %in% nsimi$seqName]<-'TRUE'

semi<-read.csv('../semi_extractable/v47_semi-extractable_genes.csv',header=T)
comb$semi_extract[comb$transcript %in% semi$x]<-'TRUE'

write.csv(comb,'XN_rnetwork_stab_exp_pattern_semiext.csv',row.names = F)



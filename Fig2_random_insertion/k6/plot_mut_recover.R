# plot results
setwd("XNM/recovery/k6/wt")
setwd("XNM/recovery/k6/mut10")
setwd("XNM/recovery/k6/mut20")
setwd("XNM/recovery/k6/mut30")
setwd("XNM/recovery/k6/mut40")
setwd("XNM/recovery/k6/mut50")

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

# read in the files 112
files <- list.files(
  path        = ".", 
  pattern     = "_recover_percentage\\.csv$", 
  full.names  = FALSE
)


if (exists('comb')) {rm(comb)}

for (f in files) {
  
  temp<-fread(f,header=T)
  
  temp$pool_id<-NULL
  temp$insert_pos<-NULL
  
  if (exists('comb')) {
    comb<-rbindlist(list(comb,temp),fill=F)
  } else {
    comb<-temp
  }
  
}

rm(temp)

#sum(comb$recover_pert==0)
#comb[which(comb$recover_pert==0),]
#table(comb$query_id[which(comb$recover_pert==0)])

xheaders <- c('XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
              'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
              'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14') 

nheaders<-paste0('NEAT1chk',c(1:22))
mheaders<-paste0('MALAT1chk',c(1:9))

xnm<-c(xheaders,nheaders,mheaders)

# calculate mean for each query_id for each pval
# and total number that is >0 for each query_id and each pval

df_summary<- comb %>%
  group_by(query_id) %>%
  summarise(across(starts_with('recover_pert'),
            list(mean=~mean(.x, na.rm=TRUE),
                 count_gt0=~sum(.x>0,na.rm=T)),
            .names='{.col}_{.fn}')
  ) %>%
  mutate(query_id=factor(query_id,levels=xnm)) %>%
  arrange(query_id)


write.csv(df_summary,'../XNM_queryrecover_mean_mut50.csv',row.names = F)


mean(recov_mean$recover_pert)
# 99.87515
median(recov_mean$recover_pert)
# 99.97989
min(recov_mean$recover_pert)
# 97.72397


##################### plot
setwd("XNM/recovery/k6")
library(ggplot2)

wt<-read.csv('XNM_queryrecover_mean_wt.csv',header=T)
wt$group<-'WT'

mut10<-read.csv('XNM_queryrecover_mean_mut10.csv',header=T)
sum(wt$query_id!=mut10$query_id)
mut10$group<-'mut10%'

mut20<-read.csv('XNM_queryrecover_mean_mut20.csv',header=T)
sum(wt$query_id!=mut20$query_id)
mut20$group<-'mut20%'

mut30<-read.csv('XNM_queryrecover_mean_mut30.csv',header=T)
sum(wt$query_id!=mut30$query_id)
mut30$group<-'mut30%'

mut40<-read.csv('XNM_queryrecover_mean_mut40.csv',header=T)
sum(wt$query_id!=mut40$query_id)
mut40$group<-'mut40%'

mut50<-read.csv('XNM_queryrecover_mean_mut50.csv',header=T)
sum(wt$query_id!=mut50$query_id)
mut50$group<-'mut50%'

comb<-rbind(wt,mut10,mut20,mut30,mut40,mut50)

range(comb$recover_pert_mean)

comb$group<-factor(comb$group,levels=c('WT','mut10%','mut20%','mut30%','mut40%','mut50%'))


rand_summary<- comb %>%
  group_by(group) %>%
  summarise(across(starts_with('recover_pert'),median))

write.csv(rand_summary,'XNM_queryrecover_randomization_median_k6.csv',row.names = F)


combperct<-comb[,c(1,2,4,6,8,10)]
combcount<-comb[,c(1,3,5,7,9,10)]

colnames(combperct)[2:5]<-c('nofilter','p0.05','p0.01','p0.001')
colnames(combcount)[2:5]<-c('nofilter','p0.05','p0.01','p0.001')

combperct_l<-pivot_longer(combperct,-c(query_id,group),names_to = 'filter',values_to = 'recover_perct_mean')
combcount_l<-pivot_longer(combcount,-c(query_id,group),names_to = 'filter',values_to = 'recover_count_gt0')

combperct_l$group<-factor(combperct_l$group,levels=c('WT','mut10%','mut20%','mut30%','mut40%','mut50%'))
combperct_l$filter<-factor(combperct_l$filter,levels=c('nofilter','p0.05','p0.01','p0.001'))

combcount_l$group<-factor(combcount_l$group,levels=c('WT','mut10%','mut20%','mut30%','mut40%','mut50%'))
combcount_l$filter<-factor(combcount_l$filter,levels=c('nofilter','p0.05','p0.01','p0.001'))


p<-ggplot(data=combperct_l,aes(x=group,y=recover_perct_mean,fill=filter))+
  geom_boxplot(aes(group = interaction(group, filter)),         # key!
               position = position_dodge2(width = 0.75, preserve = "single"),
               outlier.shape = 19) +
  # geom_boxplot(color='#4a4a4a' ,outlier.shape = NA,fill='grey') +
  # geom_point(aes(color=group),size=4,shape=19,alpha=0.85,
  #            position = position_jitter(width = 0.35,height = 0))+
  coord_cartesian(ylim=c(0,100))+
  scale_y_continuous(breaks = seq(from=0, to=100, by=20))+
  scale_x_discrete()+
  theme(
    panel.background=element_rect(fill='white'),
    panel.grid.major=element_line(color='grey',linewidth=0.3),
    legend.key.height = unit(1,'line'),
    legend.key.width  = unit(1,'line'),
    legend.title = element_text(size=30),
    legend.text=element_text(size=30),
    legend.position = 'top',
    axis.text.x=element_text(size=30,color='black',angle=45,hjust=1),
    #axis.text.x=element_text(size=30,color='black'),
    axis.text.y=element_text(size=30,color='black'),
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30,angle=90),
    strip.text = element_text(size = 30)
  ) +
  labs(x = "Queries", y = "mean recover%")+
  scale_fill_brewer(palette='Set1')
#facet_wrap(~dataset, ncol = 4, scales = "fixed")

pdf('mean_queryrecover_perct_boxplot_allmuts_k6.pdf',width=13,height=8.5)
plot(p)
dev.off()



p<-ggplot(data=combcount_l,aes(x=group,y=recover_count_gt0/10,fill=filter))+
  geom_boxplot(aes(group = interaction(group, filter)),         # key!
               position = position_dodge2(width = 0.75, preserve = "single"),
               outlier.shape = 19) +
  # geom_boxplot(color='#4a4a4a' ,outlier.shape = NA,fill='grey') +
  # geom_point(aes(color=group),size=4,shape=19,alpha=0.85,
  #            position = position_jitter(width = 0.35,height = 0))+
  coord_cartesian(ylim=c(0,100))+
  scale_y_continuous(breaks = seq(from=0, to=100, by=20))+
  scale_x_discrete()+
  theme(
    panel.background=element_rect(fill='white'),
    panel.grid.major=element_line(color='grey',linewidth=0.3),
    legend.key.height = unit(1,'line'),
    legend.key.width  = unit(1,'line'),
    legend.title = element_text(size=30),
    legend.text=element_text(size=30),
    legend.position = 'top',
    axis.text.x=element_text(size=30,color='black',angle=45,hjust=1),
    #axis.text.x=element_text(size=30,color='black'),
    axis.text.y=element_text(size=30,color='black'),
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30,angle=90),
    strip.text = element_text(size = 30)
  ) +
  labs(x = "Queries", y = "recover% > 0 count")+
  scale_fill_brewer(palette='Set1')
#facet_wrap(~dataset, ncol = 4, scales = "fixed")

pdf('mean_queryrecover_count_boxplot_allmuts_k6.pdf',width=13,height=8.5)
plot(p)
dev.off()




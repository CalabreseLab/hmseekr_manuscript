#################################################
########## p0.05
setwd("XNM/pattern")

library(Biostrings)
library(dplyr)
library(tidyr)


###################### process XIST files
# transcripts that contained a hit to Repeats A and F within its first 3kb
# followed by hits to both Repeats B1/2 
# followed by hits to D, followed by a hit to Repeat E. At p < 0.05


#### first get all transcripts that contain hits to all features
filelist05<-dir(path='./XIST_05/',pattern='*._wide.csv')

for (n in 1:length(filelist05)) {
  
  temp<-read.csv(paste0('./XIST_05/',filelist05[n]),header=T)
  
  temp<-temp[,c("seqName","total_unique_coverage")]
  
  colnames(temp)[2]<-paste0(gsub('_wide.csv','',filelist05[n]),'_tuc0.05')
  
  if (exists('comb')) {
    # only keep transcripts that has hits in both dfs
    comb<-inner_join(comb, temp, by = "seqName")
  } else {
    comb<-temp
  }
}

write.csv(comb,'transcripts_withallXISThits_p05.csv',row.names = F)
# 1689
sum(grepl('unspliced',comb$seqName)) # 1684 only 5 spliced transcripts

### go through the raw file one by one and check the pattern

xra<-read.csv('./XIST_05/XISTrA_raw.csv',header=T)

# keep the rows/transcripts in xra that exists in comb
xra<-semi_join(xra,comb,by='seqName')

# check length 
length(unique(xra$seqName))

# filter for rows that has start <=3000
xra<-xra[which(xra$Start<=3000),]

xra<-xra[,c(1,4)]

xra_collapsed <- xra %>%
  group_by(seqName) %>%
  summarise(
    Start   = min(Start,   na.rm = TRUE),
    .groups = "drop"
  )
# 490

xrf<-read.csv('./XIST_05/XISTrF_raw.csv',header=T)

# keep the rows/transcripts in xrf and exists in xra_collapsed
xrf<-semi_join(xrf,xra_collapsed,by='seqName')

# check length 
length(unique(xrf$seqName))

# filter for rows that has start <=3000
xrf<-xrf[which(xrf$Start<=3000),]

xrf<-xrf[,c(1,4)]

xrf_collapsed <- xrf %>%
  group_by(seqName) %>%
  summarise(
    Start   = min(Start,   na.rm = TRUE),
    .groups = "drop"
  )
# 140

# get the max for xra and xrf start as the mark for B12
xraf <- xra_collapsed %>%
  inner_join(xrf_collapsed, by = "seqName", suffix = c(".1", ".2")) %>%
  transmute(
    seqName,
    Start = pmax(Start.1, Start.2, na.rm = TRUE)
  )
# 140

xrb12<-read.csv('./XIST_05/XISTrB12_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
xrb12<-semi_join(xrb12,xraf,by='seqName')
# check length
length(unique(xrb12$seqName))

# Check if xrb12 has at least one start > xraf$start (for that same name)
# If yes, then among xrb1 rows with start > xraf$start, pick the smallest such start
# Output name + that chosen start

xrb12_fl <- xrb12 %>%
  inner_join(xraf, by = "seqName", suffix = c(".B12", ".AF")) %>%
  filter(Start.B12 > Start.AF) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.B12, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.B12)

# 126


xrd<-read.csv('./XIST_05/XISTrD_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
xrd<-semi_join(xrd,xrb12_fl,by='seqName')
# check length 
length(unique(xrd$seqName))

xrd_fl <- xrd %>%
  inner_join(xrb12_fl, by = "seqName", suffix = c(".D", ".B12")) %>%
  filter(Start.D > Start.B12) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.D, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.D)
# 99


xre<-read.csv('./XIST_05/XISTrE_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
xre<-semi_join(xre,xrd_fl,by='seqName')
# check length 
length(unique(xre$seqName))

xre_fl <- xre %>%
  inner_join(xrd_fl, by = "seqName", suffix = c(".E", ".D")) %>%
  filter(Start.E > Start.D) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.E, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.E)
# 88

comb_fl<-semi_join(comb,xre_fl,by='seqName')

write.csv(comb_fl,'XIST_pattern_matched_p05.csv',row.names = F)


chromlist<-c(paste0('chr',c(1:22)),'chrX','chrY')

templist<-strsplit(comb_fl$seqName,'_',fixed=T)
(tempchrom<-sapply(templist,'[[',2))
# all has common chrom


###################### process NEAT1 files
# contained a hit to NEAT1 domain #1 within its first 3kb, 
# followed by hits to domains #9 or 10, 
# followed by at least 1kb of hits to any domain from #11-16, 
# and lastly a hit do domain #22 At p < 0.05


#### first get all transcripts that contain hits to all features
filelist05<-dir(path='./NEAT1_05/',pattern='*._wide.csv')

for (n in 1:length(filelist05)) {
  
  temp<-read.csv(paste0('./NEAT1_05/',filelist05[n]),header=T)
  
  temp<-temp[,c("seqName","total_unique_coverage")]
  
  colnames(temp)[2]<-paste0(gsub('_wide.csv','',filelist05[n]),'_tuc0.05')
  
  if (exists('comb')) {
    # only keep transcripts that has hits in both dfs
    comb<-inner_join(comb, temp, by = "seqName")
  } else {
    comb<-temp
  }
}

chromlist<-c(paste0('chr',c(1:22)),'chrX','chrY')

templist<-strsplit(comb$seqName,'_',fixed=T)
tempchrom<-sapply(templist,'[[',2)
comb<-comb[(tempchrom %in% chromlist),]

write.csv(comb,'transcripts_withallNEAT1hits_p05.csv',row.names = F)
# 5182
sum(grepl('unspliced',comb$seqName)) # 4516

### go through the raw file one by one and check the pattern

n1<-read.csv('./NEAT1_05/NEAT1chk1_raw.csv',header=T)

# keep the rows/transcripts in xra that exists in comb
n1<-semi_join(n1,comb,by='seqName')

# check length 
length(unique(n1$seqName))

# filter for rows that has start <=3000
n1<-n1[which(n1$Start<=3000),]

n1<-n1[,c(1,4)]

n1_collapsed <- n1 %>%
  group_by(seqName) %>%
  summarise(
    Start   = min(Start,   na.rm = TRUE),
    .groups = "drop"
  )
# 4757

n910<-read.csv('./NEAT1_05/NEAT1chk910_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
n910<-semi_join(n910,n1_collapsed,by='seqName')
# check length 
length(unique(n910$seqName))

# Check if xrb1 has at least one start > xraf$start (for that same name)
# If yes, then among xrb1 rows with start > xraf$start, pick the smallest such start
# Output name + that chosen start

n910_fl <- n910 %>%
  inner_join(n1_collapsed, by = "seqName", suffix = c(".n910", ".n1")) %>%
  filter(Start.n910 > Start.n1) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.n910, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.n910)

# 4728

n1116<-read.csv('./NEAT1_05/NEAT1chk1116_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
n1116<-semi_join(n1116,n910_fl,by='seqName')
# check length 
length(unique(n1116$seqName))

# Check if xrb1 has at least one start > xraf$start (for that same name)
# If yes, then among xrb1 rows with start > xraf$start, pick the smallest such start
# Output name + that chosen start

n1116_fl <- n1116 %>%
  inner_join(n910_fl, by = "seqName", suffix = c(".n1116", ".n910")) %>%
  filter(Start.n1116 > Start.n910) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.n1116, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.n1116)
# 4602

# require the unique coverage > 1kb in n1116
tuc_n1116<-read.csv('./NEAT1_05/NEAT1chk1116_wide.csv',header=T)
tuc_n1116<-tuc_n1116[which(tuc_n1116$total_unique_coverage>=1000),]
# 6525
n1116_fl<-semi_join(n1116_fl,tuc_n1116,by='seqName')
# 4165


n22<-read.csv('./NEAT1_05/NEAT1chk22_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
n22<-semi_join(n22,n1116_fl,by='seqName')
# check length 
length(unique(n22$seqName))

n22_fl <- n22 %>%
  inner_join(n1116_fl, by = "seqName", suffix = c(".n22", ".n1116")) %>%
  filter(Start.n22 > Start.n1116) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.n22, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.n22)
# 2901


comb_fl<-semi_join(comb,n22_fl,by='seqName')

write.csv(comb_fl,'NEAT1_pattern_matched_p05.csv',row.names = F)

###########################################################
####### p0.01
setwd("XNM/pattern")

library(Biostrings)
library(dplyr)
library(tidyr)


###################### process XIST files
# transcripts that contained a hit to Repeats A and F within its first 3kb
# followed by hits to both Repeats B1/2 
# followed by hits to D, followed by a hit to Repeat E. At p < 0.05


#### first get all transcripts that contain hits to all features
filelist05<-dir(path='./XIST_01/',pattern='*._wide.csv')

for (n in 1:length(filelist05)) {
  
  temp<-read.csv(paste0('./XIST_01/',filelist05[n]),header=T)
  
  temp<-temp[,c("seqName","total_unique_coverage")]
  
  colnames(temp)[2]<-paste0(gsub('_wide.csv','',filelist05[n]),'_tuc0.01')
  
  if (exists('comb')) {
    # only keep transcripts that has hits in both dfs
    comb<-inner_join(comb, temp, by = "seqName")
  } else {
    comb<-temp
  }
}

write.csv(comb,'transcripts_withallXISThits_p01.csv',row.names = F)
# 10
sum(grepl('unspliced',comb$seqName)) # 9 only 1 spliced transcripts

### go through the raw file one by one and check the pattern

xra<-read.csv('./XIST_01/XISTrA_raw.csv',header=T)

# keep the rows/transcripts in xra that exists in comb
xra<-semi_join(xra,comb,by='seqName')

# check length 
length(unique(xra$seqName))

# filter for rows that has start <=3000
xra<-xra[which(xra$Start<=3000),]

xra<-xra[,c(1,4)]

xra_collapsed <- xra %>%
  group_by(seqName) %>%
  summarise(
    Start   = min(Start,   na.rm = TRUE),
    .groups = "drop"
  )
# 2

xrf<-read.csv('./XIST_01/XISTrF_raw.csv',header=T)

# keep the rows/transcripts in xrf and exists in xra_collapsed
xrf<-semi_join(xrf,xra_collapsed,by='seqName')

# check length 
length(unique(xrf$seqName))

# filter for rows that has start <=3000
xrf<-xrf[which(xrf$Start<=3000),]

xrf<-xrf[,c(1,4)]

xrf_collapsed <- xrf %>%
  group_by(seqName) %>%
  summarise(
    Start   = min(Start,   na.rm = TRUE),
    .groups = "drop"
  )
# 2

# get the max for xra and xrf start as the mark for B12
xraf <- xra_collapsed %>%
  inner_join(xrf_collapsed, by = "seqName", suffix = c(".1", ".2")) %>%
  transmute(
    seqName,
    Start = pmax(Start.1, Start.2, na.rm = TRUE)
  )
# 2

xrb12<-read.csv('./XIST_01/XISTrB12_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
xrb12<-semi_join(xrb12,xraf,by='seqName')
# check length
length(unique(xrb12$seqName))

# Check if xrb12 has at least one start > xraf$start (for that same name)
# If yes, then among xrb1 rows with start > xraf$start, pick the smallest such start
# Output name + that chosen start

xrb12_fl <- xrb12 %>%
  inner_join(xraf, by = "seqName", suffix = c(".B12", ".AF")) %>%
  filter(Start.B12 > Start.AF) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.B12, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.B12)

# 2


xrd<-read.csv('./XIST_01/XISTrD_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
xrd<-semi_join(xrd,xrb12_fl,by='seqName')
# check length 
length(unique(xrd$seqName))

xrd_fl <- xrd %>%
  inner_join(xrb12_fl, by = "seqName", suffix = c(".D", ".B12")) %>%
  filter(Start.D > Start.B12) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.D, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.D)
# 2


xre<-read.csv('./XIST_01/XISTrE_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
xre<-semi_join(xre,xrd_fl,by='seqName')
# check length 
length(unique(xre$seqName))

xre_fl <- xre %>%
  inner_join(xrd_fl, by = "seqName", suffix = c(".E", ".D")) %>%
  filter(Start.E > Start.D) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.E, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.E)
# 2

comb_fl<-semi_join(comb,xre_fl,by='seqName')

write.csv(comb_fl,'XIST_pattern_matched_p01.csv',row.names = F)


chromlist<-c(paste0('chr',c(1:22)),'chrX','chrY')

templist<-strsplit(comb_fl$seqName,'_',fixed=T)
(tempchrom<-sapply(templist,'[[',2))
# all has common chrom


###################### process NEAT1 files
# contained a hit to NEAT1 domain #1 within its first 3kb, 
# followed by hits to domains #9 or 10, 
# followed by at least 1kb of hits to any domain from #11-16, 
# and lastly a hit do domain #22 At p < 0.05


#### first get all transcripts that contain hits to all features
filelist05<-dir(path='./NEAT1_01/',pattern='*._wide.csv')

for (n in 1:length(filelist05)) {
  
  temp<-read.csv(paste0('./NEAT1_01/',filelist05[n]),header=T)
  
  temp<-temp[,c("seqName","total_unique_coverage")]
  
  colnames(temp)[2]<-paste0(gsub('_wide.csv','',filelist05[n]),'_tuc0.01')
  
  if (exists('comb')) {
    # only keep transcripts that has hits in both dfs
    comb<-inner_join(comb, temp, by = "seqName")
  } else {
    comb<-temp
  }
}

chromlist<-c(paste0('chr',c(1:22)),'chrX','chrY')

templist<-strsplit(comb$seqName,'_',fixed=T)
tempchrom<-sapply(templist,'[[',2)
comb<-comb[(tempchrom %in% chromlist),]

write.csv(comb,'transcripts_withallNEAT1hits_p01.csv',row.names = F)
# 1816
sum(grepl('unspliced',comb$seqName)) # 1770

### go through the raw file one by one and check the pattern

n1<-read.csv('./NEAT1_01/NEAT1chk1_raw.csv',header=T)

# keep the rows/transcripts in xra that exists in comb
n1<-semi_join(n1,comb,by='seqName')

# check length 
length(unique(n1$seqName))

# filter for rows that has start <=3000
n1<-n1[which(n1$Start<=3000),]

n1<-n1[,c(1,4)]

n1_collapsed <- n1 %>%
  group_by(seqName) %>%
  summarise(
    Start   = min(Start,   na.rm = TRUE),
    .groups = "drop"
  )
# 1531

n910<-read.csv('./NEAT1_01/NEAT1chk910_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
n910<-semi_join(n910,n1_collapsed,by='seqName')
# check length 
length(unique(n910$seqName))

# Check if xrb1 has at least one start > xraf$start (for that same name)
# If yes, then among xrb1 rows with start > xraf$start, pick the smallest such start
# Output name + that chosen start

n910_fl <- n910 %>%
  inner_join(n1_collapsed, by = "seqName", suffix = c(".n910", ".n1")) %>%
  filter(Start.n910 > Start.n1) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.n910, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.n910)

# 1527

n1116<-read.csv('./NEAT1_01/NEAT1chk1116_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
n1116<-semi_join(n1116,n910_fl,by='seqName')
# check length 
length(unique(n1116$seqName))

# Check if xrb1 has at least one start > xraf$start (for that same name)
# If yes, then among xrb1 rows with start > xraf$start, pick the smallest such start
# Output name + that chosen start

n1116_fl <- n1116 %>%
  inner_join(n910_fl, by = "seqName", suffix = c(".n1116", ".n910")) %>%
  filter(Start.n1116 > Start.n910) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.n1116, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.n1116)
# 1429

# require the unique coverage > 1kb in n1116
tuc_n1116<-read.csv('./NEAT1_01/NEAT1chk1116_wide.csv',header=T)
tuc_n1116<-tuc_n1116[which(tuc_n1116$total_unique_coverage>=1000),]
# 4708
n1116_fl<-semi_join(n1116_fl,tuc_n1116,by='seqName')
# 1251


n22<-read.csv('./NEAT1_01/NEAT1chk22_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
n22<-semi_join(n22,n1116_fl,by='seqName')
# check length 
length(unique(n22$seqName))

n22_fl <- n22 %>%
  inner_join(n1116_fl, by = "seqName", suffix = c(".n22", ".n1116")) %>%
  filter(Start.n22 > Start.n1116) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.n22, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.n22)
# 918


comb_fl<-semi_join(comb,n22_fl,by='seqName')

write.csv(comb_fl,'NEAT1_pattern_matched_p01.csv',row.names = F)



###########################################################
########### p0.001
setwd("XNM/pattern")

library(Biostrings)
library(dplyr)
library(tidyr)


###################### process XIST files
# transcripts that contained a hit to Repeats A and F within its first 3kb
# followed by hits to both Repeats B1/2 
# followed by hits to D, followed by a hit to Repeat E. At p < 0.05


#### first get all transcripts that contain hits to all features
filelist05<-dir(path='./XIST_001/',pattern='*._wide.csv')

for (n in 1:length(filelist05)) {
  
  temp<-read.csv(paste0('./XIST_001/',filelist05[n]),header=T)
  
  temp<-temp[,c("seqName","total_unique_coverage")]
  
  colnames(temp)[2]<-paste0(gsub('_wide.csv','',filelist05[n]),'_tuc0.001')
  
  if (exists('comb')) {
    # only keep transcripts that has hits in both dfs
    comb<-inner_join(comb, temp, by = "seqName")
  } else {
    comb<-temp
  }
}

write.csv(comb,'transcripts_withallXISThits_p001.csv',row.names = F)
# 2
sum(grepl('unspliced',comb$seqName)) # 1

### go through the raw file one by one and check the pattern

xra<-read.csv('./XIST_001/XISTrA_raw.csv',header=T)

# keep the rows/transcripts in xra that exists in comb
xra<-semi_join(xra,comb,by='seqName')

# check length 
length(unique(xra$seqName))

# filter for rows that has start <=3000
xra<-xra[which(xra$Start<=3000),]

xra<-xra[,c(1,4)]

xra_collapsed <- xra %>%
  group_by(seqName) %>%
  summarise(
    Start   = min(Start,   na.rm = TRUE),
    .groups = "drop"
  )
# 2

xrf<-read.csv('./XIST_001/XISTrF_raw.csv',header=T)

# keep the rows/transcripts in xrf and exists in xra_collapsed
xrf<-semi_join(xrf,xra_collapsed,by='seqName')

# check length 
length(unique(xrf$seqName))

# filter for rows that has start <=3000
xrf<-xrf[which(xrf$Start<=3000),]

xrf<-xrf[,c(1,4)]

xrf_collapsed <- xrf %>%
  group_by(seqName) %>%
  summarise(
    Start   = min(Start,   na.rm = TRUE),
    .groups = "drop"
  )
# 2

# get the max for xra and xrf start as the mark for B12
xraf <- xra_collapsed %>%
  inner_join(xrf_collapsed, by = "seqName", suffix = c(".1", ".2")) %>%
  transmute(
    seqName,
    Start = pmax(Start.1, Start.2, na.rm = TRUE)
  )
# 2

xrb12<-read.csv('./XIST_001/XISTrB12_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
xrb12<-semi_join(xrb12,xraf,by='seqName')
# check length
length(unique(xrb12$seqName))

# Check if xrb12 has at least one start > xraf$start (for that same name)
# If yes, then among xrb1 rows with start > xraf$start, pick the smallest such start
# Output name + that chosen start

xrb12_fl <- xrb12 %>%
  inner_join(xraf, by = "seqName", suffix = c(".B12", ".AF")) %>%
  filter(Start.B12 > Start.AF) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.B12, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.B12)

# 2


xrd<-read.csv('./XIST_001/XISTrD_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
xrd<-semi_join(xrd,xrb12_fl,by='seqName')
# check length 
length(unique(xrd$seqName))

xrd_fl <- xrd %>%
  inner_join(xrb12_fl, by = "seqName", suffix = c(".D", ".B12")) %>%
  filter(Start.D > Start.B12) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.D, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.D)
# 2


xre<-read.csv('./XIST_001/XISTrE_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
xre<-semi_join(xre,xrd_fl,by='seqName')
# check length 
length(unique(xre$seqName))

xre_fl <- xre %>%
  inner_join(xrd_fl, by = "seqName", suffix = c(".E", ".D")) %>%
  filter(Start.E > Start.D) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.E, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.E)
# 2

comb_fl<-semi_join(comb,xre_fl,by='seqName')

write.csv(comb_fl,'XIST_pattern_matched_p001.csv',row.names = F)


chromlist<-c(paste0('chr',c(1:22)),'chrX','chrY')

templist<-strsplit(comb_fl$seqName,'_',fixed=T)
tempchrom<-sapply(templist,'[[',2)
# all has common chrom


###################### process NEAT1 files
# contained a hit to NEAT1 domain #1 within its first 3kb, 
# followed by hits to domains #9 or 10, 
# followed by at least 1kb of hits to any domain from #11-16, 
# and lastly a hit do domain #22 At p < 0.05


#### first get all transcripts that contain hits to all features
filelist05<-dir(path='./NEAT1_001/',pattern='*._wide.csv')

for (n in 1:length(filelist05)) {
  
  temp<-read.csv(paste0('./NEAT1_001/',filelist05[n]),header=T)
  
  temp<-temp[,c("seqName","total_unique_coverage")]
  
  colnames(temp)[2]<-paste0(gsub('_wide.csv','',filelist05[n]),'_tuc0.001')
  
  if (exists('comb')) {
    # only keep transcripts that has hits in both dfs
    comb<-inner_join(comb, temp, by = "seqName")
  } else {
    comb<-temp
  }
}


chromlist<-c(paste0('chr',c(1:22)),'chrX','chrY')

templist<-strsplit(comb$seqName,'_',fixed=T)
tempchrom<-sapply(templist,'[[',2)
comb<-comb[(tempchrom %in% chromlist),]

write.csv(comb,'transcripts_withallNEAT1hits_p001.csv',row.names = F)
# 19
sum(grepl('unspliced',comb$seqName)) # 19

### go through the raw file one by one and check the pattern

n1<-read.csv('./NEAT1_001/NEAT1chk1_raw.csv',header=T)

# keep the rows/transcripts in xra that exists in comb
n1<-semi_join(n1,comb,by='seqName')

# check length 
length(unique(n1$seqName))

# filter for rows that has start <=3000
n1<-n1[which(n1$Start<=3000),]

n1<-n1[,c(1,4)]

n1_collapsed <- n1 %>%
  group_by(seqName) %>%
  summarise(
    Start   = min(Start,   na.rm = TRUE),
    .groups = "drop"
  )
# 6

n910<-read.csv('./NEAT1_001/NEAT1chk910_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
n910<-semi_join(n910,n1_collapsed,by='seqName')
# check length 
length(unique(n910$seqName))

# Check if xrb1 has at least one start > xraf$start (for that same name)
# If yes, then among xrb1 rows with start > xraf$start, pick the smallest such start
# Output name + that chosen start

n910_fl <- n910 %>%
  inner_join(n1_collapsed, by = "seqName", suffix = c(".n910", ".n1")) %>%
  filter(Start.n910 > Start.n1) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.n910, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.n910)

# 6

n1116<-read.csv('./NEAT1_001/NEAT1chk1116_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
n1116<-semi_join(n1116,n910_fl,by='seqName')
# check length 
length(unique(n1116$seqName))

# Check if xrb1 has at least one start > xraf$start (for that same name)
# If yes, then among xrb1 rows with start > xraf$start, pick the smallest such start
# Output name + that chosen start

n1116_fl <- n1116 %>%
  inner_join(n910_fl, by = "seqName", suffix = c(".n1116", ".n910")) %>%
  filter(Start.n1116 > Start.n910) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.n1116, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.n1116)
# 5

# require the unique coverage > 1kb in n1116
tuc_n1116<-read.csv('./NEAT1_001/NEAT1chk1116_wide.csv',header=T)
tuc_n1116<-tuc_n1116[which(tuc_n1116$total_unique_coverage>=1000),]
# 1619
n1116_fl<-semi_join(n1116_fl,tuc_n1116,by='seqName')
# 4


n22<-read.csv('./NEAT1_001/NEAT1chk22_raw.csv',header=T)

# keep the rows/transcripts exist after filtering
n22<-semi_join(n22,n1116_fl,by='seqName')
# check length 
length(unique(n22$seqName))

n22_fl <- n22 %>%
  inner_join(n1116_fl, by = "seqName", suffix = c(".n22", ".n1116")) %>%
  filter(Start.n22 > Start.n1116) %>%                          # keep only df2 starts above df1 start
  group_by(seqName) %>%
  slice_min(order_by = Start.n22, n = 1, with_ties = FALSE) %>%  # pick smallest such start
  ungroup() %>%
  transmute(seqName, Start = Start.n22)
# 3


comb_fl<-semi_join(comb,n22_fl,by='seqName')

write.csv(comb_fl,'NEAT1_pattern_matched_p001.csv',row.names = F)

####################################
# add gene biotype to the matched transcripts

setwd("XNM/pattern")

library(Biostrings) # needed to be loaded first, before rtracklayer
library(rtracklayer)
library(data.table) # better data handling
library(dplyr)

gtf_file<-'comprehensive_annotation_v47_unspliced_new_kcnq1ot1_1_22_25.gtf'

# extract info from gtf
gtf<-import(gtf_file)

# get the transcript id, gene id, strand, start and end coordinates for each exon
genes <- gtf[ mcols(gtf)$type == "gene" ]

genes <- as.data.table(genes)

add_gene_type<-function(filepath,genes) {
  
  df<-read.csv(filepath,header=T)
  
  temp<-strsplit(df$seqName,'_',fixed=T)
  
  tgeneid<-sapply(temp,'[[',6)
  
  df$gene_type<-''
  
  for (n in 1:length(tgeneid)) {
    
    tid<-tgeneid[n]
    
    df$gene_type[n]<-genes$gene_type[which(genes$gene_id==tid)]
  }
  
  write.csv(df,filepath,row.names = F)
}


add_gene_type('XIST_pattern_matched_p05.csv',genes)
add_gene_type('XIST_pattern_matched_p01.csv',genes)
add_gene_type('XIST_pattern_matched_p001.csv',genes)

add_gene_type('NEAT1_pattern_matched_p05.csv',genes)
add_gene_type('NEAT1_pattern_matched_p01.csv',genes)
add_gene_type('NEAT1_pattern_matched_p001.csv',genes)


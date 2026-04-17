# filter and process the gridsearch results
# copy over all csv output from the gridearch function
setwd("gridsearch")

filelist<-dir(pattern='*.csv',full.names = F)

transit_p<-data.frame(file=character(),qT=double(),nT=double(),
                      knum=integer(),total_n=integer(),
                      r_median=double(),r_std=double(),
                      len_median=integer(),len_std=double(),
                      top50_r_median=double(),tip50_r_std=double(),
                      top50_len_median=integer(),top50_len_std=double())

# initiate vector to store flagged file info
flagfile<-vector(mode='character',length=0)

for (file in filelist) {
  
  df<-read.csv(file,header=T)
  
  med_n<-median(df$total_n)
  
  maxmed<-max(df$top50_r_median)
  
  temp<-df[which(df$top50_r_median==maxmed),]
  
  if (nrow(temp)>1) {
    maxn<-max(temp$total_n)
    temp<-temp[which(temp$total_n==maxn),]
  }
  
  if (temp$total_n>=min(med_n,10000)) {
    temp$file<-gsub('.csv','',file)
    temp<-temp[,c(13,1:12)]
    transit_p<-rbind(transit_p,temp)
  } else {
    # if the length is short, do not write into transit_p
    # write the name into flagfile
    
    flagfile<-c(flagfile,file)
    
  }
  
  
} 

write.csv(transit_p,'optimized_transit_p.csv',quote=F,row.names = F)
write.table(flagfile,'flagfile.txt',sep='\t',row.names = F,col.names = F,quote=F)

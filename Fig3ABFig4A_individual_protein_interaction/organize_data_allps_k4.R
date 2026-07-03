# add query coordinates, pval, adjust pval, top5 RBP sig percentage
setwd("XNM")

plist<-c('p05','p01','p001')

for (p in plist) {
  
  dfname<-paste0('XNM_hmseekrhits_paired_unspliced_pval_CKnormed_',p,'.csv')
  print(dfname)
  df<-read.csv(dfname,header=T)
  
  colnames(df)[1]<-'query'
  
  combheaders <- read.csv('query_coords.csv',header=T)
  
  df$query_coords<-combheaders$query_coords
  
  stats<-read.csv('XNM_RBP_rpm_counts_CKnormed_chunks_filtered.csv',header=T)
  
  df$top5RBP<-stats$top5RBP
  
  df$top5rpm<-stats$top5rpm
  
  df$top5fc<-stats$top5fc
  
  htotalcount<-read.csv(paste0('XNM_hmseekrhits_paired_unspliced_ncount_CKnormed_',p,'.csv'),header=T)
  
  df$hit_totalcount<-htotalcount$AARS
  
  ttotalcount<-read.csv(paste0('XNM_hmseekrhits_paired_unspliced_transcript_count_CKnormed_',p,'.csv'),header=T)
  
  df$transcript_totalcount<-ttotalcount$AARS
  
  medlen<-read.csv(paste0('XNM_hmseekrhits_paired_unspliced_medlen_CKnormed_',p,'.csv'),header=T)
  
  df$medlen<-medlen$AARS
  
  df$top5RBP_pval<-''
  df$top5RBP_adjpval<-''
  # top5 RBP sig percentage
  df$top5hitpct<-NA
  
  for (n in 1:nrow(df)) {
    
    temp<-df[n,]
    
    trbp<-unlist(strsplit(temp$top5RBP,';',fixed=T))
    
    df$top5RBP_pval[n]<-paste(temp[trbp],collapse = ';')
    
    ptemp<-as.numeric(temp[trbp])
    
    ptemp.adj<-p.adjust(ptemp,method='BH')
    
    df$top5RBP_adjpval[n]<-paste0(ptemp.adj,collapse = ';')
    
    df$top5hitpct[n]<-sum(ptemp.adj<0.05)/length(ptemp.adj)
    
  }
  
  df<-df[,c(1,141,142,149,148,150,143:147,2:140)]
  
  write.csv(df,dfname,row.names = F)
  
}



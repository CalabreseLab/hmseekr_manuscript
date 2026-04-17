setwd("sigRBP_network")

# get the list of unique top5RBP that is sig in p0.05
# for the functional modules of X and N
df<-read.csv('XNmodule_hmseekrhits_paired_unspliced_pval_CKnormed_p05.csv',header=T)

for (n in 1:nrow(df)) {
  t5<-df$top5RBP[n]
  t5<-unlist(strsplit(t5,';',fixed=T))
  
  t5p<-df$top5RBP_adjpval[n]
  t5p<-as.numeric(unlist(strsplit(t5p,';',fixed=T)))
  
  if (exists('comb')) {
    comb<-rbind(comb,t5)
    combpval<-rbind(combpval,t5p)
  } else {
    comb<-t5
    combpval<-t5p
  }
  
}

comb<-as.data.frame(comb)
rownames(comb)<-df$query

combpval<-as.data.frame(combpval)
rownames(combpval)<-df$query

siglist<-comb[combpval<0.05]

siglist<-unique(siglist)

write.csv(siglist,'top5sigRBP_XNmodule.csv',row.names = F)

####### write out X and N specific
xcomb<-comb[grepl('XIST',rownames(comb)),]
xcombpval<-combpval[grepl('XIST',rownames(combpval)),]

xsiglist<-xcomb[xcombpval<0.05]
xsiglist<-unique(xsiglist)
write.csv(xsiglist,'top5sigRBP_Xmodule.csv',row.names = F)


ncomb<-comb[grepl('NEAT1',rownames(comb)),]
ncombpval<-combpval[grepl('NEAT1',rownames(combpval)),]

nsiglist<-ncomb[ncombpval<0.05]
nsiglist<-unique(nsiglist)
write.csv(nsiglist,'top5sigRBP_Nmodule.csv',row.names = F)

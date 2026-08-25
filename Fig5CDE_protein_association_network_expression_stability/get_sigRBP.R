setwd("sigRBP_network")


#######################################################
# K6


# get the list of unique top5RBP that is sig in p0.05
# for the functional modules of X and N
df<-read.csv('XNmodule_hmseekrhits_paired_unspliced_pval_CKnormed_p05_k6.csv',header=T)

if (exists('comb')) {rm(comb)}

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

#siglist<-comb[combpval<0.05]
#siglist<-unique(siglist)

siglist<-unique(unlist(comb,use.names = F))

write.csv(siglist,'top5RBP_XNmodule_k6.csv',row.names = F)

####### write out X and N specific
xcomb<-comb[grepl('XIST',rownames(comb)),]
xcombpval<-combpval[grepl('XIST',rownames(combpval)),]

#xsiglist<-xcomb[xcombpval<0.05]
#xsiglist<-unique(xsiglist)
xsiglist<-unique(unlist(xcomb,use.names = F))
write.csv(xsiglist,'top5RBP_Xmodule_k6.csv',row.names = F)


ncomb<-comb[grepl('NEAT1',rownames(comb)),]
ncombpval<-combpval[grepl('NEAT1',rownames(combpval)),]

#nsiglist<-ncomb[ncombpval<0.05]
#nsiglist<-unique(nsiglist)
nsiglist<-unique(unlist(ncomb,use.names = F))
write.csv(nsiglist,'top5RBP_Nmodule_k6.csv',row.names = F)

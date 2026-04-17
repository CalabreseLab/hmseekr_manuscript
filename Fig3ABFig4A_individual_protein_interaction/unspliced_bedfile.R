# generate bedfile for chunks and save score as negative log of the pval

# setwd set working directory properly
setwd("XNM")
library(Biostrings) # needed to be loaded first, before rtracklayer
library(rtracklayer)
library(data.table) # better data handling

# set so that no scientific format is used 
options(scipen=999)


##################################################################
# convert to a bedfile for pool seqs
##################################################################

# legnthfilter file
lenfilter<-read.csv('lengthfilter_XNM.csv',header=T)

files <- list.files(
  path = "./hits/",
  pattern = "^tsc_condE_newpool_unspliced_.*_4_viterbi_seekr\\.txt$",
  full.names = F
)

for (n in 1:length(files)) {
  
  print(paste0('files_',n))
  
  hits<-fread(paste0('./hits/',files[n]),sep='\t')
  
  chunkname<-files[n]
  chunkname<-gsub('tsc_condE_newpool_unspliced_','',chunkname)
  chunkname<-gsub('_4_viterbi_seekr.txt','',chunkname)
  
  lfidx<-which(lenfilter$query==chunkname)
  
  if(length(lfidx)>0) {
    lenmin<-lenfilter$lenmin[lfidx]
    lenmax<-lenfilter$lenmax[lfidx]
  } else {
    lenmin<-200
    lenmax<-6000
  }
  
  
  hits<-hits[which(hits$Length>lenmin & hits$Length<lenmax),]
  
  # only keep the p val <0.05 ones so that to lessen later step calculations
  hits<-hits[which(hits$seekr_pval<0.05),]
  
  hitsname<-gsub('>','',hits$seqName)
  
  hitsbed<-data.table(chrom=character(),
                      chromStart=integer(),
                      chromEnd=integer(),
                      name=character(),
                      score=numeric(),
                      strand=character())
  
  if (nrow(hits)<1) {
    
    print(paste0(chunkname,' has no hits after filter, write empty files'))
    
    outputname<-gsub('.txt','_hits.bed',files[n],fixed=T)
    fwrite(hitsbed, paste0('./bedfiles/',outputname), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
    
  } else {
    
    for (m in 1:nrow(hits)) {
      
      if(m %% 2000 ==0){print(m)}
      
      tempid<-hitsname[m]
      
      if (grepl('unspliced',tempid)) {
        
        # attach m to the gene name
        # there would be multiple hits for the same gene
        # attached m make sure each one is uniquely named
        # so that later when we split a hit into multiple segments
        # they can still be traced back
        name<-paste0(tempid,'_',m)
        
        tempid<-strsplit(tempid,'_')
        
        chrom<-tempid[[1]][2]
        strand<-tempid[[1]][5]
        score<-hits$seekr_pval[m]
        
        hstart<-as.numeric(hits$Start[m])
        hend<-as.numeric(hits$End[m])
        # hits using 1 start coords
        # 1, 20 is length 20 start at the first position
        if (strand=='+') {
          chromStart<-as.numeric(tempid[[1]][3])+hstart-1
          chromEnd<-as.numeric(tempid[[1]][3])+hend-1
        } else if (strand=='-') {
          chromEnd<-as.numeric(tempid[[1]][4])-hstart+1
          chromStart<-as.numeric(tempid[[1]][4])-hend+1
        } else {
          print('strand is wrong')
          print(n)
        }
        
        tempbed <- data.table(chrom, chromStart, chromEnd, name, score, strand)
        
        
        
      } else {
        # if gene is spliced we need function to calculate coords
        print('found spliced transcripts in the hits')
        print(n)
        print(m)
        break
        
      }
      
      
      hitsbed<-rbindlist(list(hitsbed,tempbed),fill=FALSE)
      
    }
    
    
    outputname<-gsub('.txt','_hits.bed',files[n],fixed=T)
    fwrite(hitsbed, paste0('./bedfiles/',outputname), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    
  }
    
}
  
  
  
  
  




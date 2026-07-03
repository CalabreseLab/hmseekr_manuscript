# generate bedfile for all transcripts in the unspliced pool

library(Biostrings) # needed to be loaded first, before rtracklayer
library(rtracklayer)
library(data.table) # better data handling




##################################################################
# convert to a bedfile for all transcripts in unspliced pool
##################################################################


# a function to create BED file based on info from tempbed
# create a BED file for each genes with 25 bp window for that gene
# for each row in tempbed generate a bed file with 25 bp window
# slide from start to end
# name col is: window_number of window
# score col is: 0
# strand col is: strand
# start is the start of the window
# end is the end of the window or the end of the gene, which ever is smaller
# if strand is +, slide from chromStart
# if strand is -, slide from chromEnd
# window_number is the order number of the window
create_bedfile <- function(tempbed, window_size) {
  # Ensure tempbed is a data.table
  tempbed <- as.data.table(tempbed)
  
  # Initialize a data.table to store the new rows
  result <- data.table()
  
  # Initialize window number counter
  window_number <- 1
  
  # Iterate over each row in tempbed
  for (i in 1:nrow(tempbed)) {
    # Extract the current row
    row <- tempbed[i]
    
    # Calculate the start and end positions for each window based on the strand
    if (row$strand == "+") {
      start_pos <- row$chromStart
      end_pos <- row$chromEnd
      while (start_pos < end_pos) {
        current_end <- min(start_pos + window_size - 1, end_pos)
        new_row <- data.table(
          chrom = row$chrom,
          chromStart = start_pos,
          chromEnd = current_end,
          name = paste0(row$name, '_', window_number),
          score = row$score,
          strand = row$strand
        )
        result <- rbind(result, new_row)
        start_pos <- start_pos + window_size
        window_number <- window_number + 1
      }
    } else if (row$strand == "-") {
      start_pos <- row$chromEnd
      end_pos <- row$chromStart
      while (start_pos > end_pos) {
        current_start <- max(start_pos - window_size + 1, end_pos)
        new_row <- data.table(
          chrom = row$chrom,
          chromStart = current_start,
          chromEnd = start_pos,
          name = paste0(row$name, '_', window_number),
          score = row$score,
          strand = row$strand
        )
        result <- rbind(result, new_row)
        start_pos <- start_pos - window_size
        window_number <- window_number + 1
      }
    }
  }
  
  return(result)
}


####################################
fa <- readDNAStringSet('v47_filtered_ot1fixed_unspliced.fa')

# Extract headers
hdr <- names(fa)


# set window size 
window_size<-25

# prevent scientific notion of numbers
options(scipen = 999)
# generate bedfiles for all genes in chunks
for (n in 1:length(hdr)) {
  
  if (n %% 1000 == 0) {print(n)}
  temp <- hdr[n]
    
  tempid<-strsplit(temp,'_')
  
  chrom<-tempid[[1]][2]
  name<-temp
  strand<-tempid[[1]][5]
  score<-0
  chromStart<-as.integer(tempid[[1]][3])
  chromEnd<-as.integer(tempid[[1]][4])
  
  
  tempbed <- data.table(chrom, chromStart, chromEnd, name, score, strand)
  finaltemp<-create_bedfile(tempbed, window_size)
  # this is 1 based both inclusive coords
  # bedfile is 0 based half inclusive coords
  finaltemp$chromStart<-finaltemp$chromStart-1
  fwrite(finaltemp, paste0('./bedfiles/',temp,'.bed'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}








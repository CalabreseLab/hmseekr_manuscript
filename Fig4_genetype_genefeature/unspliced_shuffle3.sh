#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-0
#SBATCH --mem=50g


module load bedtools/2.31.1


# Loop through each file name
for file in XNM/bedfiles/tsc_condE_newpool_unspliced_*.bed; do
    # Print the current filename
    echo "shuffle: $file"

    base=$(basename "$file" _hits_clean.bed)
    
    bedtools shuffle -i "$file" -g ../hg38_filtered.fa.fai -excl "$file" -incl v47_filtered_ot1fixed_unspliced_clean.bed > "./shufflebed1/${base}_shuffle1_clean.bed"

    bedtools shuffle -i "$file" -g ../hg38_filtered.fa.fai -excl "$file" -incl v47_filtered_ot1fixed_unspliced_clean.bed > "./shufflebed2/${base}_shuffle2_clean.bed"

    bedtools shuffle -i "$file" -g ../hg38_filtered.fa.fai -excl "$file" -incl v47_filtered_ot1fixed_unspliced_clean.bed > "./shufflebed3/${base}_shuffle3_clean.bed"
  
done




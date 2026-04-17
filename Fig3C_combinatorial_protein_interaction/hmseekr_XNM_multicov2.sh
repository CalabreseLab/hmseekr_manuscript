#!/bin/bash

module load bedtools/2.31.1

ls *_sorted.bam > multicov_files_list_XNM.txt

bam_files=$(cat multicov_files_list_XNM.txt | tr '\n' ' ')


for file in ../seqstosummary/bedfiles/tsc_condE_newpool_unspliced_*_4_viterbi_seekr_ck2.bed
do
    base=$(basename "$file" .bed)
    sbatch -p general --time=1-0 -n 8 -N 1 --mem=64g -o "../seqstosummary/multicov/${base}_multicov.out" --wrap="bedtools multicov -s -D -bams $bam_files -bed $file"

done

#!/bin/bash

module load bedtools

# run under where the bam files are saved

# all eclip bams from the analysis
ls *_sorted.bam > multicov_files_list_XNM.txt

bam_files=$(cat multicov_files_list_XNM.txt | tr '\n' ' ')


sbatch -p general --time=1-0 -n 8 -N 1 --mem=128g -o "XNM/XNM_chunk_bedfile_multicov.out" --wrap="bedtools multicov -s -D -bams $bam_files -bed XNM/XNM_chunk_bedfile.bed"

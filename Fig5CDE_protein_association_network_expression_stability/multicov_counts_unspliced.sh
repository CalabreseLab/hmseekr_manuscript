#!/bin/bash

module load bedtools/2.31.1

ls *.bam > multicov_files_list_XNpattern.txt

bam_files=$(cat multicov_files_list_XNpattern.txt | tr '\n' ' ')

for file in ../bedfiles/*
do
    base=$(basename "$file" .bed)
    sbatch -p general --time=3:00:00 -n 8 -N 1 --mem=64g -o "../multicov_out/${base}_multicov.out" --wrap="bedtools multicov -s -D -q 30 -bams $bam_files -bed $file"
done

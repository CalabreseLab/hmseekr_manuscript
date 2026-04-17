module load bedtools


sort -k1,1 -k2,2n lncRNA_raw.bed > lncRNA_sorted.bed
sort -k1,1 -k2,2n other_raw.bed > other_sorted.bed
sort -k1,1 -k2,2n protein_coding_raw.bed > protein_coding_sorted.bed


bedtools merge -i lncRNA_sorted.bed -s -c 4,5,6 -o distinct,distinct,distinct > lncRNA.bed
bedtools merge -i other_sorted.bed -s -c 4,5,6 -o distinct,distinct,distinct > other.bed
bedtools merge -i protein_coding_sorted.bed -s -c 4,5,6 -o distinct,distinct,distinct > protein_coding.bed




#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-0
#SBATCH --mem=150g


module load bedtools

# List of file names
xheaders=('XISTint1' 'XISTrA' 'XISTrF' 'XISTint2' 'XISTrB1' 'XISTint3' 'XISTrB2' 'XISTint4' 'XISTint5' 'XISTrD' 'XISTint6' 'XISTint7' 'XISTint8' 'XISTint9' 'XISTrE' 'XISTint10' 'XISTint11' 'XISTint12' 'XISTint13' 'XISTint14')

nheaders=$(seq -f "NEAT1chk%g" 1 22)
mheaders=$(seq -f "MALAT1chk%g" 1 9)



# Combine all headers into a single list
xnmfiles=( "${xheaders[@]}" $(echo $nheaders) $(echo $mheaders))



# Loop through each file name
for xnmfilename in "${xnmfiles[@]}"
do

    echo "Processing file: $xnmfilename"
    
    bedtools shuffle -i "./bedfiles/tsc_condE_newpool_unspliced_${xnmfilename}_4_viterbi_seekr_hits.bed" -g ../../hg38_filtered.fa.fai -excl "./bedfiles/tsc_condE_newpool_unspliced_${xnmfilename}_4_viterbi_seekr_hits.bed" -incl ./v47_filtered_ot1fixed_unspliced.bed > "./bedfiles/tsc_condE_newpool_unspliced_${xnmfilename}_4_viterbi_seekr_ck2.bed"

  
done





#!/bin/bash

module load python/3.12.1


# List of file names
# Define the xheaders array

xheaders=('XISTint1' 'XISTrA' 'XISTrF' 'XISTint2' 'XISTrB1' 'XISTint3' 'XISTrB2' 'XISTint4' 'XISTint5' 'XISTrD' 'XISTint6' 'XISTint7' 'XISTint8' 'XISTint9' 'XISTrE' 'XISTint10' 'XISTint11' 'XISTint12' 'XISTint13' 'XISTint14')

nheaders=$(seq -f "NEAT1chk%g" 1 22)
mheaders=$(seq -f "MALAT1chk%g" 1 9)



# Combine all headers into a single list
xnmfiles=( "${xheaders[@]}" $(echo $nheaders) $(echo $mheaders))



# Loop through each file name
for xnmfilename in "${xnmfiles[@]}"
do

  sbatch -p general --time=24:00:00 -n 8 -N 1 --mem=96g --wrap="hmseekr_findhits_condE -pool './v47_filtered_ot1fixed_unspliced.fa' -m './markovModels/${xnmfilename}_v47TSC/4/hmm.dict' -k 4 -name tsc_condE_newpool_unspliced_${xnmfilename} -dir './hits/' -a 'ATCG' -fa -pb"
  
done

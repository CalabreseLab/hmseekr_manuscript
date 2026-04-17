#!/bin/bash

module load python/3.12.1

# List of file names
xheaders=('XISTint1' 'XISTrA' 'XISTrF' 'XISTint2' 'XISTrB1' 'XISTint3' 'XISTrB2' 'XISTint4' 'XISTint5' 'XISTrD' 'XISTint6' 'XISTint7' 'XISTint8' 'XISTint9' 'XISTrE' 'XISTint10' 'XISTint11' 'XISTint12' 'XISTint13' 'XISTint14')

nheaders=$(seq -f "NEAT1chk%g" 1 22)
mheaders=$(seq -f "MALAT1chk%g" 1 9)

# Combine all headers into a single list
xnmfiles=( "${xheaders[@]}" $(echo $nheaders) $(echo $mheaders))



# Loop through each file name
for xnmfilename in "${xnmfiles[@]}"
do

  sbatch -p general --time=3:00:00 -n 8 -N 1 --mem=96g --wrap="hmseekr_hitseekr -hd './hits/tsc_condE_newpool_unspliced_${xnmfilename}_4_viterbi.txt' -qf './fastaFiles/${xnmfilename}.fa' -bkgf './v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -li 0 -la 1000000 -pf 1.1 -rf -1.1 -dir './hits/' -pb"

done


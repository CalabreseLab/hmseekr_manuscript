#!/bin/bash

module load python/3.12.1


# Directory to store the output
base_dir="gridsearch"

# Loop through all files containing 'int' or 'chk' in their names
for fasta_file in *{int,chk}*.fa; do
    # XIST int1239 is short and nees special li and la
    if [[ ! $fasta_file =~ int(1|2|3|9)\.fa$ ]]; then
        # Extract the filename without the extension
        filename=$(basename "$fasta_file" .fa)

        # Construct the name and dir arguments
        job_name="${filename}"
        job_dir="${base_dir}/gridsearch_${filename}/"

        # Submit the job
        sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf '${fasta_file}' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 200 -la 6000 -name '${job_name}' -dir '${job_dir}' -a 'ATCG'"
    fi
done



# manually submit for short queries as we need to change the li and la argument

sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XISTrA.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 100 -la 1800 -name 'XISTrA' -dir 'gridsearch/gridsearch_XISTrA/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XISTrB1.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 800 -name 'XISTrB1' -dir 'gridsearch/gridsearch_XISTrB1/' -a 'ATCG'"

sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XISTrB2.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 1000 -name 'XISTrB2' -dir 'gridsearch/gridsearch_XISTrB2/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XISTrD.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 200 -la 6000 -name 'XISTrD' -dir 'gridsearch/gridsearch_XISTrD/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XISTrE.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 150 -la 2800 -name 'XISTrE' -dir 'gridsearch/gridsearch_XISTrE/' -a 'ATCG'"

sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XISTrF.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 200 -la 6000 -name 'XISTrF' -dir 'gridsearch/gridsearch_XISTrF/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XISTint1.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 75 -la 1500 -name 'XISTint1' -dir 'gridsearch/gridsearch_XISTint1/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XISTint2.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 1000 -name 'XISTint2' -dir 'gridsearch/gridsearch_XISTint2/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XISTint3.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 185 -la 3000 -name 'XISTint3' -dir 'gridsearch/gridsearch_XISTint3/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XISTint9.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 145 -la 2500 -name 'XISTint9' -dir 'gridsearch/gridsearch_XISTint9/' -a 'ATCG'"

sbatch -p general --time=3-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'HOTTIP-MLL.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 200 -la 6000 -name 'HOTTIP-MLL' -dir 'gridsearch/gridsearch_HOTTIP-MLL/' -a 'ATCG'"

sbatch -p general --time=3-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'XIST-PfIMI.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 200 -la 6000 -name 'XIST-PfIMI' -dir 'gridsearch/gridsearch_XIST-PfIMI/' -a 'ATCG'"

sbatch -p general --time=3-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'MALAT1E.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 200 -la 6000 -name 'MALAT1E' -dir 'gridsearch/gridsearch_MALAT1E/' -a 'ATCG'"

sbatch -p general --time=3-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'MALAT1M.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 200 -la 6000 -name 'MALAT1M' -dir 'gridsearch/gridsearch_MALAT1M/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'ncRNA-a7.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 150 -la 2500 -name 'ncRNA-a7' -dir 'gridsearch/gridsearch_ncRNA-a7/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'FIRRE-RRD.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 1000 -name 'FIRRE-RRD' -dir 'gridsearch/gridsearch_FIRRE-RRD/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'MEG3-NRE.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 60 -la 1200 -name 'MEG3-NRE' -dir 'gridsearch/gridsearch_MEG3-NRE/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'NXF1-enChr.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 1000 -name 'NXF1-enChr' -dir 'gridsearch/gridsearch_NXF1-enChr/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'JPX-9.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 1000 -name 'JPX-9' -dir 'gridsearch/gridsearch_JPX-9/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'PVT1-22.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 1000 -name 'PVT1-22' -dir 'gridsearch/gridsearch_PVT1-22/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'NKILA-SRSFNRE.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 50 -la 1000 -name 'NKILA-SRSFNRE' -dir 'gridsearch/gridsearch_NKILA-SRSFNRE/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'e-Ccnd1.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 100 -la 2000 -name 'e-Ccnd1' -dir 'gridsearch/gridsearch_e-Ccnd1/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'e-YY1.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 90 -la 1800 -name 'e-YY1' -dir 'gridsearch/gridsearch_e-YY1/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'e-Med13l.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 1000 -name 'e-Med13l' -dir 'gridsearch/gridsearch_e-Med13l/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'e-Klf6.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 100 -la 2000 -name 'e-Klf6' -dir 'gridsearch/gridsearch_e-Klf6/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'e-Sp3.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 130 -la 2500 -name 'e-Sp3' -dir 'gridsearch/gridsearch_e-Sp3/' -a 'ATCG'"



sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'e-ID1.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 100 -la 2000 -name 'e-ID1' -dir 'gridsearch/gridsearch_e-ID1/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'PRRX2-eRNA.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 80 -la 1500 -name 'PRRX2-eRNA' -dir 'gridsearch/gridsearch_PRRX2-eRNA/' -a 'ATCG'"


sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'SEMA3C-eRNA.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 100 -la 2000 -name 'SEMA3C-eRNA' -dir 'gridsearch/gridsearch_SEMA3C-eRNA/' -a 'ATCG'"



sbatch -p general --time=2-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'Blustr-5ss.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 50 -la 1000 -name 'Blustr-5ss' -dir 'gridsearch/gridsearch_Blustr-5ss/' -a 'ATCG'"



sbatch -p general --time=3-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'HSat3_sense.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 6000 -name 'HSat3_sense' -dir 'gridsearch/gridsearch_HSat3_sense/' -a 'ATCG'"

sbatch -p general --time=3-0 -n 8 -N 1 --mem=128g --wrap="hmseekr_gridsearch -qf 'HSat3_antisense.fa' -nf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -pool 'v47_filtered_ot1fixed_combined.fa' -bkgf 'v47.lncRNA.can.500.nodup.OT1fixed.fa' -k 4 -ql 0.9,0.99,0.03 -nl 0.9,0.99,0.03 -step -fc findhits_condE -li 25 -la 6000 -name 'HSat3_antisense' -dir 'gridsearch/gridsearch_HSat3_antisense/' -a 'ATCG'"


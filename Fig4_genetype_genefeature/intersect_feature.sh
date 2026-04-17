#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-0
#SBATCH --mem=50g

module load bedtools/2.31.1

# Set the folder where your feature bed files are located
FEATURE_FOLDER="XNM/v47_gene_features"  

# for p 0.05 files
# Loop over each feature bed file
for feature in "$FEATURE_FOLDER"/*_filtered.bed; do

	echo "Processing feature file: $feature"
    # Get the base name (without path and .bed extension)
    base_feature=$(basename "$feature" _filtered.bed)
    
    # Create a directory in the current folder named after the feature
    mkdir -p ./realdata05/rd_${base_feature}

    # Define a CSV file with header for this feature
    csv_file="./realdata05/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    for hits in XNM/bedfiles/*_unspliced_*_seekr_hits_clean_filtered0.05.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./realdata05/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        # -c For each entry in A, report the number of hits in B that overlaps at least 1bp. Reports 0 for A entries that have no overlap with B.
        # -s Force “strandedness”. That is, only report hits in B that overlap A on the same strand. 
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"

        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"


    done
done


# for p0.01 files
# Loop over each feature bed file
for feature in "$FEATURE_FOLDER"/*_filtered.bed; do

    echo "Processing feature file: $feature"
    # Get the base name (without path and .bed extension)
    base_feature=$(basename "$feature" _filtered.bed)
    
    # Create a directory in the current folder named after the feature
    mkdir -p ./realdata01/rd_${base_feature}

    # Define a CSV file with header for this feature
    csv_file="./realdata01/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    for hits in XNM/bedfiles/*_unspliced_*_seekr_hits_clean_filtered0.01.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./realdata01/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        # -c For each entry in A, report the number of hits in B that overlaps at least 1bp. Reports 0 for A entries that have no overlap with B.
        # -s Force “strandedness”. That is, only report hits in B that overlap A on the same strand. 
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"

        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"


    done
done



# for p0.001 files
# Loop over each feature bed file
for feature in "$FEATURE_FOLDER"/*_filtered.bed; do

    echo "Processing feature file: $feature"
    # Get the base name (without path and .bed extension)
    base_feature=$(basename "$feature" _filtered.bed)
    
    # Create a directory in the current folder named after the feature
    mkdir -p ./realdata001/rd_${base_feature}

    # Define a CSV file with header for this feature
    csv_file="./realdata001/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    for hits in XNM/bedfiles/*_unspliced_*_seekr_hits_clean_filtered0.001.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./realdata001/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        # -c For each entry in A, report the number of hits in B that overlaps at least 1bp. Reports 0 for A entries that have no overlap with B.
        # -s Force “strandedness”. That is, only report hits in B that overlap A on the same strand. 
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"

        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"


    done
done
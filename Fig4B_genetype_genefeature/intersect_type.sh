#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-0
#SBATCH --mem=50g

module load bedtools


# Set the folder where your feature bed files are located
FEATURE_FOLDER="v47_gene_types"  

# for p 0.05 files
# Loop over each feature bed file
for feature in "$FEATURE_FOLDER"/*.bed; do

	echo "Processing feature file: $feature"
    # Get the base name (without path and .bed extension)
    base_feature=$(basename "$feature" _clean.bed)
    
    # Create a directory in the current folder named after the feature
    mkdir -p ./realtype05/rd_${base_feature}

    # Define a CSV file with header for this feature
    csv_file="./realtype05/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    # only do unspliced
    for hits in XNM/bedfiles/*_unspliced_*_seekr_hits_clean_filtered0.05.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./realtype05/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"


        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"


    done
done


# for p 0.01 files
# Loop over each feature bed file
for feature in "$FEATURE_FOLDER"/*.bed; do

    echo "Processing feature file: $feature"
    # Get the base name (without path and .bed extension)
    base_feature=$(basename "$feature" _clean.bed)
    
    # Create a directory in the current folder named after the feature
    mkdir -p ./realtype01/rd_${base_feature}

    # Define a CSV file with header for this feature
    csv_file="./realtype01/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    # only do unspliced
    for hits in XNM/bedfiles/*_unspliced_*_seekr_hits_clean_filtered0.01.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./realtype01/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"


        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"



    done
done



# for p 0.001 files
# Loop over each feature bed file
for feature in "$FEATURE_FOLDER"/*.bed; do

    echo "Processing feature file: $feature"
    # Get the base name (without path and .bed extension)
    base_feature=$(basename "$feature" _clean.bed)
    
    # Create a directory in the current folder named after the feature
    mkdir -p ./realtype001/rd_${base_feature}

    # Define a CSV file with header for this feature
    csv_file="./realtype001/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    # only do unspliced
    for hits in XNM/bedfiles/*_unspliced_*_seekr_hits_clean_filtered0.001.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./realtype001/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"


        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"



    done
done


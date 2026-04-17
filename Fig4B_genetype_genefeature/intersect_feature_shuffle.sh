#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2-0
#SBATCH --mem=96g

module load bedtools/2.31.1

# Set the folder where your feature bed files are located
FEATURE_FOLDER="XNM/v47_gene_features"  

##############################################
# for p 0.05
# Loop over each feature bed file
echo "Processing p0.05"
for feature in "$FEATURE_FOLDER"/*_filtered.bed; do

	echo "Processing feature file: $feature"
    # Get the base name (without path and .bed extension)
    base_feature=$(basename "$feature" _filtered.bed)
    
    # Create a directory in the current folder named after the feature
    mkdir -p ./shuffledata1_05/rd_${base_feature}
    mkdir -p ./shuffledata2_05/rd_${base_feature}
    mkdir -p ./shuffledata3_05/rd_${base_feature}

    # Define a CSV file with header for this feature
    csv_file="./shuffledata1_05/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"


    echo "shuffle1"
    
    # Loop over each hits bed file in the current directory
    for hits in ./shufflebed1/*_seekr_shuffle1_clean_filtered0.05.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./shuffledata1_05/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"

        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"

    done

    echo "shuffle2"
    # Define a CSV file with header for this feature
    csv_file="./shuffledata2_05/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    for hits in ./shufflebed2/*_seekr_shuffle2_clean_filtered0.05.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./shuffledata2_05/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"


        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"

    done

    echo "shuffle3"
    # Define a CSV file with header for this feature
    csv_file="./shuffledata3_05/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    for hits in ./shufflebed3/*_seekr_shuffle3_clean_filtered0.05.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./shuffledata3_05/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
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



##############################################
# for p 0.01
# Loop over each feature bed file
echo "Processing p0.01"

for feature in "$FEATURE_FOLDER"/*_filtered.bed; do

    echo "Processing feature file: $feature"
    # Get the base name (without path and .bed extension)
    base_feature=$(basename "$feature" _filtered.bed)
    
    # Create a directory in the current folder named after the feature
    mkdir -p ./shuffledata1_01/rd_${base_feature}
    mkdir -p ./shuffledata2_01/rd_${base_feature}
    mkdir -p ./shuffledata3_01/rd_${base_feature}

    # Define a CSV file with header for this feature
    csv_file="./shuffledata1_01/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"


    echo "shuffle1"
    
    # Loop over each hits bed file in the current directory
    for hits in ./shufflebed1/*_seekr_shuffle1_clean_filtered0.01.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./shuffledata1_01/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"

        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"

    done

    echo "shuffle2"
    # Define a CSV file with header for this feature
    csv_file="./shuffledata2_01/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    for hits in ./shufflebed2/*_seekr_shuffle2_clean_filtered0.01.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./shuffledata2_01/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"


        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"

    done

    echo "shuffle3"
    # Define a CSV file with header for this feature
    csv_file="./shuffledata3_01/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    for hits in ./shufflebed3/*_seekr_shuffle3_clean_filtered0.01.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./shuffledata3_01/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
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




##############################################
# for p 0.001
# Loop over each feature bed file
echo "Processing p0.001"

for feature in "$FEATURE_FOLDER"/*_filtered.bed; do

    echo "Processing feature file: $feature"
    # Get the base name (without path and .bed extension)
    base_feature=$(basename "$feature" _filtered.bed)
    
    # Create a directory in the current folder named after the feature
    mkdir -p ./shuffledata1_001/rd_${base_feature}
    mkdir -p ./shuffledata2_001/rd_${base_feature}
    mkdir -p ./shuffledata3_001/rd_${base_feature}

    # Define a CSV file with header for this feature
    csv_file="./shuffledata1_001/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"


    echo "shuffle1"
    
    # Loop over each hits bed file in the current directory
    for hits in ./shufflebed1/*_seekr_shuffle1_clean_filtered0.001.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./shuffledata1_001/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"

        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"

    done

    echo "shuffle2"
    # Define a CSV file with header for this feature
    csv_file="./shuffledata2_001/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    for hits in ./shufflebed2/*_seekr_shuffle2_clean_filtered0.001.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./shuffledata2_001/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
        # Run bedtools intersect with the -c -s option:
        bedtools intersect -a "$hits" -b "$feature" -c -s > "$output"


        # Count rows with an overlap (last column > 0)
        overlap=$(awk '($NF+0)>0 {c++} END {print c+0}' "$output")
        # Count rows without an overlap (last column == 0)
        nonoverlap=$(awk '($NF+0)==0 {c++} END {print c+0}' "$output")
        
        # Append the stats for this hits file to the CSV file
        echo "${hits_base},${overlap},${nonoverlap}" >> "$csv_file"

    done

    echo "shuffle3"
    # Define a CSV file with header for this feature
    csv_file="./shuffledata3_001/unsplice_${base_feature}_stats.csv"
    echo "file,overlap,nonoverlap" > "$csv_file"
    
    # Loop over each hits bed file in the current directory
    for hits in ./shufflebed3/*_seekr_shuffle3_clean_filtered0.001.bed; do
        hits_base=$(basename "$hits" .bed)
        # Define an output file name; adjust the naming as needed
        output="./shuffledata3_001/rd_${base_feature}/${hits_base}_${base_feature}_intersect.bed"
        
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
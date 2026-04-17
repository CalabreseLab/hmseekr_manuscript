# filter unspliced bedfiles so that score <0.05 0.01 or 0.001
# Loop over all matching files in XNM

for file in XNM/bedfiles/*_unspliced_*_seekr_hits_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.05' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.05.bed"
done


for file in XNM/bedfiles/*_unspliced_*_seekr_hits_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.01' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.01.bed"
done

for file in XNM/bedfiles/*_unspliced_*_seekr_hits_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.001' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.001.bed"
done


# when shuffling with bedtools shuffle the name and score does not change
# so the score in col5 are the same as hits
# filter this way would be the same as filtering on the score of hits

cd ../shufflebed1
for file in XNM/shufflebed1/*_seekr_shuffle1_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.05' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.05.bed"
done

for file in XNM/shufflebed1/*_seekr_shuffle1_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.01' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.01.bed"
done

for file in XNM/shufflebed1/*_seekr_shuffle1_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.001' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.001.bed"
done

cd ../shufflebed2
for file in XNM/shufflebed2/*_seekr_shuffle2_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.05' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.05.bed"
done

for file in XNM/shufflebed2/*_seekr_shuffle2_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.01' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.01.bed"
done

for file in XNM/shufflebed2/*_seekr_shuffle2_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.001' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.001.bed"
done


cd ../shufflebed3
for file in XNM/shufflebed3/*_seekr_shuffle3_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.05' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.05.bed"
done

for file in XNM/shufflebed3/*_seekr_shuffle3_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.01' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.01.bed"
done

for file in XNM/shufflebed3/*_seekr_shuffle3_clean.bed; do
    base=$(basename "$file")
    awk '$5 < 0.001' "$file" > "XNM/bedfiles/${base%.bed}_filtered0.001.bed"
done



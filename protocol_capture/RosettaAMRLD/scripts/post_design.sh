#!/bin/bash
# post_design.sh <ligand name>
# This script put together a ranked score file after a design trial. Run this script under the production folder

currLIG=$1

# Navigate to the ligand directory
cd $currLIG/

# Remove existing consolidated score file if it exists
[ -e design_scores.sc ] && rm design_scores.sc
echo "Design	total_score	interface_delta_X	LE2(REU)	SMILES" >> design_scores_header.temp

# Get a list of all score files (including those without a prefix)
score_files=$(ls *_score.sc score.sc 2>/dev/null)

# Loop over all score files
for file in $score_files
do
  # Extract relevant scores and SMILES for each entry in the score file
  while read -r line; do
    filename=$(echo "$line" | awk '{print $1}')
    interface=$(echo "$line" | awk '{print $3}')
    total=$(echo "$line" | awk '{print $2}')
    smi=$(awk -F"The best scoring ligand is " '/The best scoring ligand is /{print $2}' $filename.log)
    nHAtoms=$(obabel -:"${smi}" -osmi -d --append "atoms" | awk '{print $NF}')
    LE=$(echo "scale=3;$interface / sqrt($nHAtoms)" | bc)
    
    # Append the results to the consolidated score file
    echo "$filename	$total	$interface	$LE	$smi" >> design_scores_no_header.sc
  done < <(../../scripts/extract_scores.bash $file)
done

# Sort the consolidated file based on the fourth column (LE score)
sort -nk4 design_scores_no_header.sc >> design_scores_header.temp

# Rename the final file to design_scores.sc
mv design_scores_header.temp design_scores.sc

# Clean up temporary files
rm design_scores_no_header.sc

../../scripts/draw_top_design.py -i design_scores.sc --output ../$currLIG.png

#!/bin/bash

# Create the output directory if it doesn't exist
mkdir -p ./html

# Loop through each sample and render the .qmd file
for sample in iESC_0h iESC_6h iESC_12h iESC_24h iESC_48h iESC_72h iESC_96h iESC_120h; do
  # Create a temporary directory for each sample
  temp_dir=$(mktemp -d)

  # Render the file to the temporary directory
  quarto render QC.qmd -P sample="$sample" --output-dir "$temp_dir" -o "QC_${sample}.html"

  # Move the rendered file to the html directory
  mv "$temp_dir/QC_${sample}.html" ./html/

  # Remove the temporary directory
  rm -rf "$temp_dir"
done

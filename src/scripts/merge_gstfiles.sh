#!/bin/bash

# Directory containing the ROOT files
input_dir="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/FinalSPSPiAnalysis/Rarita/Carbon/4GeV/"
output_dir="/pnfs/genie/persistent/users/jtenavid/e4nu_files/GENIE_Files/2024Generation/FinalSPSPiAnalysis/Rarita/Carbon/4GeV/Grouped/"
num_chunks=10  # Number of merged output files

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Get all ROOT files in the directory
files=($(ls "$input_dir"*.root))

# Total number of files
total_files=${#files[@]}

# Number of files per chunk (roughly)
files_per_chunk=$(( (total_files + num_chunks - 1) / num_chunks ))

# Merge files in chunks
for ((i = 0; i < num_chunks; i++)); do
    # Calculate the range of files for this chunk
    start=$((i * files_per_chunk))
    end=$((start + files_per_chunk))
    if [ "$end" -gt "$total_files" ]; then
        end=$total_files
    fi

    # Get the list of files for this chunk
    chunk_files=("${files[@]:$start:$((end - start))}")

    # Skip empty chunks
    if [ ${#chunk_files[@]} -eq 0 ]; then
        continue
    fi

    # Prepare the output filename
    output_file="${output_dir}merged_chunk_$((i + 1)).gst.root"

    # Run hadd for this chunk
    echo "Merging files ${start} to $((end - 1)) into $output_file"
    hadd -f "$output_file" "${chunk_files[@]}"
done

echo "All files have been merged into $num_chunks chunks in $output_dir"

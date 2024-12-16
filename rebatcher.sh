#!/bin/bash

# Default folder name
BASE_FOLDER="auto_generated"

# Get suffix from user input
echo "Enter suffix (or leave empty for none): "
read SUFFIX

if [[ -n "$SUFFIX" ]]; then
    FOLDER="${BASE_FOLDER}_${SUFFIX}"
else
    FOLDER="$BASE_FOLDER"
fi

# Check if the folder exists
if [[ ! -d "$FOLDER" ]]; then
    echo "Error: Folder '$FOLDER' does not exist."
    exit 1
fi

# Loop through all .slurm files in the folder and submit them
for slurm_script in "$FOLDER"/*.slurm; do
    if [[ -f "$slurm_script" ]]; then  # Check if it's a regular file
        echo "Submitting $slurm_script"
        sbatch "$slurm_script"
    fi
done

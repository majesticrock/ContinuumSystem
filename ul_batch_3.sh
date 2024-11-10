#!/bin/bash
rm -rf auto_generated_ul_3/
mkdir -p auto_generated_ul_3
# Define the name of the input file
input_file="params/for_auto.txt"

# Read the lines of the input file into an array of strings
readarray -t NEW_VALUES < "${input_file}"
TOKEN="phonon_coupling"


for NEW_VALUE in "${NEW_VALUES[@]}"; do
  NEW_NAME=$(echo "$NEW_VALUE" | sed 's/ /_/g')
  # Loop through each line in the config file
  while read line; do
    if [[ $line == \#* ]]; then
      continue
    fi

    # Split the line into token and value
    TOKEN_NAME=$(echo "$line" | awk '{print $1}')
    TOKEN_VALUE=$(echo "$line" | cut -d' ' -f2-)

    if [[ "$TOKEN_NAME" == "$TOKEN" ]]; then
      # replace the value with the new one
      sed "s/$TOKEN_NAME $TOKEN_VALUE/$TOKEN_NAME $NEW_VALUE/" params/ul_cluster.config > auto_generated_ul_3/$NEW_NAME.config
      break
    fi
  done < params/ul_cluster.config
  cp slurm/ul_modes.slurm auto_generated_ul_3/$NEW_NAME.slurm
  sed -i "s|#SBATCH --job-name=modes|#SBATCH --job-name=ul3_$NEW_NAME|" auto_generated_ul_3/$NEW_NAME.slurm
  sed -i "s|#SBATCH --output=/home/althueser/phd/cpp/ContinuumSystem/modes_output.txt|#SBATCH --output=/home/althueser/phd/cpp/ContinuumSystem/output_ul3_$NEW_NAME.txt|" auto_generated_ul_3/$NEW_NAME.slurm
  sed -i "s|mpirun ./build_ul/ContinuumSystem params/ul_cluster.config|mpirun ./build_ul/ContinuumSystem auto_generated_ul_3/$NEW_NAME.config|" auto_generated_ul_3/$NEW_NAME.slurm

  # Execute the program
  sbatch auto_generated_ul_3/$NEW_NAME.slurm
done

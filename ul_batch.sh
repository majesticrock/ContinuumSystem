#!/bin/bash
input_file="params/for_auto.txt"
readarray -t NEW_VALUES < "${input_file}"

declare -A TOKENS
TOKENS=(
  [1]="phonon_coupling"
  [g]="phonon_coupling"
  [2]="omega_debye"
  [w]="omega_debye"
  [3]="screening"
  [l]="screening"
)

echo "Select the parameter that is to be varied:"
echo "1 / g: phonon_coupling"
echo "2 / w: omega_debye"
echo "3 / l: screening"
read -p "Enter your choice: " choice
TOKEN=${TOKENS[$choice]}

if [ -z "$TOKEN" ]; then
  echo "Invalid choice. Exiting."
  exit 1
fi

CURRENT_TIME=$(date +"%Y%m%d_%H%M%S")

rm -rf auto_generated_ul_${CURRENT_TIME}/
mkdir -p auto_generated_ul_${CURRENT_TIME}

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
      sed "s/$TOKEN_NAME $TOKEN_VALUE/$TOKEN_NAME $NEW_VALUE/" params/ul_cluster.config > auto_generated_ul_${CURRENT_TIME}/$NEW_NAME.config
      break
    fi
  done < params/ul_cluster.config
  cp slurm/ul_modes.slurm auto_generated_ul_${CURRENT_TIME}/$NEW_NAME.slurm
  sed -i "s|#SBATCH --job-name=modes|#SBATCH --job-name=${CURRENT_TIME}_$NEW_NAME|" auto_generated_ul_${CURRENT_TIME}/$NEW_NAME.slurm
  sed -i "s|#SBATCH --output=/home/althueser/phd/cpp/ContinuumSystem/modes_output.txt|#SBATCH --output=/home/althueser/phd/cpp/ContinuumSystem/ul_${CURRENT_TIME}_output_$NEW_NAME.txt|" auto_generated_ul_${CURRENT_TIME}/$NEW_NAME.slurm
  sed -i "s|mpirun ./build_cluster/ContinuumSystem params/ul_cluster.config|mpirun ./build_cluster/ContinuumSystem auto_generated_ul_${CURRENT_TIME}/$NEW_NAME.config|" auto_generated_ul_${CURRENT_TIME}/$NEW_NAME.slurm

  # Execute the program
  sbatch auto_generated_ul_${CURRENT_TIME}/$NEW_NAME.slurm
done

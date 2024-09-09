#!/bin/bash --login

# Variables
# Simulation script
# Amber topology file
# Amber coordinate file

## Setting up filepaths, working directories ##

master_path="$(pwd)"
simulation_script="$1"
system_name="$2"
job_name="${system_name}.job"

output_file="${system_name}.out"

echo $master_path
echo $simulation_script
echo $system_name
echo $output_file

## SLURM directives ##

cat > $job_name << EOF
#!/bin/bash --login
#SBATCH -o openmm_%J.out
#SBATCH -e openmm_%J.err
#SBATCH --job-name=openmm
#SBATCH --account=scw1976
#SBATCH -p s_compute_chem
#SBATCH --ntasks=10
#SBATCH --tasks-per-node=10
#SBATCH --mem-per-cpu=4000
#SBATCH --time=3-00:00

## Module loading ##

module purge
module load mamba/4.14.0-0
source activate AmberTools23

## Job info ##

echo "OPENMM job $SLURM_JOBID"
echo "INPUT STYSTEM $system_name"
echo "SIMULATION_SCRIPT $simulation_script"

## Working directory setup ##

python3 "$simulation_script" "$system_name" > "$output_file"
EOF
## Submit the job ##

chmod +x $job_name
sbatch "$job_name"
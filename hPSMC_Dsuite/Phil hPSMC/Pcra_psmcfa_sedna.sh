#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=psmcfa   ## job name
#SBATCH -e psmcfa_%j.e.txt    ## error message name
#SBATCH -o psmcfa.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 5    ## (cpu cores per task)
#SBATCH -p standard  #partitians= standard= 96G; medmem=192G; himem = 1.5TB
#SBATCH --mem=25G # max memory
#SBATCH -t 1-0   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1                   ## Run a single task
#SBATCH --array=1-3 # The total number of possible pairwise comparisons from a list of options is n(n-1)/2
#  SBATCH --array=0-$(($(wc -l < BW_consensus_genomes_list.txt) - 1))  # Array range (according to ChatGPT)
#  SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################

module load bio/angsd/0.940

psmcfa_script=/home/pmorin/scripts/misc/psmcfa_from_2_fastas.py
# The python script for this step can be downloaded from the Cahill et al. 2016 repository (https://github.com/jacahill/hPSMC/tree/master)
######################### set up python virtual environment:###############


ProjDir=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/hPSMC
CONSDIR=${ProjDir}/hPSMC_consensus_sequences
CONSLIST=${CONSDIR}/Pcra_consensus_autosome_genomes_list.txt
# note, consensus genome must be decompressed
psmcfaDIR=${ProjDir}/psmcfa_files

mkdir -p ${psmcfaDIR}

# combine consensus sequences from two species to create a pseudodiploid
# run script for all combinations of genomes in consensus genomes list (scripted with ChatGPT)
# Read file list into an array
mapfile -t file_list < ${CONSLIST}

# Get the total number of files
num_files=${#file_list[@]}

# Dynamically calculate the total number of pairs
total_pairs=$((num_files * (num_files - 1) / 2))

# Validate SLURM_ARRAY_TASK_ID (-le='less or equal'; -gt="greater than"; -ge="greater than or equal")
if [[ $SLURM_ARRAY_TASK_ID -gt $total_pairs ]]; then
    echo "Invalid SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Adjust SLURM_ARRAY_TASK_ID to start from 0
adjusted_task_id=$((SLURM_ARRAY_TASK_ID - 1))

pair_index=0
for ((file1_index=0; file1_index<num_files-1; file1_index++)); do
    for ((file2_index=file1_index+1; file2_index<num_files; file2_index++)); do
        if [[ $pair_index -eq $adjusted_task_id ]]; then
            file1=${file_list[$file1_index]}
            file2=${file_list[$file2_index]}

# get sample names from filenames
ID1=`echo $file1 | cut -f2 -d "_"`
ID2=`echo $file2 | cut -f2 -d "_"`


# Run your Python script with the file pair
source /home/pmorin/programs/python/bin/activate

python3 ${psmcfa_script} -b10 -m5 ${CONSDIR}/"$file1" ${CONSDIR}/"$file2" > ${psmcfaDIR}/hPSMC_${ID1}_${ID2}.psmcfa
            exit 0
        fi
        pair_index=$((pair_index + 1))
    done
done

# If no valid pair was found (should not happen)
echo "Error: Pair not found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
exit 1

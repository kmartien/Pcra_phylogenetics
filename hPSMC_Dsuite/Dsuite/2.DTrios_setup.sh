#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=Dsuite   ## job name
#SBATCH -e Dsuite%j.e.txt    ## error message name
#SBATCH -o Dsuite.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 1    ## (cpu cores per task)
#SBATCH -p medmem  #partitians= standard= 96G; medmem=192G; himem = 1.5TB
#  SBATCH --mem=25G # max memory
#SBATCH -t 1:00:00   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1                   ## Run a single task
#  SBATCH --array=1-3 # 
#  SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################

## ===================================
#  Variables you need to change before running!
#    1. VCF_DIR: Pathway to directory containing your split per-chromosome VCF files
#    2. SPECIESTREE: Pathway to species tree
#    3. SPECIES_SET: Pathway to species set file
#    4. OUT_DIR: Pathway to desired output directory
#    5. DSUITE_PATH: Pathway to Dsuite Build folder containing Dsuite executable
## ===================================

# User Defined Variables
VCF_DIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite/vcf_files/split-vcf
SPECIES_SET=${VCF_DIR}/"Pcra_species_sets.txt"
SPECIES_TREE=${VCF_DIR}/Pcra_Gmac_Oorc_129_largest_Busco_rearranged.newick
OUT_DIR=${VCF_DIR}/DTrios_output
#  DSUITE_PATH="~/bin/Dsuite/Build"
DSUITE_PATH=/home/pmorin/programs/Dsuite/Build

mkdir -p ${OUT_DIR}

# Output file with all the commands
COMMANDS_FILE="dtrios_cmds.txt"
JOB_SCRIPT="run_dtrios_array.sh"

# Clear any existing commands
> "$COMMANDS_FILE"

# Generate Dsuite Dtrios commands + output to command file for all chromosome vcf files in a directory (chromosomes extracted from whole genome vcf file using script 1.split_vcf_by_chrom)
for vcf in $VCF_DIR/*renamed.vcf.gz; do
    output_prefix=$(basename "$vcf" .vcf.gz)
    echo "$DSUITE_PATH/Dsuite Dtrios -o $OUT_DIR/$output_prefix --ABBAclustering -t $SPECIES_TREE $vcf $SPECIES_SET" >> "$COMMANDS_FILE"
done

# Count number of lines/commands
NUM_COMMANDS=$(wc -l < "$COMMANDS_FILE")

# Create SLURM job script
cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=Dsuite   ## job name
#SBATCH -e Dsuite%j.e.txt    ## error message name
#SBATCH -o Dsuite.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -p medmem
#SBATCH --time=4-0
#SBATCH --array=1-${NUM_COMMANDS}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40G
#SBATCH --output=$OUT_DIR/logs/Dtrios_%A_%a.out
#SBATCH --error=$OUT_DIR/logs/Dtrios_%A_%a.err

# Load necessary modules or conda environments
# module load dsuite
dsuite=/home/pmorin/programs/Dsuite
# source activate dtrios_env (python)
source /home/pmorin/programs/python/bin/activate 


# Get the command corresponding to this task ID
CMD=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $COMMANDS_FILE)

# Run the command
echo "Running: \$CMD"
eval \$CMD

EOF

# Make sure log directory exists
mkdir -p $OUT_DIR/logs

# Submit the job
sbatch "$JOB_SCRIPT"

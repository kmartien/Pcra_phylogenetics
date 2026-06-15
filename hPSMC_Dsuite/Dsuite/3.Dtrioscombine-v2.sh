#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=DsuiteCombine   ## job name
#SBATCH -e DsuiteCombine%j.e.txt    ## error message name
#SBATCH -o DsuiteCombine.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 10    ## (cpu cores per task)
#SBATCH -p medmem  #partitians= standard= 96G; medmem=192G; himem = 1.5TB
#  SBATCH --mem=25G # max memory
#SBATCH -t 1:00:00   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1                   ## Run a single task
#  SBATCH --array=1-3 # 
#  SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################

# User Defined Variables
VCF_DIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite/vcf_files/split-vcf
SPECIES_TREE=${VCF_DIR}/Pcra_Gmac_Oorc_129_largest_Busco_rearranged.newick
PREFIX="Pcra_dsuite_results"
OUT_DIR=${VCF_DIR}/DTrios_output
DSUITE_PATH="/home/pmorin/programs/Dsuite/Build"
Dtrios_env="/home/pmorin/programs/python/bin/activate"

# Output file with all the commands
COMMANDS_FILE="dtrios_cmds.txt"
JOB_SCRIPT="run_dtrioscombine.sh"


# Load necessary modules or conda environments
# module load dsuite
DSUITE_PATH=${DSUITE_PATH}

# Collects file prefixes, but removes the "_combine.txt" suffix that Dsuite adds to the output files
# and stores them in a space-separated string
file_prefixes=""
for f in "$OUT_DIR"/*_combine.txt
do 
    abs_path=$(realpath "$f") 
    trimmed="${abs_path%_combine.txt}" 
    file_prefixes+="$trimmed " 
done

# source activate dtrios_env
source ${Dtrios_env}

# Run the command
echo "Running: \$DSUITE_PATH/Dsuite DtriosCombine -o $PREFIX -t $SPECIES_TREE $file_prefixes"
$DSUITE_PATH/Dsuite DtriosCombine -o $PREFIX -t $SPECIES_TREE $file_prefixes

deactivate

# Submit the job
# sbatch "$JOB_SCRIPT"

###########
# Usage: Dsuite DtriosCombine [OPTIONS] DminFile1 DminFile2 DminFile3 ....
# Combine the BBAA, ABBA, and BABA counts from multiple files (e.g per-chromosome) and output the overall D stats,
# p-values and f4-ratio values
# 
# 	   -h, --help                              display this help and exit
# 	   -o, --out-prefix=OUT_FILE_PREFIX        (optional) the prefix for the files where the results should be written
# 											   output will be put in OUT_FILE_PREFIX_combined_BBAA.txt, OUT_FILE_PREFIX_combined_Dmin.txt, OUT_FILE_PREFIX_combined_tree.txt etc.
# 											   by default, the prefix is "out"
# 	   -n, --run-name                          (optional) run-name will be included in the output file name after the PREFIX
# 	   -t , --tree=TREE_FILE.nwk               (optional) a file with a tree in the newick format specifying the relationships between populations/species
# 											   D and f4-ratio values for trios arranged according to the tree will be output in a file with _tree.txt suffix
# 	   -s , --subset=start,length              (optional) only process a subset of the trios

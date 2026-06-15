#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=DsuiteFigs   ## job name
#SBATCH -e DsuiteFigs%j.e.txt    ## error message name
#SBATCH -o DsuiteFigs.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 1    ## (cpu cores per task)
#SBATCH -p standard  #partitians= standard= 96G; medmem=192G; himem = 1.5TB
#  SBATCH --mem=25G # max memory
#SBATCH -t 0:30:00   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1                   ## Run a single task
#  SBATCH --array=1-3 # 
#  SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################

# User Defined Variables
VCF_DIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite/vcf_files/split-vcf
SPECIES_SET=${VCF_DIR}/"Pcra_species_sets.txt"
PREFIX="Pcra_dsuite_results"
OUT_DIR=${VCF_DIR}/DTrios_output

# Output file with all the commands
COMMANDS_FILE="figures_cmds.txt"
JOB_SCRIPT="run_figures.sh"
Dsuite_scripts=/home/pmorin/programs/Dsuite
# from https://github.com/mmatschiner/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/


## Load necessary modules or conda environments
# module purge # this is a programming language; already on Sedna?
# module load GCCcore/12.2.0 # gcc/8.3.1 is on Sedna
# module load Ruby/3.2.2

# Generate plot_order file from species_set.txt file
cut -f 2 $SPECIES_SET | uniq | head -n 1000 > $OUT_DIR/plot_order.txt

# Generate D-statistic plot(D_MAX sets the upper-limit of the color scale for the plot - may need to alter based on your data)
D_MAX=0.2
ruby ${Dsuite_scripts}/plot_d.rb ${PREFIX}_combined_BBAA.txt $OUT_DIR/plot_order.txt $D_MAX ${PREFIX}_combined_BBAA_D.svg

# Generate f-4 ratio plot (F4_MAX sets the upper-limit of the color scale for the plot - may need to alter based on your data)
F4_MAX=0.1
ruby ${Dsuite_scripts}/plot_f4ratio.rb ${PREFIX}_combined_BBAA.txt $OUT_DIR/plot_order.txt $F4_MAX ${PREFIX}_combined_BBAA_f4ratio.svg

# Submit the job
# sbatch "$JOB_SCRIPT"



#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=DsuiteFbranch   ## job name
#SBATCH -e DsuiteFbranch%j.e.txt    ## error message name
#SBATCH -o DsuiteFbranch.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 1    ## (cpu cores per task)
#SBATCH -p standard  #partitians= standard= 96G; medmem=192G; himem = 1.5TB
#  SBATCH --mem=25G # max memory
#SBATCH -t 0:30:00   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1                   ## Run a single task
#  SBATCH --array=1-3 # 
#  SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>
######################################################################################

# User Defined Variables
VCF_DIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/Dsuite/vcf_files/split-vcf
SPECIES_TREE=${VCF_DIR}/Pcra_Gmac_Oorc_129_largest_Busco_rearranged.newick
SPECIES_SET=${VCF_DIR}/"Pcra_species_sets.txt"
PREFIX="Pcra_dsuite_results"
ZSCORES=$VCF_DIR/Ziphiidae_Dsuite_Zscores.txt
OUT_DIR=${VCF_DIR}/DTrios_output
DSUITE_PATH="/home/pmorin/programs/Dsuite/Build"
DTOOLS_PATH="/home/pmorin/programs/Dsuite/utils/dtools.py"
DTOOLS_PATH_z="/home/pmorin/programs/Dsuite/utils/dtools-zscore-v2.py"
DTOOLS_print="/home/pmorin/programs/Dsuite/utils/plot_zscore_pdf-v2.py"

# Output file with all the commands
# COMMANDS_FILE="fbranch_cmds.txt"
# JOB_SCRIPT="run_fbranch.sh"


# 1. Calulate f-branch statistic

$DSUITE_PATH/Dsuite Fbranch -p 0.001 --Zb-matrix $SPECIES_TREE ${PREFIX}_combined_tree.txt > ${PREFIX}_combined_Fbranch_wZscore.txt
# or, without z-scores:
$DSUITE_PATH/Dsuite Fbranch -p 0.001 $SPECIES_TREE ${PREFIX}_combined_tree.txt > ${PREFIX}_combined_Fbranch_noZscore.txt


# -p, --pthresh (default=0.01) fb scores whose associated p-value is less than 
# --Zb-matrix: (optional)  output the equivalent of fb-statistic, but with Z-scores to assess statistical significance. Add only to get combined_Fbranch.txt file, then remove from script and re-run before running python script (won't work otherwise).


# 2. Generate Zscore distribution plot with % cut-off values, and f-branch figure

source /home/pmorin/programs/python/bin/activate

### generate zscore distribution with cut-off percent

# usage: plot_zscore_pdf.py [-h] fbranch percentile
python3 ${DTOOLS_print} ${PREFIX}_combined_Fbranch_wZscore.txt 20.0 --output Z-score-plot_noZscore_20pct.pdf

	#  positional arguments:
	#   fbranch          Fbranch output file with Z-scores (-Z option)
	#   percentile     Top percentile to highlight in dataset
	#  optional arguments:
	#   -h, --help       Show this help message and exit
	#   --output         Output PNG file (default: Z-score-PDF.png) (default: Z-score-PDF.png)


### generate fbranch plot

# with z-score numbers within boxes (requires combined_Fbranch.txt with Zscores (--Zb-matrix option) in (1) above)
python3 $DTOOLS_PATH_z --zscore 15.0 --ladderize --use_distances --tree-label-size 7 --outgroup $(grep "Outgroup" $SPECIES_SET | cut -f1) -n fbranch ${PREFIX}_combined_Fbranch_wZscore.txt $SPECIES_TREE

# without z-score numbers: (requires combined_Fbranch.txt WITHOUT Zscores (--Zb-matrix option) in (1) above)
python3 $DTOOLS_PATH --ladderize --use_distances --tree-label-size 7 --outgroup $(grep "Outgroup" $SPECIES_SET | cut -f1) -n fbranch ${PREFIX}_combined_Fbranch_noZscore.txt $SPECIES_TREE

# --zscore 10.0 = minimum Zscore to show in plot (if --Zb-matrix used to generate Fbranch.txt in step 1 above).

deactivate

EOF

# Submit the job
# sbatch "$JOB_SCRIPT"

#####
# optional arguments for dtools.py:
#   -h, --help            show this help message and exit
#   -n RUN_NAME, --run-name RUN_NAME
#                         Base file name for output plots. (default: fbranch)
#   --outgroup OUTGROUP   Outgroup name in newick file. (default: Outgroup)
#   --use_distances       Use actual node distances from newick file when
#                         plotting tree. (default: False)
#   --ladderize           Ladderize the input tree before plotting. (default:
#                         False)
#   --color-cutoff COLOR_CUTOFF
#                         Set the darkest red to this f_branch value. (default:
#                         None)
#   --tree-label-size TREE_LABEL_SIZE
#                         Set the font size of the tree leaf names. (default:
#                         14)
#   --dpi DPI             Set the dpi for the output .png. (default: 150)

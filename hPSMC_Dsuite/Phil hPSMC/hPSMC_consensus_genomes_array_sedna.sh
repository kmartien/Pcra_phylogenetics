#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=consensus_gen   ## job name
#SBATCH -e consensus_gen%j.e.txt    ## error message name
#SBATCH -o consensus_gen.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE,BEGIN,END,FAIL,ALL)
#SBATCH -c 5    ## (cpu cores per task)
#SBATCH -p medmem  #partitians= standard= 96G; medmem=192G; himem = 1.5TB
# SBATCH --mem=40G # max memory
#SBATCH -t 2-0   ## walltime in mins,mins:secs,hrs:mins:secs,days-hours 
#SBATCH --array=1-3 # The total number of possible pairwise comparisons from a list of #SBATCH --ntasks=1                   ## Run a single task (from Morgan's script)
# SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################

module load bio/angsd/0.940

ProjDir=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025
OUTDIR=${ProjDir}/hPSMC_consensus_sequences
# error=log_error_files
# hPSMC_cons_genomes_method_050325
BAMLIST=Pcra_bamlist.txt
THREADS=5

mkdir -p ${OUTDIR}
mkdir -p ${error}

# create consensus sequences from BAM files (#-doFasta 2; call consensus base)
# Generate the list of files using ls command

readarray -t files < ${ProjDir}/${BAMLIST}

# Get the number of files
# num_files=${#files[@]}
cd ${ProjDir}
# Get the file for this array task
current_file=${files[$SLURM_ARRAY_TASK_ID-1]}
file=`echo $current_file | cut -f5 -d "/"`
echo ${file}
sp=`echo $file | cut -f1 -d "_"`
echo ${sp}
ID=`echo ${file} | cut -f2 -d "_"`
echo ${ID}

# List of chromosomes to include (space-separated) (from reference genome)
# check chr names in bam using samtools view -H Pcra_z0018462_dedup_noRepeats.bam | head -25
CHROMOSOMES=("NC_090296.1" "NC_090297.1" "NC_090298.1" "NC_090299.1" "NC_090300.1" "NC_090301.1" "NC_090302.1" "NC_090303.1" "NC_090304.1" "NC_090305.1" "NC_090306.1" "NC_090307.1" "NC_090308.1" "NC_090309.1" "NC_090310.1" "NC_090311.1" "NC_090312.1" "NC_090313.1" "NC_090314.1" "NC_090315.1" "NC_090316.1")

# Convert array to ANGSD format (comma-separated)
CHROM_LIST=$(IFS=,; echo "${CHROMOSOMES[*]}")

angsd -i ${current_file} -nThreads ${THREADS} -doFasta 2 -minQ 25 -minmapq 25 -uniqueonly 1 -remove_bads 1 -setMinDepthInd 5 -doCounts 1 -rf <(echo -e "${CHROMOSOMES[*]}" | tr ' ' '\n') -out ${OUTDIR}/${sp}_${ID}_hPSMC_autosome_consensus

# angsd -i ${current_file} -nThreads ${THREADS} -doFasta 2 -r "CM046080.1:1-247284360 CM046081.1:1-193415349 CM046082.1:1-185486766 CM046083.1:1-184908288 CM046084.1:1-155519127 CM046085.1:1-136865726 CM046086.1:1-125805690 CM046087.1:1-117335195 CM046088.1:1-114113211 CM046089.1:1-111292810 CM046090.1:1-106130936 CM046091.1:1-104700920 CM046092.1:1-93362332 CM046093.1:1-91985674 CM046094.1:1-91196782 CM046095.1:1-88733150 CM046096.1:1-78597266 CM046097.1:1-67049709 CM046098.1:1-66567213 CM046099.1:1-45032744" -minQ 25 -minmapq 25 -uniqueonly 1 -remove_bads 1 -setMinDepthInd 5 -doCounts 1 -out ${OUTDIR}/${sp}_${ID}_hPSMC_autosome_consensus


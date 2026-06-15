#! /bin/bash

#SBATCH --job-name=hPSMC
#SBATCH -e hPSMC_%j.e.txt
#SBATCH -o hPSMC_%j.log 
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  # (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 10    # <number of cores to ask for>
#SBATCH --mem=40G   # job memory request
#SBATCH -p standard
#SBATCH --array=1-3
#SBATCH -t 3-0   # walltime in mins or mins:secs or hrs:mins:secs. 
#SBATCH -D /scratch/pmorin/temp # <folder to change to when starting the job>
######################
#programs
module load bio/samtools/1.11
module load bio/bcftools/1.11
module load tools/gnuplot/5.4.1
PSMC=~/programs/psmc

set -eux
##############################################################

INDIR=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/hPSMC
psmcfa_DIR=${INDIR}/psmcfa_files
OUTDIR=${INDIR}/hPSMC_output
mkdir -p ${OUTDIR}
psmcfa_files=Pcra_psmcfa_files.txt # list of psmcfa files; e.g., ls -1 *.psmcfa > psmcfa_files.txt
THREADS=10

# PSMC parameters
TEMP_DIR=/scratch/pmorin/temp
MUT=4.90E-10 #mutation rate (µ/site/yr)
# Dornburg 2012 estimate for odontocetes (nuDNA) was 9.1x10^-10 sub/site/yr. Robinson et al. (vaquita paper) estimate was 4.90x10^-10 sub/site/yr.
GEN=24 # based on pilot whale. from Taylor et al. 2007, Table 1, T(r=0)=generation length under pre-disturbance conditions. 

MUTRATE=$(echo $GEN $MUT | awk '{ printf "%.12f", $1*$2 }') # mutation rate (µ/site/yr) multiplied by generation time. 

 #vaquita estimated rate = 5.83E-9/site/gen (5.83E-9/11.9yr/gen = 4.90E-10 µ/site/yr). Mutation rate for Ziphiids with ~20yr/gen = 9.80E-9 µ/site/gen.
MUT2=4.90E10_msy	# msy=mutations per site per year for printing mutation rate on plot and in file name

t=15	# default from Li and Durbin is 15
PSMC_INT="4+25*2+4+6" # default from Li and Durbin is 4+25*2+4+6; changed to 1+1+1+1+25*2+4+6 based on Hilgers et al. 2025 (in press)
#########################################################################
#########################################################################
#########################################################################
# Generate the list of files using ls command
readarray -t files < ${psmcfa_DIR}/${psmcfa_files}

# Get the number of files
# num_files=${#files[@]}
cd ${psmcfa_DIR}
# Get the file for this array task
current_file=${files[$SLURM_ARRAY_TASK_ID-1]}
short1=`echo $current_file | cut -f1 -d "."` # everything before "."
short_file=`echo $short1 | cut -f2,3 -d "_"` # two species names when filenames like "hPSMC_Mgra_Mste.psmcfa"

#########################################################################
# PSMC
ID=${short_file}

#generate the psmc file using the default settings for humans (-N25, -t20, -r5).
${PSMC}/psmc -N25 -t${t} -r5 -p ${PSMC_INT} -o ${OUTDIR}/${ID}_${MUT2}_t${t}.psmc ${psmcfa_DIR}/$current_file

#Make psmc plots and adapt the scaling using this psmc file.
nice ${PSMC}/utils/psmc_plot.pl -u ${MUTRATE} -g ${GEN} -s 10 -RM "" ${OUTDIR}/${ID}_${MUT2}_t${t}_psmc.out ${OUTDIR}/${ID}_${MUT2}_t${t}.psmc

######################
# bootstrap PSMC

# split the PSMC file
# ${PSMC}/utils/splitfa ${INDIR}/$current_file > ${OUTDIR}/split_${ID}_diploid.psmcfa
# 
# # PSMC bootstrap, multithread = 12 (-P 12)
# seq 100 | xargs -P ${THREADS} -i ${PSMC}/psmc -N25 -t${t} -r5 -b -p ${PSMC_INT} -o ${TEMP_DIR}/${ID}_round-{}.psmc ${OUTDIR}/split_${ID}_diploid.psmcfa | sh
# 
# # cat original.psmc round-*.psmc > combined.psmc
# # within each folder containing original sample psmc and bootstrap files
# cat ${OUTDIR}/${ID}_${MUT}_t${t}.psmc ${TEMP_DIR}/${ID}_round-*.psmc > ${OUTDIR}/merged_${ID}_t${t}_boot.psmc
# 
# # And then plot it: (-RM -> Id placed on plot).
# ${PSMC}/utils/psmc_plot.pl -u ${MUTRATE} -g ${GEN} -RM "" ${OUTDIR}/merged_${ID}_t${t}_boot.out ${OUTDIR}/merged_${ID}_t${t}_boot.psmc
# 
# cp ${INDIR}/*PSMCboot_SEDNA.sh ${OUTDIR}
########################################################################################
# for PSMC parameter info, see https://github.com/lh3/psmc
# the `-p' option specifies that there are 64 atomic time intervals and 28 (=1+25+1+1)
# free interval parameters. The first parameter spans the first 4 atomic time
# intervals, each of the next 25 parameters spans 2 intervals, the 27th spans 4
# intervals and the last parameter spans the last 6 time intervals. The `-p' and
# `-t' options are manually chosen such that after 20 rounds of iterations, at
# least ~10 recombinations are inferred to occur in the intervals each parameter
# spans.

### use the custom script provide in the utils folder with PSMC to convert the fastq file to a pseudo-fasta file. It basically codes regions in 100bp bins as being heterozygote or homozygote.


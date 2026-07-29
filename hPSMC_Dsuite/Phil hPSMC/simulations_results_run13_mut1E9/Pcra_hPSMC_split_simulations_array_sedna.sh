#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=hPSMC_sim   ## job name
#SBATCH -e hPSMC_sim.%j.e.txt    ## error message name
#SBATCH -o hPSMC_sim.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 2    ## (cpu cores per task)
#SBATCH -p medmem  # partitians= standard= 96G; medmem=192G; himem = 1.5TB
# SBATCH --mem=20G # max memory
#SBATCH -t 72:00:00   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1                   ## Run a single task
#SBATCH --array=1-3 # The total number of possible pairwise comparisons from a list of options is n(n-1)/2
# SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################

module load tools/gnuplot/5.4.1
module load bio/ms/201810

scriptdir=/home/pmorin/scripts/misc/
script=${scriptdir}hPSMC_quantify_split_time.py # from https://github.com/jacahill/hPSMC/tree/master
PSMC=~/programs/psmc

# Parameters:
ProjDir=/home/pmorin/projects/Pcra/hPSMC_Dsuite_June2025/hPSMC
INDIR=${ProjDir}/hPSMC_output
# PSMClist=Mesoplodon_hPSMC_psmc_files_list.txt # e.g., ls -1 *t15.psmc > hPSMC_psmc_files_list.txt
# preDivNe_files=Mesoplodon_hPSMC_preDiv_Ne.txt # Ne in integer (actual inferred Ne); same order as 
PSMClist=Pcra_PSMC_PreDivNe_list.csv # comma separated list of pseudodiploid psmc files and the pre-divergence Ne from the empirical hPSMC analysis. Value = 1E4*psmc value (e.g., psmc=5.2291, use 52291).
# ensure that csv file doesn't contain dos characters
dos2unix ${INDIR}/${PSMClist}

# PreDivNe=30000 					# pre-divergence Ne
Range_lower=0				# divergence time range lower limit; default=0
Range_upper=300000				# divergence time range upper limit; default=10000000
Sim_num=7						# number of simulations to run 
THREADS=5						# number of simulations to run in parallel, 1/cpu
MUT=1.0E-9 #mutation rate (µ/site/yr)
MUT2=1.00E9_msy
# Sim_num is calculated as upper-lower/50+1, e.g. ((1500000-50000)/50000)+1=30 (simulation for every 50,000yr in the range)
# ((300000-10000)/10000)+1=30 (simulation for every 10,000yr in the range)
# MUT: Dornburg 2012 estimate for odontocetes (nuDNA) was 9.1x10^-10 sub/site/yr. Robinson et al. (vaquita paper) estimate was 4.90x10^-10 sub/site/yr (vaquita estimated rate = 5.83E-9/site/gen (5.83E-9/11.9yr/gen = 4.90E-10 µ/site/yr)). Van Cise et al. 2025? ~1.9E-10 for odontocetes; ~1.75E-10 mysticetes
GEN=24
MUTRATE=$(echo $GEN $MUT | awk '{ printf "%.12f", $1*$2 }') # mutation rate (µ/site/yr) multiplied by generation time = µ/site/gen (used in PSMC)

######################### set up python virtual environment:###############

# Get working psmc file based on array number
NUM=$(printf %02d ${SLURM_ARRAY_TASK_ID})

file=$(head -n ${NUM} ${INDIR}/${PSMClist} | tail -n 1)
echo ${file}

# get filename
psmcfile=$(echo "$file" | cut -d ',' -f 1)
echo "$psmcfile"

Ne=$(echo "$file" | cut -d ',' -f 2)
echo "$Ne"

# get pre-divergence Ne from ${ProjDir}/${preDivNe}
# Ne=$(head -n ${NUM} ${ProjDir}/${preDivNe_files} | tail -n 1)
# echo $Ne
# PreDivNe=$(echo "scale=2; ${Ne}/10000" | bc) # converts Ne to Ne/10000 for python script.
#echo ${PreDivNe}

# filesname format = ${species1}_${species2}_4.90E10_msy_t15.psmc
species1=`echo $psmcfile | cut -f1 -d "_"`
echo $species1
species2=`echo $psmcfile | cut -f2 -d "_"`
echo $species2

OUTDIR=${ProjDir}/${species1}_${species2}_mut${MUT}_simulations
mkdir -p ${OUTDIR}

# empirical psmc output file
emp=${species1}_${species2}_${MUT2}_t15.psmc
echo $emp

# Run Python script with the file pair
source /home/pmorin/programs/python/bin/activate

# write bash script to run simulations
python ${script} -H ${scriptdir} -P ${PSMC}/psmc -N ${Ne} -l ${Range_lower} -u ${Range_upper} -s ${Sim_num} -p ${THREADS} -o ${OUTDIR}/hPSMC_sim_ ${INDIR}/${emp} > ${OUTDIR}/${species1}_${species2}_simulations.sh

# execute bash file to run simulations
bash ${OUTDIR}/${species1}_${species2}_simulations.sh

# -N is the predivergence Ne (see below; default 10,000yr)
# -l = lower limit of time range around divergence
# -u = upper limit of time range around divergence
# -s = number of simulations for 1/50,000 years

#From this output you will get a series of of .psmc files for each simulated interval (e.g. 200000, 250000, 300000 ....) and you plot the psmc file as you did before. 
# species1.species2.300000.ms_sim.psmc, etc.

#########################################################################################
# generate plot data snd plot simulated hPSMC for each simulation.
cd ${OUTDIR}
for file2 in *.ms_sim.psmc 
do 
nice ${PSMC}/utils/psmc_plot.pl -u ${MUTRATE} -g ${GEN} -s 10 -RM ${species1}"_"${species2} ${OUTDIR}/${species1}_${species2}_${file2}.out ${file2}
done 

# modify plot output (.psmc.out.0.txt) files to add the time block to the data
for file in *out.0.txt 
do 
echo ${file}
ID=`echo ${file} | cut -f1 -d "."`
echo ${ID}
block=`echo ${ID} | cut -f5 -d "_"`
echo ${block}
awk '{print $1,$2,$3='${block}'}' ${file} > ${OUTDIR}/${file}.plot.me.txt  
done 

#Combine all the runs for one sample
cat ${OUTDIR}/*.plot.me.txt > ${OUTDIR}/${species1}_${species2}_mut${MUT}_simulations.combined.txt

# add empirical hpsmc output to simulations.combined.txt file
for file in ${INDIR}/${species1}_${species2}*psmc.out.0.txt
do 
echo ${file}
awk '{print $1,$2,$3="hpsmc"}' ${file} > ${OUTDIR}/${species1}_${species2}_hpsmc.plot.me.txt
done 

cat ${OUTDIR}/${species1}_${species2}_hpsmc.plot.me.txt >> ${OUTDIR}/${species1}_${species2}_mut${MUT}_simulations.combined.txt

# use R script "hPSMC_plot_combined_simulations.R"

######## determining input parameters:
# 1) pre-divergence Ne from PSMC file
# The third column is the Ne, so in this case, I would choose 1.50 (which is 15,000 because this column is Ne * 10,000)
# 0	3285755.73581677	0.006313	0.000000	0.000000
# 17136.4587150424	3285755.73581677	0.006838	0.000000	0.000000
# 35698.0339269068	8388165.25794203	0.002901	0.000000	0.000000
# 55800.2198813559	8388165.25794203	0.003142	0.000000	0.000000
# 77572.9476048729	10526950.380935	0.002712	0.000000	0.000000
# 101156.460110169	10526950.380935	0.002937	0.000000	0.000000
# 126697.187601695	3496482.36206071	0.009578	0.000000	0.000000
# 154362.184256356	3496482.36206071	0.010374	0.000000	0.000000
# 184324.691442797	28347.9963749137	1.385910	0.000002	0.000001
# 216778.574502119	28347.9963749137	1.501104	0.000002	0.000002
# 251925.948364407	1326.12676020286	34.754659	0.000055	0.000053
# 289995.739123941	1326.12676020286	37.641471	0.000060	0.000072
# 331229.247258475	108.914672995763	496.203763	0.000790	0.001049
# etc.
# 
#  	OPTIONS:
#  	-h print help message with user options
#  	-o OUT, --out=OUT     output directory for simulations and prefix all files for the run, default="./hPSMC_sim_"
#  	-N NE, --Ne=NE        The ancestral population size to simulate, default=10,000
#  	-l LOWER, --lower=LOWER		lower bound for simulations, the most recent divergence time to be simulated
#  	-u UPPER, --upper=UPPER		upper bound for simulations, the most ancient divergence time to be simulated
#  	-s SIMS, --sims=SIMS  the number of simulations to conduct, simulations will evenly split between high and low, minimum value=2, minimum meaningful value=3
	# one simulation per 50,000 years: To figure out this out take the difference between u and l (1500000 - 200000) and divide it by 50,000. Then add 1 to the output. In this case that is 26 +1 = 27.
#  	-p PAR, --parallel=PAR		Number of simulations to run simultaneously
#  	-P PSMC, --PSMC=PSMC  If the psmc executable is not in your path give it's location, default = "psmc"
#  	-m MS, --ms=MS        If the ms executable is not in your path give it's location, default = "ms"
#  	-H HPSMC, --hPSMC=HPSMC
#  	If the hPSMC directory is not in your path give it's location, NOTE:  Just the directory not the script.  default = "./"
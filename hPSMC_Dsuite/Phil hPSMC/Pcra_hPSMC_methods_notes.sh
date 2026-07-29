hPSMC methods notes

1. generate consensus genomes from bam files (deduplicated, repeat masked).
	a) hPSMC method = hPSMC_consensus_genomes_sedna.sh. Uses angsd, but different parameters than the method used for extracting consensus genomes for phylogenomics.
# angsd -i ${BAM} -nThreads ${THREADS} -doFasta 2 -minQ $baseQ -minmapq $mapQ -uniqueonly 1 -remove_bads 1 -setMinDepthInd 5 -doCounts 1 -out ${BW_CON}/${SP}_${sample}_aligned_sorted_unique 
	# baseQ=25
	# mapQ=25
hPSMC_Dsuite_June2025]$ sbatch hPSMC_consensus_genomes_array_sedna.sh
Submitted batch job 2312929, 6/6/25 (done, ~25min)# 3 SWFSC Pcra assemblies to Pcra ref genome
# autosomes only (chr 1-21), repeat-masked bam file

2. Generate the pseudo-diploid genomes.
	a) hPSMC method uses python script psmcfa_from_2_fastas.py from Cahill et al. 2016 repository (https://github.com/jacahill/hPSMC/tree/master). Bash script = 
sbatch Pcra_psmcfa_sedna.sh
Submitted batch job 2315389, 6/6/25 (done, <20min)

#######################################################

3) Run hPSMC script (array script for PSMC in parallel): hPSMC_array_SEDNA.sh
# requires mutation rate and taxa-specific generation time. 

sbatch Pcra_hPSMC_array_SEDNA.sh
Submitted batch job 2316482, 6/6/25 (done, 2.5 days)

4) run hPSMC simulations

# prediv-NE set to ~32000
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2341634, 6/12/25

# re-run with prediv-Ne set to ~18000
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2422147, 6/28/25

# re-run with prediv-Ne set to ~9000
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2527317, 7/11/25

# re-run with prediv-Ne set to ~4500
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2535700, 7/13/25

# re-run with prediv-Ne set to ~2200
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2573171, 7/18/25
# two failed because I left a decimal point in the pre-div Ne estimate file. The one that finished was in the right time range, but Ne was too low, so the Ne range didn't overlap.

# re-run with prediv-Ne set to ~3500
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2611707

# re-run with prediv-Ne set to ~4500
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2626839

# re-run with prediv-Ne set to ~4800
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2637061

# back to run 6 prediv-Ne's (9500-9700). The Ne's are too low, so that the range of Ne on the plots falls below the exponential increase part of the psmc plot. 
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2675891, 7/28/25

# the minimum pre-div Ne now overlapps, but the exponential increase of the PSMC is outside of the range of the simulations. 
# conclusion is that simulations don't work at this time scale unless the Ne is very small, not based on empirical PSMC analysis. 

##############
# Re-running original hPSMC to make sure the mutation rate is correct (4.90E-10 instead of 9.10E-10). After that, run modified simulation script that correctly grabs the empirical hPSMC file. 

sbatch Pcra_hPSMC_array_SEDNA.sh
Submitted batch job 2687664, 7/30/25

# simulations
# Pre-div Ne's ~35000, based on hPSMC (double check)
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2720565, 8/4/25
# results: back to original plot pattern, with simulations around 1M year, and empirical data ~100,000yr.

# simulations using mutation rate from Cahill et al, 1.0x10-9, and 7 intervals of 50kyr from 0 to 300,000yr. 
sbatch Pcra_hPSMC_split_simulations_array_sedna.sh
Submitted batch job 2736964, 8/6/25






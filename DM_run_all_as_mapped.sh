#!/bin/bash
#SBATCH -A snic2022-5-255
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 6-01:00:00
#SBATCH -J DM_mergedMM_after40


#this script is to generate DamMet F file for all samples mapped with Asian Elephants

cd /crex/proj/snic2022-6-144/nobackup/YUXIN/DamMet/raw_data_AsRef

#declare -a elephant_arr=("RG.E.maximus.T13.bam" "RG.E.maximus.T16.bam" "RG.E.maximus.T18.bam" 
#	"RG.L.africana.T11.bam" 
#	"RG.M.americanum.ERR2260503.bam" 
#	"RG.merged.E.maximus.T13_T13_18.bam" 
#	"RG.merged.M.primigenius.E467.L163.L164.P005.M6.bam"
#	"RG.M.primigenius.E467.bam" "RG.M.primigenius.L163.bam" "RG.M.primigenius.L164.bam" "RG.M.primigenius.M6.bam" "RG.M.primigenius.P005.bam"
#	"RG.P.antiquus.ERR2260504.bam")
declare -a elephant_arr=("RG.merged.M.primigenius.E467.L163.L164.P005.M6.bam")


declare -a chrom_arr 
mapfile -t chrom_arr < /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/chromsomes_list.txt

#mapfile -t chrom_arr < <(tail -n 8 /crex/proj/snic2022-6-144/nobackup/YUXIN/scripts/chromsomes_list.txt)

declare -a NCPGs_arr=(25 50 10 100)

for sample in "${elephant_arr[@]}"
do
#do a for loop for BAM files and its corresponding output prefix 
	BAM="/crex/proj/snic2022-6-144/nobackup/METHYLATION_GENOMES/$sample"
	# N positions to estimate deamination rates from each end
	MAX_POS=30
	# path to reference genome
	FASTA=/crex/proj/snic2022-6-144/nobackup/TOM/MAMMOTH/REF/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.fna
	# minimum mapping quality
	minMapQ=1
	# minimum base quality
	minBaseQ=1
	# global prior on the fraction of methylated CpG's
	M=0.75
	# CpGs to exclude due to low coverage
	# EXCLUDE= None
	# N cpgs to include per window for estimating $f$, usually varied from 10-100
	#NCPG=10

	for chrom in "${chrom_arr[@]}"
	do
	# specify which chromosome to analyze
		CHROM="$chrom"
		echo "$chrom"

		for NCPG in "${NCPGs_arr[@]}"
		do 
			
			#output = samplename + chromosome + nCPGparatmeter
			sample_name=${sample%%".bam"*}
			OUT=As_"$sample_name"_"$chrom"_"$NCPG"

			#echo $sample_name
			echo $OUT
			/home/yuxind/DamMet/DamMet estDEAM -b ${BAM} -r ${FASTA} -c ${CHROM} -q ${minMapQ} -Q ${minBaseQ} -P ${MAX_POS} -M ${M} -O ${OUT}  -skip -verbose

			/home/yuxind/DamMet/DamMet estF -b ${BAM} -r ${FASTA} -c ${CHROM} -q ${minMapQ} -Q ${minBaseQ} -P ${MAX_POS} -M ${M} -O ${OUT} -N ${NCPG} -verbose -skip_empty_cpg
		done
	done
done


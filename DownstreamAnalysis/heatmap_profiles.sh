#!/bin/bash
#SBATCH --time=02:30:00
#SBATCH --account=def-jdilwort
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=HEATMAP_PROFILES
#SBATCH --output=%x-%j.out

module load python/3.9 scipy-stack/2021a
source ngsENV/bin/activate

echo "\nComputing matrix...\n"
computeMatrix scale-regions --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
	--samplesLabel Sample_1 Sample_2 Sample_3 \
	-R Reference_Files/mm10/gencode.vM23.annotation.gtf \
       -S Sample_1.bw \
	Sample_2.bw \
	Sample_3.bw \
       --skipZeros \
       -o samples_matrix.gz \
       --outFileSortedRegions samples_regions.bed \
       --numberOfProcessors max

echo "\nPlotting heatmap...\n"
plotHeatmap -m samples_matrix.gz -out samples_heatmap.png --colorList 'white,#00BFC4' 'white,#00BFC4' 'white,#F8766D'

echo "\nPlotting profile...\n"
plotProfile -m samples_matrix.gz --perGroup -out samples_profile.png --colors '#00BFC4' '#00BFC4' '#F8766D'

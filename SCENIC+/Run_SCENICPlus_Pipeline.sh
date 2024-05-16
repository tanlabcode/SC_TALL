#!/bin/bash	
#SBATCH -c 64
#SBATCH --mem=700G
#SBATCH -t 7-0:0:0
#SBATCH --output=%x.%j.out

####conda activate ScenicPlus

snakemake --cores 64



#!/bin/bash	
#SBATCH -c 4
#SBATCH --mem=100G
#SBATCH -t 7-0:0:0
#SBATCH --output=%x.%j.out

module load BEDTools/2.30.0-GCC-11.3.0

REGION_BED="/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL/SCENICPlus/Data_SCENICplus/ATAC_Region_Names.bed"
GENOME_FASTA="/mnt/isilon/tan_lab/sussmanj/Single_Cell_Tools/ScenicPlus/Genome_Files/hg38.fa"
CHROMSIZES="/mnt/isilon/tan_lab/sussmanj/Single_Cell_Tools/ScenicPlus/Genome_Files/hg38.chrom.sizes"
DATABASE_PREFIX="T_ALL_40"
SCRIPT_DIR="/mnt/isilon/tan_lab/sussmanj/Single_Cell_Tools/ScenicPlus/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        hg38_T_ALL_1kb_bg_padding.fa \
        1000 \
        yes

#!/bin/bash	
#SBATCH -c 40
#SBATCH --mem=600G
#SBATCH -t 7-0:0:0
#SBATCH --output=%x.%j.out

####conda activate create_cistarget_databases

CBDIR="/mnt/isilon/tan_lab/sussmanj/Single_Cell_Tools/ScenicPlus/Motif_Collection/v10nr_clust_public/singletons"

#Change the FASTA file 
FASTA_FILE="/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL/SCENICPlus/Database_Files/hg38_T_ALL_1kb_bg_padding.fa"

MOTIF_LIST="/mnt/isilon/tan_lab/sussmanj/Single_Cell_Tools/ScenicPlus/Motif_Collection/v10nr_clust_public/motifs.txt"
SCRIPT_DIR="/mnt/isilon/tan_lab/sussmanj/Single_Cell_Tools/ScenicPlus/create_cisTarget_databases"
DATABASE_PREFIX="T_ALL"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 40
    


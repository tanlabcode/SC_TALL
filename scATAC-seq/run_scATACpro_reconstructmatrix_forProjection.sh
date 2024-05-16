#!/bin/bash

#sampleIDS=($(ls /mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/RawDataOrganize)) ## for all ETP-ALL samples

#sampleIDS=(T_ALL_PASZMC_scATAC T_ALL_PAUHWY_scATAC T_ALL_PAUUVD_scATAC T_ALL_PAVMWF_scATAC T_ALL_PAWGWD_scATAC T_ALL_PATDFE_scATAC T_ALL_PAUKIZ_scATAC T_ALL_PAUUZY_scATAC T_ALL_PAVSEI_scATAC T_ALL_PATENL_scATAC T_ALL_PAUMAV_scATAC T_ALL_PAUXIY_scATAC T_ALL_PAVSRU_scATAC T_ALL_PAWIIR_scATAC T_ALL_PAUMXB_scATAC T_ALL_PAUYJE_scATAC T_ALL_PAVTCV_scATAC T_ALL_PAWRJP_scATAC T_ALL_PAREAT_scATAC T_ALL_PATEVG_scATAC T_ALL_PAUNDK_scATAC T_ALL_PAVFKN_scATAC T_ALL_PAVTXP_scATAC T_ALL_PAWRLZ_scATAC T_ALL_PARWPU_scATAC T_ALL_PATIPB_scATAC T_ALL_PAUNZE_scATAC T_ALL_PAVINC_scATAC T_ALL_PAVUFK_scATAC T_ALL_PAXHIJ_scATAC T_ALL_PASIGC_scATAC T_ALL_PATMYZ_scATAC T_ALL_PAUPTX_scATAC T_ALL_PAVJCH_scATAC T_ALL_PAVVVF_scATAC T_ALL_PASKMG_scATAC T_ALL_PATPKZ_scATAC T_ALL_PAUPVR_scATAC T_ALL_PAVLII_scATAC T_ALL_PAVVVK_scATAC T_ALL_PASWWT_scATAC T_ALL_PATTDP_scATAC T_ALL_PAURIX_scATAC T_ALL_PAVLKA_scATAC T_ALL_PAVYSA_scATAC T_ALL_PASZKM_scATAC T_ALL_PAUFAM_scATAC T_ALL_PAURXZ_scATAC T_ALL_PAVLUN_scATAC T_ALL_PAVYVY_scATAC T_ALL_7767-143_scATAC T_ALL_7767-150_scATAC T_ALL_7767-1832_scATAC T_ALL_7767-209_scATAC T_ALL_7767-2101_scATAC T_ALL_7767-2142_scATAC T_ALL_7767-2144_scATAC T_ALL_7767-2167_scATAC T_ALL_7767-2976_scATAC T_ALL_7767-31_scATAC)

sampleIDS=(T_ALL_7767-2101_scATAC T_ALL_7767-2144_scATAC)

clength=${#sampleIDS[@]}

job0='job_Slurm'
echo "#!/bin/bash" > $job0
echo "#SBATCH -c 4"  >> $job0
echo "#SBATCH --mem-per-cpu=24G"  >> $job0
echo "#SBATCH -t 99:99:99"  >> $job0
echo "module load singularity"  >> $job0

for (( i=0; i<${clength}; i++ ));
do
    jobk=${job0}_reconstructmatrix_${sampleIDS[$i]} 
    cp $job0 $jobk
    echo "singularity exec --bind /mnt -H /mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/scATACpro/${sampleIDS[$i]} --cleanenv /mnt/isilon/tan_lab/yuw1/run_scATAC-pro/NB/scatac-pro_1.3.1.sif scATAC-pro -s reConstMtx -i /mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/scATACpro/mergedPeaks_HD_Reference/peaks/merged_peaks.bed,/mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/scATACpro/${sampleIDS[$i]}/output/summary/${sampleIDS[$i]}.fragments.txt,/mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/scATACpro/${sampleIDS[$i]}/output/filtered_matrix/MACS2/FILTER/barcodes.txt,/mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/scATACpro/${sampleIDS[$i]}/output/filtered_matrix/MACS2/FILTER/reConstruct_matrix_forProjection -c ${sampleIDS[$i]}_configure_user.txt">> $jobk 
    mv $jobk ${sampleIDS[$i]}/
    cd ${sampleIDS[$i]}/
    sbatch $jobk
    cd ../
done

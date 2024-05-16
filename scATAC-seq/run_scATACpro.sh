#!/bin/bash

#sampleIDS=($(ls /mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/RawDataOrganize)) ## for all ETP-ALL samples
#sampleIDS=(T_ALL_PAVSEI_scATAC T_ALL_PAVYVY_scATAC T_ALL_7767-209_scATAC T_ALL_PAWGWD_scATAC T_ALL_PAVTCV_scATAC T_ALL_PAWIIR_PDX_scATAC T_ALL_PAVTXP_scATAC T_ALL_PAWIIR_scATAC T_ALL_PAVUFK_scATAC T_ALL_PAWRJP_scATAC T_ALL_PAVLKA_scATAC T_ALL_PAVVVF_scATAC T_ALL_PAWRLZ_scATAC T_ALL_PAVVVK_scATAC T_ALL_PAXHIJ_scATAC T_ALL_PAVYSA_scATAC T_ALL_Thymus_Age5_scATAC)
sampleIDS=(T_ALL_PAVVVF_scATAC T_ALL_PAWRLZ_scATAC)

clength=${#sampleIDS[@]}

job0='job'
echo "#!/bin/bash" > $job0
echo "#$ -cwd" >> $job0
echo "#$ -j y" >> $job0
echo "#$ -pe smp 4"  >> $job0
echo "#$ -l h_vmem=64G"  >> $job0
echo "#$ -l mem_free=64G"  >> $job0
echo "module load singularity"  >> $job0

for (( i=0; i<${clength}; i++ ));
do
    jobk=${job0}_ETPALL_scATACpro_${sampleIDS[$i]} 
    cp $job0 $jobk
    cp configure_user.txt ${sampleIDS[$i]}_configure_user.txt
    sed -i "s/control/${sampleIDS[$i]}/g" ${sampleIDS[$i]}_configure_user.txt
    mkdir ${sampleIDS[$i]}
    mv ${sampleIDS[$i]}_configure_user.txt ${sampleIDS[$i]}/
    echo "singularity exec --bind /mnt -H /mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/scATACpro/${sampleIDS[$i]} --cleanenv /mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/scATACpro/scatac-pro_1.2.1.sif scATAC-pro -s process -i /mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/RawDataOrganize/${sampleIDS[$i]} -c ${sampleIDS[$i]}_configure_user.txt -o /mnt/isilon/tan_lab/chenc6/ETP_ALL/scATAC-Seq/scATACpro/${sampleIDS[$i]}/output">> $jobk 
    mv $jobk ${sampleIDS[$i]}/
    cd ${sampleIDS[$i]}/
    qsub $jobk
    cd ../
done

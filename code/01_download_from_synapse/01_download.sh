#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N download
#$ -o logs/download.txt
#$ -e logs/download.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

#   Performed from the transfer node, interactively, for maximal speed:
#       jhpce-transfer01.jhsph.edu

cd /dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
log_path=processed-data/synapse_upload/03-download/07-download_v5.log
out_dir=raw-data/psychENCODE/version5
credential_path=/users/neagles/synapse_credentials.yaml

mkdir -p $out_dir

{
module load synapse/2.6.0

synapse \
    -p $(grep "^token:" $credential_path | cut -d " " -f 2 | tr -d '"') \
    get \
        -q "SELECT * FROM syn49290566.3" \
        --downloadLocation $out_dir

echo "Done at $(date)."
} > $log_path 2>&1

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/

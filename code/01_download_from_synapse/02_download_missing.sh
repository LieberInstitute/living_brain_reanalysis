#!/bin/bash
#$ -cwd
#$ -l rnet,mem_free=10G,h_vmem=10G,h_fsize=300G
#$ -pe local 4
#$ -N download_missing
#$ -o logs/download_missing.txt
#$ -e logs/download_missing.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the Synapse module
module load synapse/2.7.2

## List current modules for reproducibility
module list

## RNAseq FASTQ files
## Note: this includes LBP_individual_metadata.csv and
## LBP_biospeciment_metadata.csv which are duplicated with the metadata files.
synapse \
    get \
        -q "SELECT * FROM syn11346063.39 WHERE ( \"study\" HAS ( 'LBP' ) ) AND ( \"specimenID\" LIKE '%LBPSEMA4BRAIN639' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN730' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN176' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN285' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN469' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN043' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN141' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN085' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN607' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN401' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN348' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN478' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN402' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN778' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN623' )" \
        --downloadLocation ../../raw-data/psychENCODE/v2023-07-31/RNAseq

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/

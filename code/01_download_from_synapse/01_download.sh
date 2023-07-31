#!/bin/bash
#$ -cwd
#$ -l rnet,mem_free=10G,h_vmem=10G,h_fsize=300G
#$ -pe local 4
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

## Load the Synapse module
module load synapse/2.7.2

## List current modules for reproducibility
module list

## Study metadata files
synapse \
    get \
        -q "SELECT * FROM syn11346063.39 WHERE ( \`resourceType\` = 'metadata' ) AND ( ( \"study\" HAS ( 'LBP' ) ) )" \
        --downloadLocation ../../raw-data/psychENCODE/v2023-07-14/metadata

## RNAseq FASTQ files
## Note: this includes LBP_individual_metadata.csv and
## LBP_biospeciment_metadata.csv which are duplicated with the metadata files.
synapse \
    get \
        -q "SELECT * FROM syn11346063.39 WHERE ( \"study\" HAS ( 'LBP' ) )" \
        --downloadLocation ../../raw-data/psychENCODE/v2023-07-14/RNAseq

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/

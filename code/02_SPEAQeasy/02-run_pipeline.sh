#!/bin/bash
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=800G
#$ -o ../../processed-data/03_SPEAQeasy_Y2/SPEAQeasy_output.log
#$ -e ../../processed-data/03_SPEAQeasy_Y2/SPEAQeasy_output.log
#$ -cwd
#$ -N run_pipeline

#  Get absolute path to 'joelPTSD_year2' repo
base_dir=$(git rev-parse --show-toplevel)

module load nextflow/20.01.0
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow $base_dir/code/SPEAQeasy_Y2/main.nf \
    --sample "paired" \
    --reference "hg38" \
    --strand "reverse" \
    --strand_mode "declare" \
    --annotation "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/SPEAQeasy/Annotation" \
    -with-report "${base_dir}/processed-data/03_SPEAQeasy_Y2/execution_reports/02-run_pipeline.html" \
    -w "${base_dir}/processed-data/03_SPEAQeasy_Y2/work" \
    --input "${base_dir}/processed-data/03_SPEAQeasy_Y2" \
    --output "${base_dir}/processed-data/03_SPEAQeasy_Y2/pipeline_output" \
    -profile jhpce \
    -resume

#   Log successful runs on non-test data in a central location
bash $base_dir/code/SPEAQeasy_Y2/scripts/track_runs.sh $base_dir/processed-data/03_SPEAQeasy_Y2/SPEAQeasy_output.log

#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging).
echo "Generating per-sample logs for debugging..."
bash $base_dir/code/SPEAQeasy_Y2/scripts/generate_logs.sh $base_dir/processed-data/03_SPEAQeasy_Y2/SPEAQeasy_output.log

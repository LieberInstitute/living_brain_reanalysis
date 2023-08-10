#!/bin/bash
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=800G
#$ -o ../../processed-data/02_SPEAQeasy/SPEAQeasy_output.log
#$ -e ../../processed-data/02_SPEAQeasy/SPEAQeasy_output.log
#$ -cwd
#$ -N run_pipeline

#  Get absolute path to the repo origin
base_dir=$(git rev-parse --show-toplevel)

module load nextflow/22.10.7
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow $base_dir/code/02_SPEAQeasy/SPEAQeasy/main.nf \
    --sample "paired" \
    --reference "hg38" \
    --strand "reverse" \
    --strand_mode "declare" \
    --annotation "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/SPEAQeasy/Annotation" \
    -with-report "${base_dir}/processed-data/02_SPEAQeasy/execution_reports/02-run_pipeline.html" \
    -w "${base_dir}/processed-data/02_SPEAQeasy/work" \
    --input "${base_dir}/processed-data/02_SPEAQeasy/" \
    --output "${base_dir}/processed-data/02_SPEAQeasy/pipeline_output" \
    --experiment "living_brain_reanalysis" \
    -profile jhpce \
    -resume

#   Log successful runs on non-test data in a central location
bash $base_dir/code/02_SPEAQeasy/SPEAQeasy/scripts/track_runs.sh $base_dir/processed-data/02_SPEAQeasy/SPEAQeasy_output.log

#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging).
echo "Generating per-sample logs for debugging..."
bash $base_dir/code/02_SPEAQeasy/SPEAQeasy/scripts/generate_logs.sh $base_dir/processed-data/02_SPEAQeasy/SPEAQeasy_output.log

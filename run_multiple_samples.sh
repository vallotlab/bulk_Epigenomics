#!/bin/bash
sample_sheet=$1
output_dir=$2
params=$3

# Constants paths

script_dir=/data/users/gjouault/GitLab/bulk_Epigenomics/
image=/data/users/gjouault/Singularity/bulk_Epigenomics/bulkEpigenomics.sif
bind_directory=/data/
cores=20

cd $script_dir

# Copy the pipeline to the output directory & delete logs if already existing
mkdir -p $output_dir
cp -rf Scripts $output_dir/
cp -rf annotations $output_dir/
cp -f run_bulk_Epigenomics.sh $output_dir/
cp -f Snakefile_bulk_Epigenomics.py $output_dir/
cp -f species_design_configs.tsv $output_dir/
cp -f CONFIG_TEMPLATE.yaml $output_dir/
cp -f $sample_sheet $output_dir/
rm -rf $output_dir/logs/*.log

# Create sample configuration file
python3 Scripts/sample2json.py -i $sample_sheet -o $output_dir/ -c $output_dir/species_design_configs.tsv 


# Launch the pipeline using Singularity Image

echo "singularity  exec --bind ${bind_directory}:${bind_directory} --bind ${output_dir}:/mnt/ $image /bin/bash -c \"/mnt/run_bulk_Epigenomics.sh ${cores} ${params}\" " | qsub -l nodes=1:ppn=$cores,mem=70gb -N $(basename ${output_dir})


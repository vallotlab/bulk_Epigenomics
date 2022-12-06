## Bulk Epigenomics Pipeline
## Copyleft 2018 Institut Curie
## Author(s): Pacôme Prompsy
## Contact:  pacome.prompsy@curie.fr
## This software is distributed without any guarantee under the terms of the CECILL License
## See the LICENCE file for details

## Given a input string
## Determine if the input string contains two fastqs or one
## Pass the input to bowtie given options

## Logs
echo -e "Running get_fastq_func ..."
echo

inputs=("${snakemake_input[0]}")
bowtie2_option="${snakemake_params[bowtie2_option]}"
bowtie2_index="${snakemake_params[bowtie2_index]}"
output="${snakemake_output[0]}"
threads="${snakemake[threads]}"

echo "Aligning sample ${snakemake_wildcards[sample]} with bowtie2" 2> "${snakemake_log[0]}"

# Determine if 1 or 2 fastq files
if [[ ${#inputs[@]} -eq 1 ]]; then
    cmd_in="-U ${inputs[0]}"
    fi
    
    elif [[ ${#inputs[@]} -eq 2 ]]; then
        cmd_in="-1 ${inputs[0]} -2 ${inputs[1]}"
        fi
        
        # Map the reads with bowtie 
        echo "bowtie2 ${bowtie2_option} -x ${bowtie2_index} ${cmd_in} -p ${threads}  | samtools view -F4 -Sb > ${output}" 2> "${snakemake_log[0]}"
        
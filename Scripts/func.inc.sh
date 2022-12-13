## Bulk Epigenomics Pipeline
## Copyleft 2018 Institut Curie
## Author(s): Pacôme Prompsy
## Contact:  pacome.prompsy@curie.fr
## This software is distributed without any guarantee under the terms of the CECILL License
## See the LICENCE file for details

## Given a input string
## Determine if the input string contains two fastqs or one
## Pass the input to bowtie given options
bowtie2_func()
{
   ## Logs
   echo -e "Running bowtie2_func ..."
   echo
   
   inputs=($1)
   bowtie2_options=$2
   bowtie2_index=$3
   output=$4
   threads=$5
   
   echo $inputs
   echo $bowtie2_options
   echo $bowtie2_index
   echo $output
   echo $threads
   echo
   
   # Determine if 1 or 2 fastq files
   if [[ ${#inputs[@]} -eq 1 ]]; then
           cmd_in="-U ${inputs[0]}"
   fi
   if [[ ${#inputs[@]} -eq 2 ]]; then
           cmd_in="-1 ${inputs[0]} -2 ${inputs[1]}"
   fi
  
  # Map the reads with bowtie 
  bowtie2 ${bowtie2_options} -x ${bowtie2_index} ${cmd_in} -p ${threads}  | samtools view -F4 -Sb > ${output}
   
}

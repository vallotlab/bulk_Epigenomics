## Bulk Epigenomics sample2json.py
## Copyleft 2018 Institut Curie
## Author(s): Pacôme Prompsy
## Contact:  pacome.prompsy@curie.fr
## This software is distributed without any guarantee under the terms of the CECILL License
## See the LICENCE file for details

# Given a Sample Sheet and a Config file, create a JSON file that is readable
# by the Snakemake pipeline to give each sample the correct 'parameters':
# fastq paths, genome / pdx, broad or sharp mark, has control or not...

# Importing libraries
import json
import os
import csv
import re
from os.path import join
import argparse
from collections import defaultdict
import subprocess
import argparse
import pandas as pd

# Customizable constants
broad_marks = ['h3k36me3', "h3k9me3", "h3k27me3"]

# Parsing argument
parser = argparse.ArgumentParser(description='Translate Sample sheet into json file.')

parser.add_argument('-i', '--sampleSheet', type=str, required=True,
                   help='The sample sheet to translate')
parser.add_argument('-o', '--outputDirectory', type=str,required=True,
                   help='The output dir where to create sample.json')
parser.add_argument('-c', '--config', type=str, required=True,
                   help='The configuration template in .tsv containing parameters based on genome and mark.')
parser.add_argument('-x', '--chromatinIndexing', required=False, default=False,
                   action='store_const', const=True, help='If chromatinIndexing dataset, sample sheet contains index in 5th column.')

args = parser.parse_args()
print(args.sampleSheet + " and " +  args.outputDirectory)

if args.chromatinIndexing:
	print("Chromatin Indexing dataset : " + str(args.chromatinIndexing))

# Reading the species_design config csv file containing the setup for various
# 'design'
config_df = pd.read_csv(args.config, sep = "\t", index_col=False)
sample_sheet_df = pd.read_csv(args.sampleSheet, sep = ",", header = None)
print(sample_sheet_df)
# Create list of parameters for each sample
fastq_read1 = {}
fastq_read2 = {}
marks_type = {}
genomes = {}
is_control = {}
controls = {}
bowtie2_indexes = {}
second_species_bowtie2_indexes = {}
chromosome_files = {}
peak_callers = {}
peak_calling_options = {}
peak_merging_options = {}
bamCoverage_options = {}

# For each sample, fill the parameters according to the name and the config file
with open(args.sampleSheet, "r") as f:
	reader = csv.reader(f, delimiter = ',')
	for row in reader:

		sample_name = row[2].strip()
		dataset = row[0].strip()
		genome = row[3].strip()
		# Retrieve the mark from the last part of the name
		# has to bead XXXXXXXXXX_H3K27me3 or XXXXXX_EZH2 or XXXXXXX_input
		mark = re.sub(".*_", "", sample_name).lower()
		
		
		if dataset=="custom_SE" :
		  input_dir = row[1].strip()
		  fastq_file_path_R1 = subprocess.Popen(["find -L " + input_dir + " -name '%s'.R1.fastq.gz " % sample_name], stdout=subprocess.PIPE, shell=True)
		  (out, err) = fastq_file_path_R1.communicate()
		  fastq_file_path_R1 = out.decode('utf-8').strip()
		  fastq_file_path_R2 = ""
      
		elif dataset=="custom_PE" :
		  input_dir = row[1].strip()
		  fastq_file_path_R1 = subprocess.Popen(["find -L " + input_dir + " -name '%s'.R1.fastq.gz " % sample_name], stdout=subprocess.PIPE, shell=True)
		  fastq_file_path_R2 = subprocess.Popen(["find -L " + input_dir + " -name '%s'.R2.fastq.gz " % sample_name], stdout=subprocess.PIPE, shell=True)
		  (out, err) = fastq_file_path_R1.communicate()
		  fastq_file_path_R1 = out.decode('utf-8').strip()
		  (out, err) = fastq_file_path_R2.communicate()
		  fastq_file_path_R2 = out.decode('utf-8').strip()
      
		elif args.chromatinIndexing:
			fastq_name = row[1].strip()
			index = row[4].strip()
			fastq_file_path_R1= args.outputDirectory + "/FASTQ/" + sample_name + "." + index + ".R1.fastq"
			fastq_file_path_R2= args.outputDirectory + "/FASTQ/" + sample_name + "." + index + ".R2.fastq"
		else:
		  fastq_name = row[1].strip()
		  fastq_file_path_R1=subprocess.Popen(["find -L /data/kdi_prod/dataset_all/" + dataset + "/export/user/ -name '%s'.R1.fastq.gz ! -path '*after_trimming*' " % fastq_name], stdout=subprocess.PIPE, shell=True)
		  fastq_file_path_R2=subprocess.Popen(["find -L /data/kdi_prod/dataset_all/" + dataset + "/export/user/ -name '%s'.R2.fastq.gz ! -path '*after_trimming*' " % fastq_name], stdout=subprocess.PIPE, shell=True)
		  (out, err) = fastq_file_path_R1.communicate()
		  fastq_file_path_R1=out.decode('utf-8').strip()
		  (out, err) = fastq_file_path_R2.communicate()
		  fastq_file_path_R2=out.decode('utf-8').strip()

		if args.chromatinIndexing:
                        sample_name=sample_name + "." + index
                        if sample_name in fastq_read1 :
                                fastq_read1[sample_name].append(fastq_file_path_R1)
                                fastq_read2[sample_name].append(fastq_file_path_R2)
                        else:
                                fastq_read1[sample_name] = [fastq_file_path_R1]
                                fastq_read2[sample_name] = [fastq_file_path_R2]
		else:
			try:
                                fastq_read1[sample_name] = fastq_file_path_R1
                                fastq_read2[sample_name] = fastq_file_path_R2
			except KeyError:
				fastq_read1[sample_name] = fastq_file_path_R1
				
		genomes[sample_name] = genome
		
		if (mark in broad_marks):
		    mark_type = "broad"
		else:
		     mark_type = "sharp"
		     
		marks_type[sample_name] = mark_type

		# Check if the sample is a control
		sample_common_name = re.sub("_" + mark, "", sample_name, flags=re.IGNORECASE)
		if (mark == "input"):
		    is_control[sample_name] = True
		else:
		    is_control[sample_name] = False
		    if(sample_common_name + "_input" in str(sample_sheet_df[2])):
		        controls[sample_name] = sample_common_name + "_input"
		        has_control = True
		    else:
		        controls[sample_name] = False
		        has_control = False
		
		design_index = (config_df.genome_assembly == genome) & (config_df.mark_type == mark_type) & (config_df.has_control == has_control)
		bowtie2_indexes[sample_name] = str(config_df.bowtie2_index[design_index].iloc[0])
		second_species_bowtie2_indexes[sample_name] = str(config_df.second_species_bowtie2_index[design_index].iloc[0])
		chromosome_files[sample_name] = str(config_df.chromosome_file[design_index].iloc[0])
		peak_callers[sample_name] = str(config_df.peak_caller[design_index].iloc[0])
		peak_calling_options[sample_name] = str(config_df.peak_calling_option[design_index].iloc[0])
		peak_merging_options[sample_name] = str(config_df.peak_merging_option[design_index].iloc[0])
		bamCoverage_options[sample_name] = str(config_df.bamCoverage_option[design_index].iloc[0])



				
JSON = {"fastq.r1": fastq_read1, "fastq.r2": fastq_read2, "mark_type": marks_type, "genome": genomes,
"bowtie2_index": bowtie2_indexes, "second_species_bowtie2_index": second_species_bowtie2_indexes,
"chromosome_file": chromosome_files, "peak_caller": peak_callers, 
"peak_calling_option": peak_calling_options, "peak_merging_option": peak_merging_options,
"bamCoverage_option": bamCoverage_options, "control": controls, "is_control": is_control}
print(JSON)
js = json.dumps(JSON, indent = 4, sort_keys=True)
open( args.outputDirectory + '/config_samples.json', 'w').writelines(js)

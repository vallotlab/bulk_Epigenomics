# Imports
shell.executable("/bin/bash")
import pandas as pd
import numpy as np
import csv
import os

# Global Config and per sample config reading
configfile: "CONFIG_TEMPLATE.yaml"

config_sample = json.load(open("config_samples.json"))

NON_PDX_SAMPLES = []
PDX_SAMPLES = []
MARKS_TYPE = []
SHARP_MARK_SAMPLES_HG38 = []
BROAD_MARK_SAMPLES_HG38 = []
SHARP_MARK_SAMPLES_MM10 = []
BROAD_MARK_SAMPLES_MM10 = []
SAMPLES_HG38 = []
SAMPLES_MM10 = []
GENOMES = []

# For each PDX sample, need to add as output the mouse bam file
for samp in sorted(config_sample["genome"].keys()):
  mark_type = config_sample["mark_type"][samp]
  genome = config_sample["genome"][samp]
  if(genome == "pdx"):
    genome = "hg38"
  is_control = config_sample["is_control"][samp]
  is_pdx = config_sample["second_species_bowtie2_index"][samp]

  MARKS_TYPE = MARKS_TYPE + [mark_type]
  GENOMES = GENOMES + [genome]
  
  if(is_control == False):
    if (is_pdx != "False"):
      PDX_SAMPLES = PDX_SAMPLES + [samp]
    else :
      NON_PDX_SAMPLES = NON_PDX_SAMPLES + [samp]
    if (genome == "mm10"):
        SAMPLES_MM10 = SAMPLES_MM10 + [samp]
        if (mark_type == "sharp"):
            SHARP_MARK_SAMPLES_MM10 = SHARP_MARK_SAMPLES_MM10 + [samp]
        else:
            BROAD_MARK_SAMPLES_MM10 = BROAD_MARK_SAMPLES_MM10 + [samp]
    else :
        SAMPLES_HG38 = SAMPLES_HG38 + [samp]
        if (mark_type == "sharp"):
            SHARP_MARK_SAMPLES_HG38 = SHARP_MARK_SAMPLES_HG38 + [samp]
        else:
            BROAD_MARK_SAMPLES_HG38 = BROAD_MARK_SAMPLES_HG38 + [samp]

MARKS_TYPE = np.unique(np.array(MARKS_TYPE))
GENOMES = np.unique(np.array(GENOMES))
COMBINATIONS = []
if (len(BROAD_MARK_SAMPLES_HG38) > 0):
  COMBINATIONS = COMBINATIONS + ["04_csv/consensus_peaks_matrix_broad_hg38.csv"]
if (len(SHARP_MARK_SAMPLES_HG38) > 0):
  COMBINATIONS = COMBINATIONS + ["04_csv/consensus_peaks_matrix_sharp_hg38.csv"]
if (len(BROAD_MARK_SAMPLES_MM10) > 0):
  COMBINATIONS = COMBINATIONS + ["04_csv/consensus_peaks_matrix_broad_mm10.csv"]
if (len(SHARP_MARK_SAMPLES_MM10) > 0):
  COMBINATIONS = COMBINATIONS + ["04_csv/consensus_peaks_matrix_sharp_mm10.csv"]

print("NON_PDX_SAMPLES")
print(NON_PDX_SAMPLES)
print()
print("PDX_SAMPLES")
print(PDX_SAMPLES)
print()
print("SHARP_MARK_SAMPLES_HG38")
print(SHARP_MARK_SAMPLES_HG38)
print()
print("BROAD_MARK_SAMPLES_HG38")
print(BROAD_MARK_SAMPLES_HG38)
print()
print("SHARP_MARK_SAMPLES_MM10")
print(SHARP_MARK_SAMPLES_MM10)
print()
print("BROAD_MARK_SAMPLES_MM10")
print(BROAD_MARK_SAMPLES_MM10)
print()
print("SAMPLES_MM10")
print(SAMPLES_MM10)
print()
print("SAMPLES_HG38")
print(SAMPLES_HG38)
print()
print("MARKS_TYPE")
print(MARKS_TYPE)
print()
print("GENOMES")
print(GENOMES)
print()
print("COMBINATIONS")
print(COMBINATIONS)
print()

# Getters for sample options
def get_input_fastq_r1(wildcards):
    r1 = config_sample["fastq.r1"][wildcards.sample]
    return(r1)
def get_input_fastq_r2(wildcards):
    r2 = config_sample["fastq.r2"][wildcards.sample]
    return(r2)
def get_bowtie2_index(wildcards):
    ret = config_sample["bowtie2_index"][wildcards.sample]
    return(ret)
def get_second_species_bowtie2_index(wildcards):
    ret = config_sample["second_species_bowtie2_index"][wildcards.sample]
    return(ret)
def get_chromosome_bed(wildcards):
    ret = config_sample["chromosome_file"][wildcards.sample]
    return(ret)
def get_peak_caller(wildcards):
    ret = config_sample["peak_caller"][wildcards.sample]
    return(ret)
def get_control(wildcards):
    ret = config_sample["control"][wildcards.sample]
    return(ret)
def get_peak_calling_option(wildcards):
    ret = config_sample["peak_calling_option"][wildcards.sample]
    return(ret)
def get_peak_merging_option(wildcards):
    ret = config_sample["peak_merging_option"][wildcards.sample]
    return(ret)
def get_bamCoverage_option(wildcards):
    ret = config_sample["bamCoverage_option"][wildcards.sample]
    return(ret)

  
def get_input_dedup_bam(wildcards):
    is_pdx = get_second_species_bowtie2_index(wildcards)
    if(is_pdx == "False"):
      ret = "01_mapping/" + wildcards.sample + ".sorted.bam"
    else:
      ret = "01_mapping/" + wildcards.sample + ".primary.only.bam"
    return(ret)


def get_input_peak_calling_IP(wildcards):
    ret = "02_dedup_bam/" + wildcards.sample + ".bam"
    print("get_input_peak_calling_IP")
    print(ret)
    return(ret)
  
def get_input_peak_calling_input(wildcards):
    has_control = get_control(wildcards)
    print("get_input_peak_calling_input")
    print(has_control)
    if(str(has_control) != "False"):
      ret = "02_dedup_bam/" + has_control + ".bam"
    else:
      ret = "02_dedup_bam/" + wildcards.sample + ".bam"
    print(ret)
    return(ret)

def get_input_peak_consensus(wildcards):
    mark_type = wildcards.mark_type
    genome = wildcards.genome
    if(genome == "mm10"):
        if(mark_type == "sharp"):
            l = SHARP_MARK_SAMPLES_MM10
        else:
            l = BROAD_MARK_SAMPLES_MM10
    else:
        if(mark_type == "sharp"):
            l = SHARP_MARK_SAMPLES_HG38
        else:
            l = BROAD_MARK_SAMPLES_HG38  
    ret = []
    for samp in l:
        ret = ret + ["03_peaks/" + samp + ".peaks.bed"]
    
    return(ret)

def get_input_bam_counting_peaks(wildcards):
    mark_type = wildcards.mark_type
    genome = wildcards.genome
    if(genome == "mm10"):
        if(mark_type == "sharp"):
            l = SHARP_MARK_SAMPLES_MM10
        else:
            l = BROAD_MARK_SAMPLES_MM10
    else:
        if(mark_type == "sharp"):
            l = SHARP_MARK_SAMPLES_HG38
        else:
            l = BROAD_MARK_SAMPLES_HG38  
    ret = []
    for samp in l:
        ret = ret + ["02_dedup_bam/" + samp + ".bam"]
    
    return(ret)

def get_input_bam_counting_bins(wildcards):
    genome = wildcards.genome
    if(genome == "mm10"):
            l = SAMPLES_MM10
    else:
            l = SAMPLES_HG38
    ret = []
    for samp in l:
        ret = ret + ["02_dedup_bam/" + samp + ".bam"]
    
    return(ret)
  
# Snakemake Rules
rule all:
    input:
        COMBINATIONS,
        expand("04_csv/matrix_{bins}_{genome}.csv", bins = ["10k", "50k", "100k"], genome = GENOMES),
        expand("05_coverage/{sample}.bw", sample = NON_PDX_SAMPLES + PDX_SAMPLES)
        
rule bowtie2_map:
    input:
        get_input_fastq_r1,
        get_input_fastq_r2
    output:
        temp("01_mapping/{sample}.bam")
    log:
      "logs/bowtie2_map_{sample}.log"
    threads:
        5
    params:
        mem="30g",
        bowtie2_option = config['bowtie2_option'],
        bowtie2_index = get_bowtie2_index
    shell:
        """
        source ./Scripts/func.inc.sh 2> {log}
        bowtie2_func "{input}" "{params.bowtie2_option}" "{params.bowtie2_index}" "{output}" "{threads}" 2>> {log}
        """

rule bowtie2_map_secondary:
    input:
        get_input_fastq_r1,
        get_input_fastq_r2
    output:
        temp("01_mapping/{sample}.secondary.bam")
    log:
      "logs/bowtie2_map_secondary_{sample}.log"
    threads:
        5
    params:
        mem="30g",
        bowtie2_option = config['bowtie2_option'],
        bowtie2_index = get_second_species_bowtie2_index,
    shell:
        """
        source ./Scripts/func.inc.sh 2> {log}
        bowtie2_func "{input}" "{params.bowtie2_option}" "{params.bowtie2_index}" "{output}" "{threads}" 2>> {log}
        """
        
rule sort_bam:
    input:
        "01_mapping/{sample}.bam"
    output:
        "01_mapping/{sample}.sorted.bam"
    log:
      "logs/sort_bam_{sample}.log"
    threads:
        5
    params:
        mem="20g",
        chromosome_bed = get_chromosome_bed,
        is_pdx = get_second_species_bowtie2_index
    shell:
        """
        mkdir -p tmp/
        if [[ {params.is_pdx} == "False" ]]; then 
          sambamba view -f bam -L {params.chromosome_bed} {input} | sambamba sort --tmpdir tmp/ -t {threads} /dev/stdin -o {output} 2> {log}
        else
          sambamba view -f bam -L {params.chromosome_bed} {input} | sambamba sort --tmpdir tmp/ -N -t {threads} /dev/stdin -o {output} 2> {log}
        fi
        """

rule sort_secondary_bam:
    input:
        "01_mapping/{sample}.secondary.bam"
    output:
        "01_mapping/{sample}.secondary.nsorted.bam"
    log:
      "logs/sort_bam_secondary_{sample}.log"
    threads:
        5
    params:
        mem="20g"
    shell:
        """
        mkdir -p tmp/
        sambamba sort --tmpdir tmp/ -N -t {threads} {input} -o {output} 2>> {log}
        """

rule filter_out_secondary:
    input:
        primary="01_mapping/{sample}.sorted.bam",
        secondary="01_mapping/{sample}.secondary.nsorted.bam"
    output:
        "01_mapping/{sample}.primary.only.bam"
    log:
        "logs/{sample}_filtering_mouse"
    threads:5
    params:
        mem="32g",
        time="02:00:00"
    shell:
        """
        mkdir -p tmp/
        bamcmp -1 {input.primary} -2 {input.secondary} -a {output}.unsorted -A {output}.temp -t {threads} -n -s mapq 2> {log}
        sambamba sort --tmpdir tmp/ -t {threads} {output}.unsorted -o {output} 2>> {log}
        rm {output}.temp
        """
        
rule dedup_bam:
    input:
        get_input_dedup_bam
    output:
        "02_dedup_bam/{sample}.bam"
    log:
        "logs/{sample}_dedup_bam"
    params:
        dup="logs/{sample}.dupstats",
        mem = "32g",
        time = "01:00:00",
        mark_duplicate_option = config['mark_duplicate_option']
    threads:5
    shell:
        """
        picard MarkDuplicates INPUT={input} OUTPUT={output} METRICS_FILE={params.dup} {params.mark_duplicate_option}
        sambamba index -t {threads} {output} 2> {log}
        """


rule peak_calling:
    input:
        ip=get_input_peak_calling_IP,
        control=get_input_peak_calling_input
    output:
        "03_peaks/{sample}.peaks.bed"
    log: "logs/{sample}_peak_calling"
    threads:1
    params:
        mem="10g",
        time="00:20:00",
        peak_caller=get_peak_caller,
        peak_calling_option=get_peak_calling_option,
        peak_merging_option=get_peak_merging_option,
        control=get_control
    shell:
        """
        if [[ {params.peak_caller} == "zerone" ]]; then
          zerone {params.peak_calling_option} -0 {input.control} -1 {input.ip} | cut -f1-3 | bedtools merge -d {params.peak_merging_option} -i /dev/stdin > {output} 2> {log}
        fi
        if [[ {params.peak_caller} == "macs2" ]]; then
          macs2 callpeak {params.peak_calling_option} -t {input.ip} -n {output} 2> {log}
	  rm -f {output}_peaks.gappedPeak
	  cut -f1-3 {output}_peaks.*Peak | bedtools sort -i /dev/stdin/ | bedtools merge -d {params.peak_merging_option} -i /dev/stdin > {output} 2>> {log}
        fi
        """

rule consensus_annotation:
    input:
        get_input_peak_consensus
    output:
        "03_peaks/consensus_peaks_{mark_type}_{genome}.bed"
    log: "logs/consensus_annotation_{mark_type}_{genome}"
    threads:1
    params:
        mem="5g",
        time="00:30:00"
    shell:
        "cat {input} | bedtools sort -i /dev/stdin | mergeBed -i /dev/stdin | bedtools merge -d 1000 -i /dev/stdin > {output}"

rule run_pysam_consensus:
    input:
        samples=get_input_bam_counting_peaks,
        annot="03_peaks/consensus_peaks_{mark_type}_{genome}.bed",
    output:
        "04_csv/consensus_peaks_matrix_{mark_type}_{genome}.csv"
    log: "logs/run_pysam_peaks_{mark_type}_{genome}"
    threads:1
    params:
        mem="32g",
        time="03:00:00"
    shell:
        """
        genome=$(echo {input.annot} | sed 's/.*_\|.bed//g')
        python3.6 Scripts/pysam_normalization.py -g $genome -b "{input.samples}" -a {input.annot} -o {output}  2> {log}
        """

rule run_pysam_bins:
    input:
        samples=get_input_bam_counting_bins,
        annot="annotations/{genome}.{bins}.bed"
    output:
        "04_csv/matrix_{bins}_{genome}.csv"
    log: "logs/run_pysam_bins_{bins}_{genome}"
    threads:1
    params:
        mem="32g",
        time="03:30:00"
    shell:
        """
        genome=$(echo {input.annot} | sed 's/.*_\|.bed//g')
        python3.6 Scripts/pysam_normalization.py -g $genome -b "{input.samples}" -a {input.annot} -o {output}  2> {log}
        """

rule run_bamcoverage:
    input:
        "02_dedup_bam/{sample}.bam"
    output:
        "05_coverage/{sample}.bw"
    log: "logs/{sample}_bamCoverage"
    threads: 5
    params:
        mem="32g",
        time="02:00:00",
        bamCoverage_option=get_bamCoverage_option
    shell:
        """
        bamCoverage -b {input} -o {output} --numberOfProcessors {threads} {params.bamCoverage_option} 2> {log}
        """

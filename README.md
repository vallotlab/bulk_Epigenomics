# Bulk Epigenomics Pipeline - from FASTQ to Count Matrices
  
  
Snakemake pipelines to process FASTQ files from bulk epigenomics such as ChIP-seq, CUT&Tag, Chromatin Indexing. It can take as input Human samples (hg38), mouse samples (mm10) or PDX samples (hg38 for the tumor, mm10 for the mouse Tumor Micro Environment).  
All the tools needed are embedded in a Singularity Environment, allowing you to 
run the pipeline in a containarized environment (see [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)).

  
  
    
![](Pipeline_schematic.png?raw=true "Pipeline DAG")

    
   
*Command Line:*  
```
bash run_multiple_samples.sh ../sample_sheets/SampleSheet_CutTag_Human_Tumors.tsv $kdi/ChIP_seq/Test_bulkEpigenomics_CutTag_hT/
```

  
## Set up & Requirements
  
  

#### 1. Download pipeline

In order to set up the bulk Epigenomics pipeline, first download the github repository
to a directory of your choice: 

```
git clone git@github.com:vallotlab/bulk_Epigenomics
```
  
  
#### 2. Download Singularity & Singularity Image
  
  
Then, download the Singularity Image at *link comming soon*, containing all 
the tools needed for each step of the package. This means you do not need any
additional installation except : 

 * singularity (see [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html))
 * python3 ([Python](https://www.python.org/downloads/)) and pandas python package ([pandas](https://pandas.pydata.org/))  
  
#### 3. Build Bowtie2 index & indicate path in design file

You need to have a bowtie2 index of either Human (hg38) or Mouse (mm10) genomes (see [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)).
In the **species_design_configs.tsv**, you'll need to modify all the **bowtie2_index** and
**second_species_bowtie2_index** columns with the prefix towards the bowtie2 indexes.
  
  
#### 4. Modify the launching script with your specific paths
  
  
You finally need to modify the **run_multiple_samples.sh** script changing:  

 * script_dir=~/GitLab/bulk_Epigenomics/  -> Path towards the downloaded repository
 * image=~/Singularity/bulk_Epigenomics/bulkEpigenomics.sif -> Path towards the downloaded Image
 * bind_directory=/data/  -> Root directory of the directory where the FASTQ files are located. This directory will be mounted in the container.
 * cores=20  -> Number of cores you want to use

You are now set up and can move towards creating your sample sheet !
  
  
## Launching the pipeline
  
  
Now copy and modify the 'SampleSheet_test_PE.csv' sample sheet for paired-end data or the
'SampleSheet_test_SE.csv' sample sheet for single-end data.

*Note : the SampleSheet_template.csv is formatted for use on the Institut Curie HPC, and can be used to run the 
pipeline on output FASTQs from the KDI*

Now launch the pipeline with the following command:

```
bash run_multiple_samples.sh ../sample_sheets/SampleSheet_CutTag_Human_Tumors.tsv $kdi/ChIP_seq/Test_bulkEpigenomics_CutTag_hT/
```




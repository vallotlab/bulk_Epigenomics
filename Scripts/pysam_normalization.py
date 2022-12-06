import argparse
import os
import os.path
import sys
import time
import subprocess
import itertools
import numpy as np
import pandas as pd
import re
import pysam
import timeit
import ntpath

def pysam_parse_option():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam files", dest="bam_files", help="Path to the .bam files to count")
    parser.add_argument("-g", "--organism genome", dest="genome", help="genome (mm10 / hg38) used")
    parser.add_argument("-a", "--annotation file", dest="annot", help="Annotation bed file (windows, peaks) to use")
    parser.add_argument("-o", "--output", dest="out_file", help="CSV File to write the count matrix")
    args = parser.parse_args()
    bam_files = args.bam_files
    annot = args.annot
    genome = args.genome
    out_file = args.out_file
    bam_files = list(bam_files.split(" "))
    print(bam_files)
    return(bam_files,genome, annot,out_file)

def matrix(bam_files, genome, annot_file):
  
    df = pd.read_csv(annot_file, sep="\t",names=["Chromosome","Begin","End"])
    print(bam_files)
    print(len(bam_files))
    for i in range(0,len(bam_files)):
        file = bam_files[i]
        name = os.path.basename(file)
        name = name.replace(".bam", "")
        print(file)
        print(name)
        if os.path.isfile(file)==False: sys.exit(file + "doesn't exist")
        file_sam=pysam.AlignmentFile(file,"rb")
        
        def pysam_count(x,y,z):
            return(file_sam.count(x,y,z))
        
        start = timeit.default_timer()
        print("Counting IP reads for "+ name +" and normalizing")
        
        df[name]=list(map(pysam_count,df['Chromosome'],df['Begin'],df['End']))
        stop = timeit.default_timer()
        time = str(round(stop-start))
        print("Done in " + time + " seconds")                
    return(df)


def export(matrice, annot_file, out_file):
    basename = ntpath.basename(annot_file)
    name = os.path.splitext(basename)[0]
    matrice.to_csv(out_file, index=0)


def main():
    (bam_files, genome,annot,out_file) = pysam_parse_option()
    if os.path.isfile(annot)==False: sys.exit("annot file doesn't exist")  
    name = os.path.splitext(ntpath.basename(annot))[0]
    print(name + " process started at " + time.strftime("%c"))
    print("Creating normalized matrixes from bam files...")
    data = matrix(bam_files, genome, annot)
    print('Exporting results to csv...')
    export(data, annot, out_file)
    print("Done")

if __name__ == "__main__":
    main()


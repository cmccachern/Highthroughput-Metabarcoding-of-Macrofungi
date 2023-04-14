import os                #command-line like functions, for operating system interface (finding files)
import subprocess        #recommended way of running command line programs from within python
import numpy as np       #for math and arrays
import shutil            #used to move files around
import sys               #helps with reading and writing onto text files
from tqdm import tqdm    #a progress bar for for-loops (lets you see progress in actively running loops)
from time import time    #use to track time and measure how long chunks are taking to run
import shutil            #for moving files from directory to directory
import csv               #writing and reading csv
import pandas as pd
 
##    FOR DNA
 
from Bio import Seq                 #reading in and manipulating sequence data
from Bio import SeqIO               #for reading and writing fasta/qs
from Bio.SeqRecord import SeqRecord #creating sequence records that are objects and not just strings
from Bio.SeqUtils import GC         #for calculating GC content
 
##    FOR FIGURES
 
import matplotlib.pyplot as plt          #basic plotting tools that let you do most of what you'll need to do, set up as a shortcut
import matplotlib                        #get ALL of matplotlib for fancier tools
                                         #add specific other matplotlib imports as needed

def summarize_fastqs(files, outfile):
    """
    Define a function for getting fastq file summaries
    Warning, it does math with quality scores so it's slow.
    input a list of files (full or relative path)
    output a tab delimited file of summary stats 
    about the sequences in the files
    """
    #Initialize summary stat table, one line per fastq file
    perfile = {}
    #Looking at the whole file
    perfile['path'] = []
    perfile['file'] = []
    perfile['n sequences'] = []
    perfile['total bases'] = []
    #Looking at the lengths of the sequences
    perfile['mean sequence length'] = []
    perfile['standard deviation sequence length'] = []
    perfile['median sequence length'] = []
    perfile['minimum sequence length'] = []
    perfile['maximum sequence length'] = []
    #Looking at the pool of individual base quality scores
    perfile['mean total base quality score'] = []
    perfile['standard deviation total base quality score'] = []
    perfile['min total base quality score'] = []
    perfile['max total base quality score'] = []
    #Looking at the read average quality scores
    perfile['mean read average quality score'] = []
    perfile['standard deviation read average quality score'] = []
    perfile['min read average quality score'] = []
    perfile['max read average quality score'] = []
    
    #Loop through files to collect information
    for ffile in tqdm(files):
        lengths = [] #read lengths
        qualities = []#qual values for all bases
        mean_qualities = []#read mean qual averages
        for rec in SeqIO.parse(ffile, 'fastq'):
            lengths.append(len(str(rec.seq)))
            seqqualities = rec.letter_annotations["phred_quality"] 
            qualities.extend(seqqualities)
            mean_qualities.append(np.mean(seqqualities))
            
        #Calculate and store summary stats. 
        perfile['path'].append(ffile)
        perfile['file'].append(ffile.split("/")[-1])
        perfile['n sequences'].append(len(lengths))
        perfile['total bases'].append(np.sum(lengths))
        
        #Looking at the lengths of the sequences
        perfile['mean sequence length'].append(np.mean(lengths))
        perfile['standard deviation sequence length'].append(np.std(lengths))
        perfile['median sequence length'].append(np.median(lengths))
        perfile['minimum sequence length'].append(np.min(lengths))
        perfile['maximum sequence length'].append(np.max(lengths))
        
        #Looking at the pool of individual base quality scores
        perfile['mean total base quality score'].append(np.mean(qualities))
        perfile['standard deviation total base quality score'].append(np.std(qualities))
        perfile['min total base quality score'].append(np.min(qualities))
        perfile['max total base quality score'].append(np.max(qualities))
        
        #Looking at the read average quality scores
        perfile['mean read average quality score'].append(np.mean(mean_qualities))
        perfile['standard deviation read average quality score'].append(np.std(mean_qualities))
        perfile['min read average quality score'].append(np.min(mean_qualities))
        perfile['max read average quality score'].append(np.max(mean_qualities))
        
    df = pd.DataFrame(perfile)
    df.to_csv(outfile)

def locate_input_files(input_dir):
    #Locate raw input files (includes already split replicates)
    allrawdir = os.listdir(input_dir)
    forward_fastqs = []
    reverse_fastqs = []
    for ffile in allrawdir:
        if ffile.endswith("_R1.fastq"):
            forward_fastqs.append(os.path.join(input_dir, ffile))
        elif ffile.endswith("_R2.fastq"):
            reverse_fastqs.append(os.path.join(input_dir, ffile))
 
    #Sort to make sure pairs are in the same order in these lists
    forward_fastqs = sorted(forward_fastqs)
    reverse_fastqs = sorted(reverse_fastqs)

    return forward_fastqs, reverse_fastqs

def summarize_fastq_directory(input_dir,
                              forward_out_file="1_Raw_forward_fastqs_summary.txt",
                              reverse_out_file="1_Raw_reverse_fastqs_summary.txt"):
    
    forward_fastqs, reverse_fastqs = locate_input_files(input_dir)

    #Report back on findings (and read to make sure sensible)
    print(f"Found {len(forward_fastqs)} forward files with {len(set(forward_fastqs))} unique names, like {forward_fastqs[0]}")
    print(f"Found {len(reverse_fastqs)} reverse files with {len(set(reverse_fastqs))} unique names, like {reverse_fastqs[0]}")

    # Summarize Files: Evaluates all the reads in your files and gives summary data.
    summarize_fastqs(forward_fastqs, forward_out_file)
    summarize_fastqs(reverse_fastqs, reverse_out_file)

def run_cutadapt(input_dir, trimmed_dir, cutadapt_logfile="2_cutadapt_noprimers.txt"):
    # Run CutAdapt
    # This will remove primers and spacers from the 3' end of the sequence reads.

    forward_fastqs, reverse_fastqs = locate_input_files(input_dir)

    #Define Primers

    #Forward primer (ITS7f) - without spacer (on 5' end, distal to cuts), with R for degenerate base (cutadapt fine with that)
    forward_primer = Seq.Seq("GTGARTCATCGAATCTTTG")            #converting to sequence object for handy complementing
    forward_primer_complement = forward_primer.reverse_complement()
    
    #Reverse primer (ITS4) - without spacer (on 5' end) and without 5 Ns for barcode region (between spacer and primer) 
    reverse_primer = Seq.Seq("TCCTCCGCTTATTGATATGC")
    reverse_primer_complement = reverse_primer.reverse_complement()
    
    #Make copies of those in simple string form for feeding to cutadapt
    fprimer = str(forward_primer)
    fprimer_rc = str(forward_primer_complement)
    rprimer = str(reverse_primer)
    rprimer_rc = str(reverse_primer_complement)

    #Save stdout to a file along the way - cutadapt prints lots of reporting information here
    with open(cutadapt_logfile, "w") as stdouthandle:
        
        #Loop over pairs
        for ffile in tqdm(forward_fastqs):
            
            rfile = ffile.replace("_R1.fastq","_R2.fastq")
            foutfastq = ffile.replace("R1.fastq","R1_noprimer.fastq")
            routfastq = rfile.replace("R2.fastq","R2_noprimer.fastq")
            
            #Compose command
            #-g, -G removes from 5' end for f and r reads
            #-a, -A removes given primer from 3' end for f and r reads
            cmd = f"cutadapt -a {rprimer_rc} -A {fprimer_rc} -m 200 -o {foutfastq} -p {routfastq} {ffile} {rfile}"
            
            subprocess.call(cmd, stdin=None, stdout=stdouthandle, stderr=subprocess.STDOUT, shell=True)        

            # Move file to output directory
            shutil.move(foutfastq, trimmed_dir)
            shutil.move(routfastq, trimmed_dir)

#def create_no_primers_folder():
#    # ### No Primers Folder
#    # Make a folder to put your primer-free sequences in.
#
#    #Make a directory for files without primers & adapters, cutadapt results:
#    get_ipython().run_line_magic('mkdir', 'NoPrimers')
#
#    TOTAL = os.listdir()                  #take everything in the current directory and call it "TOTAL"
#    MAIN = f"{MAINdir}/"                 #name the full directory path to the main directory "MAIN"
#    NoPrimer = f"{TRIMMEDdir}/"          #name the destination directory "NoPrimer"
#    for file in TOTAL:                    #for-loop regarding new variable "file1" in current directory
#        if ("_noprimer.fastq" in file):   #select files with "_TPaired.fastq" in their name (these were the files, forward and reverse, that did not lose their pair in trimming/filtering)
#            src = MAIN+file               #define the source location of the file in question
#            dst = NoPrimer+file           #define the destination location of the file in question
#            shutil.move(src,dst)          #move the file from source to destination


if __name__ == "__main__":
    
    results_dir = os.path.expanduser("~/Files/dev/fungi/Highthroughput-Metabarcoding-of-Macrofungi/results")
    
    MAINdir = os.path.expanduser("~/Files/dev/fungi/Highthroughput-Metabarcoding-of-Macrofungi/data/20221213_2022A_MarCorNAMA/SplitBySample")
    trimmed_dir = f"{results_dir}/NoPrimers"
    if not os.path.exists(trimmed_dir):
        os.mkdir(trimmed_dir)


    summarize_fastq_directory(input_dir=MAINdir,
                              forward_out_file=f"{results_dir}/1_Raw_forward_fastqs_summary.txt",
                              reverse_out_file=f"{results_dir}/1_Raw_reverse_fastqs_summary.txt")

    run_cutadapt(MAINdir, trimmed_dir, cutadapt_logfile=os.path.join(results_dir, "2_cutadapt_noprimers.txt"))

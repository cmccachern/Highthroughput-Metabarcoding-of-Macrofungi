import os

from run_cutadapt import summarize_fastq_directory, run_cutadapt
from run_dada2 import run_dada2

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


    inpath = os.path.join("results_dir", "NoPrimers")
    inpath2 = os.path.join("results_dir", "NoPrimers/Filter_Trim")
 
    run_dada2(inpath, inpath2, results_dir)

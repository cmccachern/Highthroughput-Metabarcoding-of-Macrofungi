import os
import subprocess


def run_dada2(inpath, inpath2, outpath):
    #inpath = os.path.expanduser("~/Files/dev/fungi/Highthroughput-Metabarcoding-of-Macrofungi/results/NoPrimersTest")
    #inpath2 = os.path.expanduser("~/Files/dev/fungi/Highthroughput-Metabarcoding-of-Macrofungi/results/NoPrimersTest/Filter_Trim")
    #outpath = os.path.expanduser("~/Files/dev/fungi/Highthroughput-Metabarcoding-of-Macrofungi/results")
    
    
    subprocess.run(["Rscript", os.path.join(os.path.dirname(__file__), "run_dada2.r"), inpath, inpath2, outpath], check=True)


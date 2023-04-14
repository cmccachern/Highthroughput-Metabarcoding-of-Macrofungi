#---
#title: "DADA2 pipeline4 _ QR"
#author: "Jessie Berta-Thompson"
#date: "3/31/2021"
#output: html_document
#---
#
#The R steps of dada2 sequence processing pipeline. 
#
#DADA2 pipeline following tutorials here:
#https://benjjneb.github.io/dada2/tutorial.html
#https://benjjneb.github.io/dada2/ITS_workflow.html
#
#and working off dada2_pipeline.Rmd from pipeline 3 - work will be similar but with more permissive filter parameters.
#
#Dataset:
#Russula illumina preliminary data (186 specimens)
#
#Data prep:
#pipeline4.ipynb (but not using those trims, so that I can show Gary dada2 trimming and because I'm still trying to find trimmomatic parameters that fix poor quality samples; this will be a clean head-to-head with pipeline 3 on dada2 filter parameters). 
#
#Set up R environment 

library(BiocManager)
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)


#Filter & trim primer-free reads for quality

#Define variable to hold working directory with fastq files after primer removal (cutadapt pre-processing)
inpath <- path.expand("~/Files/dev/fungi/Highthroughput-Metabarcoding-of-Macrofungi/results/NoPrimersTest")
inpath2 <- path.expand("~/Files/dev/fungi/Highthroughput-Metabarcoding-of-Macrofungi/results/NoPrimersTest/Filter_Trim")

outpath <- "/Users/careymccachern/Files/dev/fungi/Highthroughput-Metabarcoding-of-Macrofungi/results"

#Locate the input files
cutFs <- sort(list.files(inpath, pattern = "_R1_noprimer.fastq", full.names = TRUE))
cutRs <- sort(list.files(inpath, pattern = "_R2_noprimer.fastq", full.names = TRUE))

print("Located input files (forward, reverse).")
print(length(cutFs))
print(length(cutRs))

head(cutFs)
head(cutRs)

# Extract sample names from file names - this will be used to label data as we work

# Define a function to get sample name from file name
get.sample.name <- function(fname){
  specimen <- strsplit(basename(fname), "_R")[[1]][1]
}

#Get sample names from file names
sample.names <- unname(sapply(cutFs, get.sample.name))

#Take a look at the first few
head(sample.names)

#For data trimming, this command will:
#- maxN: Remove all Ns (I don't think there are any)
#- maxEE: Remove reads with more than 3 expected errors based on quality scores
#- truncQ: Remove portions of reads from the first instance of quality score 2
#- minLen: Remove reads shorter than 200 bp after trimming
#- compress: if they are compressed (".gz") or not
#- rm.phix: Screen for and remove phiX sequence
#
#To follow progress, check that 2_filterAndTrim was created, 
#Define a function to make filenames for this step
make.filtered.nameR1 <- function(inputname){
    first_part <- strsplit(inputname, "_R1_", fixed = "True")[[1]][1]
    filtered_name <- paste(first_part, "_R1_filtered.fastq", sep = "")
}

make.filtered.nameR2 <- function(inputname){
    first_part <- strsplit(inputname, "_R2_", fixed = "True")[[1]][1]
    filtered_name <- paste(first_part, "_R2_filtered.fastq", sep = "")
}

#Assign file paths for where to put output filtered reads
filtFs <- file.path(inpath, "Filter_Trim", sapply(basename(cutFs), make.filtered.nameR1))
filtRs <- file.path(inpath, "Filter_Trim", sapply(basename(cutRs), make.filtered.nameR2))

#Make sure naming worked
head(filtFs)
head(filtRs)

#Run filtering & trimming function
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 200, compress = FALSE, rm.phix = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE


#Output is a summary count table (and files are created on disk)
head(out)

#To check progress of filtering, in file browser, look for creation of Filter_Trim and it's contents (you can watch new files appear).

#Save summary counts of pre/post filter reads to file

#print("What type of R object is filtering summary output in?")
#typeof(out)
#class(out)

#Save summary of filtering to file (and test some file saving features - default has weird column labeling)
#Add filename as a proper column (currently names which save poorly)
withrownamesascolumn = data.frame("Filename"=rownames(out),out)
#print(withrownamesascolumn['Filename'])

#Add samples names (same order made from file names above) as another labeling column
withsamplecolumn = data.frame("Sample"=sample.names,withrownamesascolumn)
#head(withsamplecolumn)

#Write to tab-delimited file
write.table(withsamplecolumn, file=paste(outpath, "4_filtering_summary.txt", sep=""), row.names=FALSE, col.names=TRUE, sep = "\t", quote = FALSE)

#Take a quick look at summary file. 


#Filtering step slow. If already run, can alternatively access filtered fastqs like this:
#Define variable to hold working directory with fastq files after primer removal (cutadapt pre-processing)

#Locate the input files
filtFs <- sort(list.files(inpath2, pattern = "R1_filtered.fastq", full.names = TRUE))
filtRs <- sort(list.files(inpath2, pattern = "R2_filtered.fastq", full.names = TRUE))

print("Located filtered files (forward, reverse).")
print(length(filtFs))
print(length(filtRs))

#Run read error profiling function

#Learn errors (run machine learning thing)
print("Running learn errors function on forward reads.")
errF <- learnErrors(filtFs, multithread=TRUE)

print("Running learn errors function on reverse reads.")
errR <- learnErrors(filtRs, multithread=TRUE)

#Plot error profiles (inspect)
print("Creating plots for error models.")
ferrorplot <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste(outpath, "5_error_models_forward_reads_plot.pdf", sep=""), plot = ferrorplot)
rerrorplot <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste(outpath, "5_error_models_reverse_reads_plot.pdf", sep=""), plot = rerrorplot)

#Save plots of error models to file
#ggsave("error_models_forward_reads_plot.pdf", plot = ferrorplot)
#ggsave("error_models_reverse_reads_plot.pdf", plot = rerrorplot)

#Run sample inference

#Run sample inference model   (find unique and de-noise)
print("Running sample inference model, forward reads.")
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
head(dadaFs)
#multithread may need to be FALSE for my computer, but check.

print("Running sample inference model, reverse reads.")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
head(dadaRs)

#Create sequence tables for saving these intermediate results to file
seqtabF <-makeSequenceTable(dadaFs)
seqtabR <-makeSequenceTable(dadaRs)

#Merge forward and reverse read data

# Merge forward and reverse reads     ("overlapping")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
print(mergers[[1]])

#Create a sequence table from the merged table
seqtab <- makeSequenceTable(mergers)
print(dim(seqtab))

# Inspect distribution of sequence lengths
print(table(nchar(getSequences(seqtab))))

## chimera- sometimes in PCR, the front end and back end of a read may have physically come from different biological sequences ((see ancient Greek half-n-half animal thingies))
## bimera- two pieces

# Chimera Removal
print("Filtering to remove chimeras")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="per-sample", multithread=FALSE, verbose=TRUE)
## de novo = from scratch
## future experiment: "per-sample" instead of "consensus" for Bimera 'method='

print("Dimensions of sequence table after chimera removal")
dim(seqtab.nochim)
print("Fraction of total abundance counts remaining after chimera removal (whole dataset)")
sum(seqtab.nochim)/sum(seqtab)
print("Fraction of unique sequences remaining after chimera removal (whole dataset")
dim(seqtab.nochim)[[2]]/dim(seqtab)[[2]]

#Summarize process
print(dadaFs)

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
print(track)

#Save that to file. 
track_labelled = data.frame("Sample"=rownames(track),track)

#Write to tab-delimited file
write.table(track_labelled, file=paste(outpath, "6_readcount_summary.txt", sep=""), row.names=FALSE, col.names=TRUE, sep = "\t", quote = FALSE)
 
# Extract sample names again from filtered file names - to account for lost specimen after the filter and trim step.

#Get sample names from file names (from filtF this time instead of cutF)
sample.names2 <- unname(sapply(filtFs, get.sample.name))

#Take a look at the first few
head(sample.names2)
length(sample.names2)




#Save results (sequence tables) to file

#Add better labels (output format funny)
seqtabF_labelled = data.frame("Filename" = rownames(seqtabF), seqtabF)
seqtabF_labelled  = data.frame("Specimen" = sample.names2, seqtabF_labelled)
write.table(seqtabF_labelled , file=paste(outpath, "7a_inference_table_Forward.txt", sep=""), row.names=FALSE, col.names=TRUE, sep = "\t", quote = FALSE)

seqtabR_labelled  = data.frame("Filename" = rownames(seqtabR), seqtabR)
seqtabR_labelled  = data.frame("Specimen" = sample.names2, seqtabR_labelled)
write.table(seqtabR_labelled, file=paste(outpath, "7b_inference_table_Reverse.txt", sep=""), row.names=FALSE, col.names=TRUE, sep = "\t", quote = FALSE)

seqtab_labelled = data.frame("Filename" = rownames(seqtab), seqtab)
seqtab_labelled = data.frame("Specimen" = sample.names2, seqtab_labelled)
write.table(seqtab_labelled, file=paste(outpath, "7c_inference_table_merged.txt", sep=""), row.names=FALSE, col.names=TRUE, sep = "\t", quote = FALSE)

seqtab.nochim_labelled = data.frame("Filename" = rownames(seqtab.nochim), seqtab.nochim)
seqtab.nochim_labelled = data.frame("Specimen" = sample.names2, seqtab.nochim_labelled)
write.table(seqtab.nochim_labelled, file=paste(outpath, "7d_inference_table_chimerasremoved.txt", sep=""), row.names=FALSE, col.names=TRUE, sep = "\t", quote = FALSE)

#Save summary of filtering to file (and test some file saving features - default has weird column labeling)
#Add filename as a proper column (currently names which save poorly)
withrownamesascolumn = data.frame("Filename"=rownames(out),out)
#print(withrownamesascolumn['Filename'])

#Add samples names (same order made from file names above) as another labelling column
withsamplecolumn = data.frame("Sample"=sample.names2,withrownamesascolumn)
#head(withsamplecolumn)

#Write to tab-delimited file
write.table(withsamplecolumn, file=paste(outpath, "8_filtering_summary.txt", sep=""), row.names=FALSE, col.names=TRUE, sep = "\t", quote = FALSE)

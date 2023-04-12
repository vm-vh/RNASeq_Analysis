# !/bin/bash
# This script checks the qualitiy of our fastq files and performs an alignment to the tomato cDNA transcriptome reference with Kallisto.
# To run this 'shell script' you will need to open your terminal and navigate to the directory where this script resides on your computer.
# This should be the same directory where you fastq files and reference fasta file are found.
# Change permissions on your computer so that you can run a shell script by typing: 'chmod 777 readMapping.sh' (without the quotes) at the terminal prompt 
# Then type './readMapping.sh' (without the quotes) at the prompt.  
# This will begin the process of running each line of code in the shell script.

# first use fastqc to check the quality of our fastq files:
# fastqc *.gz -t 4

# next, we want to build an index from our reference fasta file 
# I get my reference plant transcriptome files from here: https://plants.ensembl.org/info/data/ftp/index.html
kallisto index -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.fa

# now map reads to the indexed reference host transcriptome
# use as many 'threads' as your machine will allow in order to speed up the read mapping process.
# note that we're also including the '&>' at the end of each line
# this takes the information that would've been printed to our terminal, and outputs this in a log file that is saved in /data/course_data

kallisto quant -i Solanum_lycopersicum.SL3.0.cdna.all.index -o k_SRR17655898  --single -l 180 -s 20 SRR17655898.trim.fastq.gz

kallisto quant -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index -o k_SRR --single -l 180 -s 20 SRR.trim.fastq.gz

# summarize fastqc and kallisto mapping results in a single summary html using MultiQC
# multiqc -d . 

echo "Finished"

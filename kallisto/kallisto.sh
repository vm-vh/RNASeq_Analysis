# !/bin/bash
# This script checks the qualitiy of the fastq files and performs an alignment to the cDNA transcriptome reference with Kallisto.

# build an index from the reference fasta file 
# reference plant transcriptome file: https://plants.ensembl.org/Coffea_canephora/Info/Index.html
kallisto index -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.fa

# individual stems
kallisto quant -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index -o k_SRR17673033 --single -l 180 -s 20 SRR17673033.trim.fastq.gz &> kallisto.log
kallisto quant -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index -o k_SRR17673032 --single -l 180 -s 20 SRR17673032.trim.fastq.gz &>> kallisto.log
kallisto quant -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index -o k_SRR17673031 --single -l 180 -s 20 SRR17673031.trim.fastq.gz &>> kallisto.log
kallisto quant -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index -o k_SRR17673030 --single -l 180 -s 20 SRR17673030.trim.fastq.gz &>> kallisto.log
# individual leaves
kallisto quant -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index -o k_SRR17673027 --single -l 180 -s 20 SRR17673027.trim.fastq.gz &>> kallisto.log
kallisto quant -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index -o k_SRR17673026 --single -l 180 -s 20 SRR17673026.trim.fastq.gz &>> kallisto.log
kallisto quant -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index -o k_SRR17673025 --single -l 180 -s 20 SRR17673025.trim.fastq.gz &>> kallisto.log
kallisto quant -i Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.index -o k_SRR17673024 --single -l 180 -s 20 SRR17673024.trim.fastq.gz &>> kallisto.log

echo "Finished"

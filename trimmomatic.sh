# !/bin/bash

# trimmomatic
# needs to be run in the same folder as both the data and trimmomatic
java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar SE SRR17673024.fastq.gz SRR17673024.trim.fastq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:70 -phred33
java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar SE SRR17673025.fastq.gz SRR17673025.trim.fastq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:70 -phred33
java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar SE SRR17673026.fastq.gz SRR17673026.trim.fastq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:70 -phred33
java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar SE SRR17673027.fastq.gz SRR17673027.trim.fastq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:70 -phred33

java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar SE SRR17673030.fastq.gz SRR17673030.trim.fastq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:70 -phred33
java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar SE SRR17673031.fastq.gz SRR17673031.trim.fastq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:70 -phred33
java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar SE SRR17673032.fastq.gz SRR17673032.trim.fastq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:70 -phred33
java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar SE SRR17673033.fastq.gz SRR17673033.trim.fastq.gz ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:70 -phred33

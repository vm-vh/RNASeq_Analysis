
# This script is for importing RNASeq data into R as well as filtering and normalizing it.

library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science
library(tximport) # package to import Kallisto results into R
library(ensembldb) # helps with ensembl
library(biomaRt) # annotations
library(edgeR) # package for differential expression analysis (used for the DGEList object and for normalization methods)
library(matrixStats)
library(cowplot)


# ----------- Import and format data -----------------------------------------------------------------------------------------------
myMart <- useMart(biomart="plants_mart", host="https://plants.ensembl.org") # access ensemble plant genomes
cof.anno <- useMart(biomart="plants_mart", dataset = "ccanephora_eg_gene", host="https://plants.ensembl.org") # load ensembl annotations for c. canephora

#formatting
Tx.cof <- getBM(attributes=c('ensembl_transcript_id',
                             'ensembl_gene_id', 'description'),
                              mart = cof.anno)
Tx.cof <- as_tibble(Tx.cof) 
Tx.cof <- dplyr::rename(Tx.cof, target_id = ensembl_transcript_id, gene_name = ensembl_gene_id)
Tx.cof <- dplyr::select(Tx.cof, "target_id", "gene_name")

# import a study design file
targets <- read_tsv("study_design.txt")

# import kallisto transcript counts
path <- file.path("kallisto_output",targets$sra_accession, "abundance.tsv") # set file paths to mapped data
all(file.exists(path)) # check to see if path worked, should be TRUE before proceeding
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx.cof, # mapps of transcript IDs to gene IDs 
                     txOut = FALSE, # data represented at gene level instead of at transcript level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE) 



# ----------- Filtering and normalization -----------------------------------------------------------------------------------------------

# visualization of abundances
myTPM <- Txi_gene$abundance #abundances in transcripts per million
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM)) # compute standard statistics for the abundance data
ggplot(myTPM.stats) +
  aes(x = SD, y = MED) +
  geom_smooth(method=lm) +
  geom_hex(show.legend = T, bins=25) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized data",
       caption="DIYtranscriptomics - Spring 2020") +
  theme_bw()

# visualization of counts
sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts) # creates a digital gene expression list object
log2.cpm <- cpm(myDGEList, log=TRUE) # use 'cpm' function to get counts per million
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") # formatting
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, 
                                  cols = L42:S31, 
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable storing the data

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "samples",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# filtering
table(rowSums(myDGEList$counts==0)==8) # genes not expressed in all samples
cpm <- cpm(myDGEList) # convert to counts per million
keepers <- rowSums(cpm>1) >= 4 # removes genes with a total count below 1 CPM over any 4 samples
                              # 4 was chosen as the number of samples in the smallest group of comparison (4) 
myDGEList.filtered <- myDGEList[keepers,] 
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE) # convert to cpm and log2 scale
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID") # format
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df,
                                           cols = L42:S31, 
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "samples",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# normalization
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM") # calculate normalization factors
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE) # convert to cpm and log2 scale
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID") # format
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = L42:S31, 
                                                names_to = "samples", # name of that new variable
                                                values_to = "expression") # name of new variable storing the data

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "samples",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3) # plot all 3 violin plots

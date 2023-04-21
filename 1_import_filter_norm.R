# This script is for importing RNASeq data into R as well as filtering and normalizing it.

Sys.unsetenv("R_LIBS_USER")
.libPaths(paste(getwd(), "temp/RLibrary", sep="/"))

library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science
library(tximport) # package for getting Kallisto results into R
library(ensembldb) # helps deal with ensembl
library(biomaRt) # annotations
library(edgeR) # package for differential expression analysis, here only used for the DGEList object and for normalization methods
library(matrixStats)
library(cowplot)

# step 1: import studydesign and kallisto transcript counts and turn into list with gene IDs through ensmbl annotations
targets <- read_tsv("study_design.txt")# read in study design

path <- file.path("kallisto_output",targets$sra_accession, "abundance.tsv") # set file paths to mapped data
all(file.exists(path)) # check to see if path worked

myMart <- useMart(biomart="plants_mart", host="https://plants.ensembl.org") # marts for plant genomes
cof.anno <- useMart(biomart="plants_mart", dataset = "ccanephora_eg_gene", host="https://plants.ensembl.org") # load ensembl annotations for c. canephora
# cof.attributes <- listAttributes(cof.anno) 
Tx.cof <- getBM(attributes=c('ensembl_transcript_id',
                             'ensembl_gene_id', 'description'),
                mart = cof.anno) # turn it into a table
Tx.cof <- as_tibble(Tx.cof) 
Tx.cof <- dplyr::rename(Tx.cof, target_id = ensembl_transcript_id, gene_name = ensembl_gene_id) # rename the two columns from biomart
Tx.cof <- dplyr::select(Tx.cof, "target_id", "gene_name") # set transcript ID as the first column

# import Kallisto transcript counts
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx.cof, # Mapping of transcript IDs to gene IDs 
                     txOut = FALSE, # data represented at gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE) 



# step 2: filter and normalize data, tidy up gene expression data
# visualization of abundances
myTPM <- Txi_gene$abundance
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_smooth(method=lm) +
  geom_hex(show.legend = T, bins=25) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized data",
       caption="DIYtranscriptomics - Spring 2020") +
  theme_bw()

sampleLabels <- targets$sample

# looking at counts
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE) # use 'cpm' function to get counts per million
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, 
                                  cols = L42:S31, 
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

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
table(rowSums(myDGEList$counts==0)==8)
# FALSE: 23184 have read counts
# TRUE: 2390 have NO read counts

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1) >= 4 #user defined - number of samples in the smallest group of comparison (4)
# FALSE: 18668 have read counts above the cutoff
myDGEList.filtered <- myDGEList[keepers,]
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
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
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# normalization
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = L42:S31, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)
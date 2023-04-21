# This script is for clustering RNASeq data, 
# Please run the import_filter_norm.R script BEFORE this one to ensure all necessary variables are defined.

Sys.unsetenv("R_LIBS_USER")
.libPaths(paste(getwd(), "temp/RLibrary", sep="/"))

# load all required libraries
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science
library(tximport) # package for getting Kallisto results into R
library(ensembldb) # helps deal with ensembl
library(biomaRt) # annotations

library(edgeR) # package for differential expression analysis, here only used for the DGEList object and for normalization methods
library(matrixStats)
library(cowplot)

library(DT) # interactive and searchable tables of GSEA results
library(gt)
library(plotly)

library(limma)

library(gplots) #for heatmaps
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
# library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
# library(enrichplot) # great for making the standard GSEA enrichment plots
# library(qusage) # Quantitative Set Analysis for Gene Expression
library(heatmaply)

#step 3: data exploration
group <- targets$group
group <- factor(group)

# hierarchical clustering
distance <- dist(t(log2.cpm.filtered.norm), method = "euclidean") 
clusters <- hclust(distance, method = "complete")
c_euc <- plot(clusters, labels=sampleLabels)

distance <- dist(t(log2.cpm.filtered.norm), method = "maximum") 
clusters <- hclust(distance, method = "complete") 
c_max <- plot(clusters, labels=sampleLabels)

# PCA
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) # the percentage variance explained by each PC
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)

# calculate averages for each gene
mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    leaf.AVG = (L42 + L41 + L32 + L31)/4, 
                    stem.AVG = (S42 + S41 + S32 + S31)/4,
                    #now make columns comparing each of the averages above that you're interested in
                    LogFC = (stem.AVG - leaf.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)

datatable(mydata.df[,c(1,10:12)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         #dom = "Blfrtip", 
                         #buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))

myplot <- ggplot(mydata.df) +
  aes(x=leaf.AVG, y=stem.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("stem vs. leaf") +
  theme_bw()

ggplotly(myplot)


# step 4 - identify differentially expressed genes (DEGs) and differential transcript usage (DTU)
group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Model mean-variance trend and fit linear model to data
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE) #model the mean-variance relationship
fit <- lmFit(v.DEGList.filtered.norm, design) # fit a linear model to your data
contrast.matrix <- makeContrasts(infection = stem - leaf, # TODO: change!!
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix) # extract the linear model fit
ebFit <- eBayes(fits) # get bayesian stats for your linear model fit

# TopTable to view DEGs
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")
gt(myTopHits.df)

vplot <- ggplot(myTopHits) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", myTopHits$geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Coffea canephora",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:9), digits=2)

# step X modules
library(tidyverse)
library(gplots)
library(RColorBrewer)
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))

modulePick <- 2 
myModule_up <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_up <- hclust(as.dist(1-cor(t(myModule_up), method="pearson")), method="complete") 

heatmap.2(myModule_up, 
          Rowv=as.dendrogram(hrsub_up), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

modulePick <- 1 
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_down <- hclust(as.dist(1-cor(t(myModule_down), method="pearson")), method="complete") 

heatmap.2(myModule_down, 
          Rowv=as.dendrogram(hrsub_down), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

# step 5
gost.res_up <- gost(rownames(myModule_up), organism = "ccanephora", correction_method = "fdr")
gostplot(gost.res_up, interactive = T, capped = T) #set interactive=FALSE to get plot for publications
gost.res_down <- gost(rownames(myModule_down), organism = "ccanephora", correction_method = "fdr")
gostplot(gost.res_down, interactive = T, capped = T) #set interactive=FALSE to get plot for publications

hs_gsea_c2 <- msigdbr(species = "Coffea canephora", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 

# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)

# view results as an interactive table
datatable(myGSEA.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)
# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res, 
          geneSetID = 47, #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = myGSEA.res$Description[47]) #can also turn off this title

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()
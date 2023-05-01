
# This script is for clustering RNASeq data, identifying differentially expressed genes and preforming GO enrichment analysis.
# Please run the import_filter_norm.R script BEFORE this one to ensure all necessary variables are defined.


# load all required libraries
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science
library(DT) # interactive tables
library(plotly) # interactive plots
library(gt) # A layered 'grammar of tables'
library(limma) # package for differential gene expression using linear modeling
library(edgeR) # package for differential expression analysis
library(gplots) 
library(gprofiler2) # tools for accessing the GO enrichment results using g:Profiler web resources
library(RColorBrewer)


# ----------- Hierarchical clustering -----------------------------------------------------------------------------------------------
distance_euc <- dist(t(log2.cpm.filtered.norm), method = "euclidean")  # compute distances based on the normalized and filtered data
clusters_euc <- hclust(distance_euc, method = "complete") # cluster the data based on distance
c_euc <- plot(clusters_euc, labels=sampleLabels) # plot the dendrogram

distance_max <- dist(t(log2.cpm.filtered.norm), method = "maximum") # compute distances based on the normalized and filtered data
clusters_max <- hclust(distance_max, method = "complete") # cluster the data based on distance
c_max <- plot(clusters_max, labels=sampleLabels) # plot the dendrogram


# ----------- PCA -----------------------------------------------------------------------------------------------
group <- factor(targets$group)
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T) # preform PCA

pc.per <- round((pca.res$sdev^2)/sum((pca.res$sdev^2))*100, 1) # the percentage variance explained by each PC (sdev^2 captures these eigenvalues from the PCA result)
pc.per.df <- data.frame(PC= paste0("PC",1:8), pc.per) # frmat as a dataframe
pc.per.df %>%
  ggplot(aes(x=PC,y=pc.per, group=1))+
  geom_point(size=8)+
  geom_line()+
  labs(title="PCA") +
  xlab("Principal Components") +
  ylab("% of Variance Explained")


colors <- c("green", "darkgreen") # define the colors for each group
pca.res.df <- pca.res$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)

pca.pivot <- pivot_longer(pca.res.df, 
                          cols = PC1:PC4, 
                          names_to = "PC", 
                          values_to = "loadings") 

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) +
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  scale_fill_manual(values = colors) + # use a manual color scale
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  scale_color_manual(values = colors) + # use a manual color scale
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)


# ----------- Gene expression averages -----------------------------------------------------------------------------------------------
mydata.df <- mutate(log2.cpm.filtered.norm.df, 
                    leaf.AVG = (L42 + L41 + L32 + L31)/4, # calculate averages for leaf
                    stem.AVG = (S42 + S41 + S32 + S31)/4, # calculate averages for stem
                    LogFC = (stem.AVG - leaf.AVG)) %>% # create columns comparing each of the averages above
             mutate_if(is.numeric, round, 2)
write_tsv(mydata.df[,c(1,10:12)], "avg.tsv")

datatable(mydata.df[,c(1,10:12)], # view as a table
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

myplot <- ggplot(mydata.df) +
  aes(x=leaf.AVG, y=stem.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("stem vs. leaf") +
  theme_bw()

ggplotly(myplot)


# ----------- Identification of differentially expressed genes (DEGs) --------------------------------------------------------------------------------------------
group <- factor(targets$group)
design <- model.matrix(~0 + group) # create a design matrix based on the tissues
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE) # model the mean-variance relationship
fit <- lmFit(v.DEGList.filtered.norm, design) # fit a linear model to the data
contrast.matrix <- makeContrasts(differance = stem - leaf, # create a contrast matrix based on the tissue differences
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix) # extract the linear model fit
ebFit <- eBayes(fits) # get bayesian stats for the linear model fit
write.fit(ebFit, file="lmfit_results.txt")

# View DEGs
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=18668, sort.by="logFC") # number is set to all DEGs
myTopHits.df <- myTopHits %>% as_tibble(rownames = "geneID")
myTopHitsgt <- data.frame(gt(myTopHits.df))
write_tsv(myTopHitsgt,"myTopHits.txt")

# creates a volcano plot of the DEGs
vplot <- ggplot(myTopHits) + 
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", rownames(myTopHits))) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="green", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="darkgreen", linewidth=1) +
  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="green") +
  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="darkgreen") +
  labs(title="Volcano plot",
       subtitle = "C. canephora",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1) # output -1 or 1: t-statistic for gene is classified as significant (0 = not significant)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,] # E: numeric matrix of normalized expression values on the log2 scale
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
write_tsv(diffGenes.df,"DiffGenes.tsv")

datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in c. canephora',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:9), digits=2)


# ----------- Gene Ontology (GO) enrichment using gProfiler2 -----------------------------------------------------------------------------------------------
myTop100Hits <- topTable(ebFit, adjust ="BH", coef=1, number=100, sort.by="logFC") # pick the top 100 genes by logFC value for carrying out GO enrichment analysis
gost.res <- gost(rownames(myTop100Hits), organism = "ccanephora", correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res, interactive = T, capped = F) # interactive manhattan plot of enriched GO terms

#save gost plot and table
publish_gostplot(
   gostplot(gost.res, interactive = F, capped = F), # static gostplot
   highlight_terms = c("GO:0005975", "GO:0048046", "GO:0030145", "GO:0003824", "GO:0045735"), # highlight top 7 lowest p-values
   filename = "gostplot100.png",
   width = NA,
   height = NA)

publish_gosttable(
   gost.res$result[order(gost.res$result$p_value, decreasing = F),],
   highlight_terms = NULL,
   use_colors = TRUE,
   show_columns = c("source", "term_name", "term_size", "intersection_size"),
   filename = "gosttable100.pdf",
   ggplot=TRUE)

# create tissue modules based on pearson correlation
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") # cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") # cluster columns by pearson correlation
module.assign <- cutree(clustRows, k=2) 

modulePick <- 2 # pick the leaf module
myModule_leaf <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
length(rownames(myModule_leaf))

gost.res.leaf <- gost(rownames(myModule_leaf), organism = "ccanephora", correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res.leaf, interactive = T, capped = F) # interactive Manhattan plot of enriched GO terms

# save gost plot and table
publish_gostplot(
  gostplot(gost.res.leaf, interactive = F, capped = F), # static gostplot
  highlight_terms = c("GO:0009536", "GO:0009579", "GO:0009507", "GO:0015979", "GO:0034357"), # highlight top 5
  filename = "gostplot_leaf.png",
  width = NA,
  height = NA)

publish_gosttable(
  gost.res.leaf$result[order(gost.res.leaf$result$p_value, decreasing = F),],
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "gosttable_leaf.pdf",
  ggplot=TRUE)


modulePick <- 1 # pick the stem module
myModule_stem <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
length(rownames(myModule_stem))

gost.res.stem <- gost(rownames(myModule_stem), organism = "ccanephora", ordered_query = T, correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res.stem, interactive = T, capped = F) # interactive Manhattan plot of enriched GO terms

# save gost plot and table
publish_gostplot(
   gostplot(gost.res.stem, interactive = F, capped = F), # static gostplot
   highlight_terms = c("GO:0008017", "GO:0015631", "GO:0008092", "GO:0003777", "GO:0003774"), # highlight top 7
   filename = "gostplot_stem.png",
   width = NA,
   height = NA)

publish_gosttable(
   gost.res.stem$result[order(gost.res.stem$result$p_value, decreasing = F),],
   highlight_terms = NULL,
   use_colors = TRUE,
   show_columns = c("source", "term_name", "term_size", "intersection_size"),
   filename = "gosttable_stem.pdf",
   ggplot= T)


# ----------- Heatmaps -----------------------------------------------------------------------------------------------
myheatcolors <- rev(brewer.pal(name="Spectral", n=11)) # define a color palette
module.color <- rainbow(length(unique(module.assign)), start=0.3, end=0.2) 
module.color <- module.color[as.vector(module.assign)] # assign module colors

# generate a heatmap for all DEGs
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))

# generate a heatmap for the leaf module
modulePick <- 2 
hrsub_leaf <- hclust(as.dist(1-cor(t(myModule_leaf), method="pearson")), method="complete") 
heatmap.2(myModule_leaf, 
          Rowv=as.dendrogram(hrsub_leaf), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

# generate a heatmap for the stem module
modulePick <- 1
hrsub_stem <- hclust(as.dist(1-cor(t(myModule_stem), method="pearson")), method="complete") 
heatmap.2(myModule_stem, 
          Rowv=as.dendrogram(hrsub_stem), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

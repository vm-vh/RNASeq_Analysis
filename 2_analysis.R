# This script is for clustering RNASeq data, 
# Please run the import_filter_norm.R script BEFORE this one to ensure all necessary variables are defined.

Sys.unsetenv("R_LIBS_USER")
.libPaths(paste(getwd(), "temp/RLibrary", sep="/"))

# load all required libraries
library(tidyverse) # you're familiar with this from the past two lectures
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR) # package for differential expression analysis
library(reshape2)
library(heatmaply) # for heatmaps
library(gplots) 
library(gprofiler2) # tools for accessing the GO enrichment results using g:Profiler web resources
library(RColorBrewer)

#step 3: data exploration

# hierarchical clustering
distance_euc <- dist(t(log2.cpm.filtered.norm), method = "euclidean") 
clusters_euc <- hclust(distance_euc, method = "complete")
c_euc <- plot(clusters_euc, labels=sampleLabels)

distance_max <- dist(t(log2.cpm.filtered.norm), method = "maximum") 
clusters_max <- hclust(distance_max, method = "complete") 
c_max <- plot(clusters_max, labels=sampleLabels)

# PCA
group <- targets$group
group <- factor(group)
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.per <- round((pca.res$sdev^2)/sum((pca.res$sdev^2))*100, 1) # the percentage variance explained by each PC (sdev^2 captures these eigenvalues from the PCA result)
pc.per.df <- data.frame(PC= paste0("PC",1:8), pc.per)
pc.per.df %>%
  ggplot(aes(x=PC,y=pc.per, group=1))+
  geom_point(size=8)+
  geom_line()+
  labs(title="PCA") +
  xlab("Principal Components") +
  ylab("% of Variance Explained")
  
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
# change to paired?
# design2 <- model.matrix(~block + group)
# colnames(design2) <- levels(group)

# Model mean-variance trend and fit linear model to data
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE) # model the mean-variance relationship
fit <- lmFit(v.DEGList.filtered.norm, design) # fit a linear model to your data
contrast.matrix <- makeContrasts(differance = stem - leaf,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix) # extract the linear model fit
ebFit <- eBayes(fits) # get bayesian stats for your linear model fit
# write.fit(ebFit, file="lmfit_results.txt")

# TopTable to view DEGs
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=18668, sort.by="logFC")
myTopHits.df <- myTopHits %>% as_tibble(rownames = "geneID")
myTopHitsgt <- data.frame(gt(myTopHits.df))
# write_tsv(myTopHitsgt,"myTopHits.txt")

vplot <- ggplot(myTopHits) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", rownames(myTopHits))) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  # annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  # annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "C. canephora",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1) # output -1 or 1: t-statistic for gene is classified as significant
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,] # E = numeric matrix of normalized expression values on the log2 scale
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in c. canephora',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:9), digits=2)
# write_tsv(diffGenes.df,"DiffGenes.txt")


# step 5

# Gene Ontology (GO) enrichment using gProfiler2
# using myTopHits
myTop100Hits <- topTable(ebFit, adjust ="BH", coef=1, number=100, sort.by="logFC") # pick the top genes for carrying out GO enrichment analysis

# myTop5Hits_logFC <- topTable(ebFit, adjust ="BH", coef=1, number=5, sort.by="logFC")
# myTop5Hits_logFC
# myTop5Hits_p <- topTable(ebFit, adjust ="BH", coef=1, number=5, sort.by="P")
# myTop5Hits_p

gost.res <- gost(rownames(myTop100Hits), organism = "ccanephora", correction_method = "fdr") # run GO enrichment analysis
ordered.gost.res <- gost.res$result[order(gost.res$result$p_value, decreasing = F),][1:7]
gostplot(gost.res, interactive = T, capped = F) # interactive manhattan plot of enriched GO terms

# save gost plot and table
# publish_gostplot(
#    gostplot(gost.res, interactive = F, capped = F), # static gostplot
#    highlight_terms = c("GO:0005975", "GO:0048046", "GO:0030145", "GO:0003824", "GO:0045735"), # highlight top 7 lowest p-values
#    filename = "gostplot100.png",
#    width = NA,
#    height = NA)
# 
# publish_gosttable(
#    gost.res$result[order(gost.res$result$p_value, decreasing = F),],
#    highlight_terms = NULL,
#    use_colors = TRUE,
#    show_columns = c("source", "term_name", "term_size", "intersection_size"),
#    filename = "gosttable100.pdf",
#    ggplot=TRUE)


# step X

clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") # cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2) # assign 

modulePick <- 2 
myModule_up <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
length(rownames(myModule_up)) # 3355

gost.res.up <- gost(rownames(myModule_up), organism = "ccanephora", correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res.up, interactive = T, capped = F) # interactive manhattan plot of enriched GO terms

# save gost plot and table
# publish_gostplot(
#   gostplot(gost.res.up, interactive = F, capped = F), # static gostplot
#   highlight_terms = c("GO:0009536", "GO:0009579", "GO:0009507", "GO:0015979", "GO:0034357"), # highlight top 5
#   filename = "gostplot_up.png",
#   width = NA,
#   height = NA)
# 
# publish_gosttable(
#   gost.res.up$result[order(gost.res.up$result$p_value, decreasing = F),],
#   highlight_terms = NULL,
#   use_colors = TRUE,
#   show_columns = c("source", "term_name", "term_size", "intersection_size"),
#   filename = "gosttable_up.pdf",
#   ggplot=TRUE)

modulePick <- 1 
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
length(rownames(myModule_down)) # 3530

gost.res.down <- gost(rownames(myModule_down), organism = "ccanephora", ordered_query = T, correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res.down, interactive = T, capped = F) # interactive manhattan plot of enriched GO terms

# save gost plot and table
# publish_gostplot(
#    gostplot(gost.res.down, interactive = F, capped = F), # static gostplot
#    highlight_terms = c("GO:0008017", "GO:0015631", "GO:0008092", "GO:0003777", "GO:0003774"), # highlight top 7
#    filename = "gostplot_down.png",
#    width = NA,
#    height = NA)
# 
# publish_gosttable(
#    gost.res.down$result[order(gost.res.down$result$p_value, decreasing = F),],
#    highlight_terms = NULL,
#    use_colors = TRUE,
#    show_columns = c("source", "term_name", "term_size", "intersection_size"),
#    filename = "gosttable_down.pdf",
#    ggplot= T)



# heatmaps
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
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

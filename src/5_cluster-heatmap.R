#------------------------------------------------------------
#                RNA seq analysis workflow
#                    6_cluster_heatmap
#                   Ryan Berger 6-5-18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: Make cluster heatmap from differential expression data
#  Input: Differential expression results object (res), rlog transformed (rld)
#  Output: heatmap showing fold change across samples
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------

## Set directory
## setwd('~/Desktop/RNAseq stuff/HIOs-rnaseq/src')

getwd() ## Should be '../HIOs-rnaseq/src'

#------------------------------------------------------------
#                    6_cluster_heatmap

# Heatmap of top variance genes
library(genefilter)
load('../results/DESeq2/8h-exp-samples/dds.Rdata')

## rlog transformation of dds 
# rld <- rlog(dds, blind = FALSE) 
# save(rld, file = "../results/DESeq2/8h-exp-samples/rld.Rdata")

## If rlog already saved
load("../results/DESeq2/8h-exp-samples/rld.Rdata")
head(assay(rld))

## Change format of code name for labels in heatmap
colData(rld)$code_name <- gsub('Styphimurium','STM',colData(rld)$code_name) %>% 
  gsub('Senteritidis', 'SE', .) %>% 
  gsub('SEnt', 'SE', .)
colData(rld)$code_name  

#------------------------------------------------------------

## Make heatmap of averages of all 4 replicates

num1 <- 35  # number of genes in heatmap. (e.g. top n variance genes)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), num1) # Define top variance genes
mat <- assay(rld)[topVarGenes, ]  # Matrix of the top variance genes from the set
mat <- mat - rowMeans(mat)  # expression is relative to mean of all samples. (fold over mean of all samples)

## Calculate average of each sample across 4 replicates
new.df <- as.data.frame(mat)
colnames(new.df)
new.df$PBS <- apply(new.df[,1:4], 1, mean)
new.df$STM <- apply(new.df[,5:8],1,mean)
new.df$SE <- apply(new.df[,9:12],1,mean)
new.df$`PBS+PMN` <- apply(new.df[,13:16],1,mean)
new.df$`STM+PMN` <- apply(new.df[,17:20],1,mean)
new.df$`SE+PMN` <- apply(new.df[,21:24],1,mean)
new.df <- as.matrix(new.df[,c(25,28,26,27,29,30)])

pheatmap(new.df, 
         cutree_rows = 3,
         cutree_cols = 3,
         treeheight_row = 10,
         treeheight_col = 10,
         cluster_cols = T,
         show_colnames = T,
         fontsize = 8.5)


## save pdf of plot
heatmap <- recordPlot() # Run after plot is finished
pdf('../img/heatmap.pdf', width = 4, height = 5, onefile = FALSE)
heatmap
dev.off()

#------------------------------------------------------------

## Old format of heatmap with all samples as individual boxes

# Cluster heatmap of all samples in dds
num1 <- 35  # number of genes in heatmap. (e.g. top n variance genes)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), num1) # Define top variance genes
mat <- assay(rld)[topVarGenes, ]  # Matrix of the top variance genes from the set
mat <- mat - rowMeans(mat)  # expression is relative to mean of all samples. (fold over mean of all samples)
anno <- as.data.frame(colData(rld)$code_name)  # sample IDs for annotation
colnames(anno) <- 'Sample'
rownames(anno) <- colData(rld)$short_name


## Manual colors for column IDs
ann_colors <- brewer.pal(8, 'Paired')[3:8]
names(ann_colors) <- unique(colData(rld)$code_name)[c(1,4,2,5,3,6)] ## Change legend orders
ann_colors <- list(Sample = ann_colors)


# Heatmap
pheatmap(mat, 
         annotation_col = anno,
         cutree_rows = 3,
         cutree_cols = 3,
         treeheight_row = 25,
         treeheight_col = 25,
         cluster_cols = T,
         show_colnames = F,
         annotation_colors = ann_colors)

?pheatmap
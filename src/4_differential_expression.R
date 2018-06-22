#------------------------------------------------------------
#                RNA seq analysis workflow
#                3_differential_expression
#                   Ryan Berger 6--18
#
#   RNA seq workflow starting with kallisto aligned reads
#
#  Purpose: compare gene expression between two samples
#  Input: DESeq dataset (dds)
#  Output: results (res) and dataframe (res.df) 
#
#
# Source: https://www.bioconductor.org/help/workflows/rnaseqGene/
#------------------------------------------------------------

## Set directory
## setwd('~/Desktop/RNAseq stuff/HIOs-rnaseq/src')

getwd() ## Should be '../HIOs-rnaseq/src'
results.dir <- "../results/DESeq2/8h-exp-samples"
#------------------------------------------------------------
#                3_differential_expression

## Load data
load('../results/DESeq2/8h-exp-samples/dds.Rdata')
dds <- DESeq2::DESeq(dds)

## Set up function for making and saving differential expression files

diff_expression <- function(dds, sample1,sample2){
  contrast <- c('new_code', sample1, sample2) # factor name, numerator condition, denominator condition
res <- results(dds, contrast = contrast)

## Add annotation - symbol and entrezID
library(AnnotationDbi)
library(org.Hs.eg.db)

# Add column for gene symbol
res$symbol <- rownames(res)
# Add column for gene Entrez ID
res$entrez <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = 'ENTREZID',
                     keytype = 'SYMBOL',
                     multiVals = 'first')
# Add column for gene name
res$name <- mapIds(org.Hs.eg.db,
                   keys = rownames(res),
                   column = 'GENENAME',
                   keytype = 'SYMBOL',
                   multiVals = 'first')
# Make results dataframe
res.df <- as.data.frame(res)
res.df <- res.df[order(res.df$padj),]

# Save output file
write.csv(res.df, file = paste0(results.dir,'/diff_expression-',sample1,'-over-', sample2,'.csv'))
print(paste(sample1, 'over', sample2, 'done'))
}

## Run differential expression for all sample+PMN over sample
unique(colData(dds)$new_code)

diff_expression(dds, 'STM+PMNs', 'STM')
diff_expression(dds, 'SE+PMNs', 'SE')
diff_expression(dds, 'PBS+PMNs', 'PBS')

## Clean up
rm(list = ls())
#------------------------------------------------------------
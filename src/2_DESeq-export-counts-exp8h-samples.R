## R script to export kallisto results to matrix with DESeq2
## Ryan Berger
## 6-4-18
## -----------------------------------------------------------------------------

getwd()  ## Should be ../HIOs-rnaseq/src

## Differential expression of kallisto results with DESeq2
kallisto.results.dir <- "../results/kallisto-Run_1822/"

## create directory to deposit results
results.dir <- "../results/DESeq2/8h-exp-samples"
dir.create(path = results.dir, recursive = TRUE)

## read in table with sample metadata
## load using the UM core provided sample submission form
samples <- readr::read_csv(file = "../results/kallisto-Run_1822/sample_key3.csv") %>% 
  filter(hr == 8 & code_name != 'HIOs+PMNs' & code_name != 'PMN')

# Set up path to read files into tximport
files <- file.path(kallisto.results.dir, samples$file_name, 'abundance.tsv')
# Add sample IDs to files
names(files) <- samples$short_name
files

## check that all files are found
if (all(file.exists(files)) == FALSE) {
    print("kallisto files not found")
} else{
  print('all kallisto files are found')
}

## associate transcripts with gene IDs
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

# Create a dataframe of transcripts from ensembl human genome
Tx <- transcripts(edb,
                  columns = c(listColumns(edb , "tx"), "gene_name"),
                  return.type = "DataFrame")
head(Tx)
# assign columns 1 (transcript ID) and 9 (gene name) to dataframe tx2gene
tx2gene <- Tx[,c(1,9)]
head(tx2gene)

## import kallisto data and generate count dataframe (dds)
## http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
library(readr)
txi <- tximport::tximport(files,
                          type = "kallisto",
                          tx2gene = tx2gene)

## export abundance counts
write.csv(txi$abundance, file = file.path(results.dir, "complete_dataset_txi.csv"))

library(DESeq2)
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

dds <- DESeq2::DESeqDataSetFromTximport(txi,
                                        colData = samples,
                                        design = ~ new_code)
## pre-filter out counts < 1
dds <- dds[rowSums(counts(dds)) > 1, ]

## write out normalized expression counts
dds <- DESeq2::estimateSizeFactors(dds)
ddscounts <- DESeq2::counts(dds, normalized = TRUE)

## write expression matrix to file
write.csv(ddscounts, file =  file.path(results.dir, "complete-dataset_DESeq2-normalized-counts.csv"))
save(dds, file = file.path(results.dir, "dds.Rdata"))

## clear working memory
rm(list = ls())

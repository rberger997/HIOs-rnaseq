## R Script to generate principle component analysis of HIO response
## Ryan Berger
## 6-4-18
## Adapted from script by David R. Hill and Ryan Berger
## -----------------------------------------------------------------------------


# Load packages
library(magrittr)
library(ggplot2)
library(RColorBrewer)
source("ggplot2-themes.R")

## Load dds object
load('../results/DESeq2/all-8h-samples/dds.Rdata')

## Perform variance stabilizing transformations on dds using rlog
## use argument blind = FALSE when multiple replicates are present
rld <- DESeq2::rlog(dds, blind = FALSE)

## Save rlog object
save(rld, file = '../results/DESeq2/all-8h-samples/rld.Rdata')

## Load rlog object
load('../results/DESeq2/all-8h-samples/rld.Rdata')


## Set order of cell line id
colData(rld)$cell_line_id <- factor(colData(rld)$cell_line_id)
colData(rld)$cell_line_id <- factor(colData(rld)$cell_line_id, 
                                    levels(colData(rld)$cell_line_id)[c(1,3,2,4)])

## Set order of injection
colData(rld)$injection_id <- factor(colData(rld)$injection_id)
colData(rld)$injection_id <- factor(colData(rld)$injection_id,
                                    levels(colData(rld)$injection_id)[c(1,2,4,3)])

# PCA plot (use just only for PC variance estimates)
pca <- plotPCA(rld, intgroup = c('injection_id', 'cell_line_id'))

# Get PCA data
pca.df <- plotPCA(rld, intgroup = c('cell_line_id', 'injection_id'), returnData = TRUE)

## Make PCA plot
plot <- ggplot(data = pca.df, aes(x = PC1, y = PC2))+ 
geom_point(shape = 21, stroke = 1, 
           aes(fill = as.factor(cell_line_id), 
               color = as.factor(injection_id)), 
           size = 10) +
  theme1 + 
  scale_y_continuous(limits = c(-25,25), breaks = seq(-25, 25, by = 10)) +
  scale_fill_brewer(palette = "Set1", name = 'Cell type') +
  scale_color_brewer(palette = "Greys", name = 'Injection', direction = 1) +    
  theme(legend.position = "right") +
  geom_hline(yintercept = 0,
             size = 1, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0,
             size = 1, linetype = "dashed", color = "grey70") +
  coord_fixed(ratio = 1) +
  xlab(pca$labels$x) + #pull variance estimates from al. plotPCA call
  ylab(pca$labels$y) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +  # Move y axis
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))  # Move x axis

plot

## save png of plot
png(filename = "../img/pca.png",
    width = 1600, height = 500)
print(plot)
dev.off()

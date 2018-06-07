## R script to generate facet volcano plot of each sample+PMN over sample 
## differential expression
## Ryan Berger
## 6-7-18

## Modified from original version by David R. Hill

##          Volcano plot
## -----------------------------------------------------------------------------
## Set directory
setwd('~/Desktop/RNAseq stuff/HIOs-rnaseq/src')

getwd() ## Should be '../HIOs-rnaseq/src'
results.dir <- "../results/DESeq2/8h-exp-samples"
## -----------------------------------------------------------------------------
library(magrittr)
## data import
data.dir <- "../results/DESeq2/8h-exp-samples"

## read in individual diff expression files
pbs <- read.csv(file = file.path(data.dir,
                                         "diff_expression-PBS+PMNs-over-PBS.csv"))
pbs$samp <- 'PBS'

stm <- read.csv(file = file.path(data.dir,
                                 "diff_expression-Styphimurium+PMNs-over-Styphimurium.csv"))
stm$samp <- 'STM'
se <- read.csv(file = file.path(data.dir,
                                 "diff_expression-SEnt+PMNs-over-Senteritidis.csv"))
se$samp <- 'SE'

## combine into single dataframe
data <- rbind(pbs, stm, se) %>% 
    dplyr::rename(SYMBOL = X)

## create status catergory for assigning colors
data$status <- ifelse(data$padj > 10^-5 | is.na(data$padj), "a",
                    ifelse(data$log2FoldChange > 0, "b", "c"))
## sort by status
data <- data[order(data$status),]


## pre-reqs
library(ggplot2)
source("ggplot2-themes.R")

## Set order of samples
data$samp <- factor(data$samp, levels = c('PBS', 'STM', 'SE'))

## set up facet plot
plot <- ggplot(data = data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(shape = 21, size = .25, aes(fill = status, color = status)) +
  facet_grid(samp ~ .) +
  ## set axis limits and labels
  xlim(c(-4.25,4.25)) +
  #ylim(c(0,25)) +
  xlab(expression(paste("log"[2],"FC(Sample+PMN/Sample)"))) +
  ylab(expression(paste("-log"[10],"(adj. P-value)"))) +
  ## colors and theme
  scale_fill_manual(values = c("grey70", color.set[1], color.set[2])) +
  scale_color_manual(values = c("grey70", color.set[1], color.set[2])) +
  theme(
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 8.5),
    axis.title = element_text(size = 8.5),
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid = element_line(size = .25),
    panel.border = element_rect(fill = NA,
                                color = "grey85",
                                size = 1)				    
  )
plot

## save pdf of plot
pdf('../img/volcano-facet.pdf', width = 3, height = 6, onefile = FALSE)
plot
dev.off()


## -----------------------------------------------------------------------------
## David's PNG format
## set up facet plot
plot <- ggplot(data = data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(shape = 21, size = 2, aes(fill = status, color = status)) +
    facet_grid(. ~ samp) +
    ## set axis limits and labels
    xlim(c(-2,2)) + ylim(c(0,10)) +
    xlab(expression(paste("log"[2],"FC(Sample+PMN/Sample)"))) +
    ylab(expression(paste("-log"[10],"(adj. P-value)"))) +
    ## colors and theme
    scale_fill_manual(values = c("grey70", color.set[1], color.set[2])) +
    scale_color_manual(values = c("grey70", color.set[1], color.set[2])) +
    theme(
        strip.text = element_text(size = 48),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 24),
	legend.position = "none",
	panel.background = element_blank(),
        panel.border = element_rect(fill = NA,
                                    color = "grey85",
                                    size = 1)				    
    )

## open plot device
#png(filename = "../img/volcano-plot.png",
#    width = 1000, height = 1000)
#print(plot)
#dev.off()
                               

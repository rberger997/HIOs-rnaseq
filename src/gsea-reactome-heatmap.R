## R script to generate a plot of GSEA results over time
## David R. Hill
## -----------------------------------------------------------------------------

library(magrittr)

## load in datasets
results <- readr::read_csv(file = "../results/GSEA/combined_GSEA_results.csv") %>%
    dplyr::rename(hr = comparison)
results <- filter(results, database == 'GO')
results <- filter(results, pvalue < 0.05 & abs(NES) > 2.5)

## save list of pathways used in the plot
write.csv(unique(results$Description), file = "../results/GSEA/plot_pathways.csv")

## load in categorized pathways and join to results
results <- readr::read_csv(file = "../results/GSEA/plot_pathways_categorized.csv") %>%
    dplyr::rename(Description = x) %>%
    dplyr::right_join(results, by = 'Description')

results$Category <- 'a'
## Sort pathways for more coherent plotting
results <- results[order(results$hr, results$NES),]
results$Description <- factor(results$Description,
                              levels = unique(results$Description))

## setup the plot
library(ggplot2)
source("ggplot2-themes.R")
plot <- ggplot(data = results[results$Category != "Other",],
               aes(x = as.factor(hr),
                   y = Description,
                   fill = NES)) +
    scale_fill_distiller(name = "NES ",
                         palette = "RdYlBu",
                         na.value = "#2C6CAD") +
    geom_tile(color = "black") +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    xlab("") + ylab("") +
    coord_fixed(ratio = 1) +
    #theme1 +
    theme(axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 32),
          strip.text.y = element_text(size = 20,
                                      angle = 360),
          strip.background = element_rect(fill = "grey90"),
          legend.position = "bottom",
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 24),
          panel.background = element_rect(fill = "white",
                                          color = "white"))

## open dvice for plotting
png(filename = "../img/gsea-reactome-heatmap.png", height = 1000, width = 1200)
print(plot)
dev.off()


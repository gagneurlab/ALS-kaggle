suppressPackageStartupMessages({
    library(OUTRIDER)
    library(umap)
    library(dplyr)
    library(ggplot2)
})

ods <- readRDS(snakemake@input$ods)
ods <- estimateSizeFactors(ods)

ct <- as.data.table(colData(ods))
ct[, case_control := sapply(sampleID, function(x) strsplit(x, '\\.')[[1]][1])]
table(ct$case_control)

center_log_sF_counts <- function(ODS){
    sizeFactors <- sizeFactors(ODS)
    counts <- counts(ODS, normalized = FALSE)
    sFcounts <- t(t(counts)/sizeFactors)
    logSFcounts <- log(sFcounts + 1)
    centerCounts <- logSFcounts - mean(logSFcounts)
}

centeredCounts <- center_log_sF_counts(ods)
um <- umap(t(centeredCounts), random_state = 5,  n_components = 3)

layout <- as.data.table(um$layout)
layout <- cbind(ct, layout)

ggplot(layout, aes(V1, V2, color = case_control)) +
    geom_point() + 
    theme_bw(base_size = 14) +
    theme(legend.position = 'bottom')
ggsave(snakemake@output$umap12)

ggplot(layout, aes(V1, V3, color = case_control)) + 
    geom_point() + 
    theme_bw(base_size = 14) +
    theme(legend.position = 'bottom')
ggsave(snakemake@output$umap13)

ggplot(layout, aes(V2, V3, color = case_control)) +
    geom_point() + 
    theme_bw(base_size = 14) +
    theme(legend.position = 'bottom')
ggsave(snakemake@output$umap23)

suppressPackageStartupMessages({
    library(data.table)
    library(OUTRIDER)
    library(umap)
    library(dplyr)
    library(ggplot2)
})

anno_col <- snakemake@wildcards$annotation_column
sample_annotation <- fread(snakemake@input$sample_annotation)
sample_annotation[, Participant_ID := make.names(sample_annotation$Participant_ID)]

ods <- readRDS(snakemake@input$ods)
ods <- estimateSizeFactors(ods)

ct <- as.data.table(colData(ods))
ct <- merge(ct, sample_annotation, by.x = 'sampleID', by.y = 'Participant_ID')
ct[, annotation_column := get(anno_col)]

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

ggplot(layout, aes(V1, V2, color = annotation_column)) +
    geom_point() + 
    labs(color = anno_col) + 
    theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom')
ggsave(snakemake@output$umap12)

ggplot(layout, aes(V1, V3, color = annotation_column)) + 
    geom_point() + 
    labs(color = anno_col) + 
    theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom')
ggsave(snakemake@output$umap13)

ggplot(layout, aes(V2, V3, color = annotation_column)) +
    geom_point() + 
    labs(color = anno_col) + 
    theme_bw(base_size = 14) + 
    theme(legend.position = 'bottom')
ggsave(snakemake@output$umap23)

suppressPackageStartupMessages({
    library(data.table)
    library(OUTRIDER)
    library(ggplot2)
})

anno_groups <- snakemake@params$anno_groups
ods <- readRDS(snakemake@input$ods)
sample_annotation <- fread(snakemake@input$sample_annotation)
sample_annotation <- sample_annotation[, c('Participant_ID', anno_groups), with = FALSE]
sample_annotation[, Participant_ID := make.names(Participant_ID)]

col_data <- colData(ods)
col_data <- merge(
    col_data,
    sample_annotation,
    by.x = 'sampleID',
    by.y = 'Participant_ID',
    all.x = TRUE
)
colData(ods) <- col_data
colnames(ods) <- col_data$sampleID

png(snakemake@output$power)
plotPowerAnalysis(ods)
dev.off()

png(snakemake@output$heatmap_raw, width = 1000, height = 1000)
plotCountCorHeatmap(ods, normalized = FALSE, colGroups = anno_groups, colColSet = "Set1")
dev.off()

png(snakemake@output$heatmap_fit, width = 1000, height = 1000)
plotCountCorHeatmap(ods, normalized = TRUE, colGroups = anno_groups, colColSet = "Set1")
dev.off()

png(snakemake@output$gene_vs_sample, width = 1000, height = 1000)
plotCountGeneSampleHeatmap(ods, normalized = TRUE, nGenes = 1000, colGroups = anno_groups, colColSet = "Set1")
dev.off()

png(snakemake@output$per_sample)
plotAberrantPerSample(ods)
dev.off()

png(snakemake@output$qqplot)
plotQQ(ods, global = TRUE)
dev.off()


plotEncDimSearch(ods) + theme_bw()
ggsave(snakemake@output$dim_search)

plotExpressedGenes(ods)
ggsave(snakemake@output$expr_genes)

res <- results(ods, all = TRUE)
ggplot(res, aes(pValue)) + geom_histogram(bins = 100) + theme_classic()
ggsave(snakemake@output$pvalues)

cnts <- data.table(
    observed = as.vector(counts(ods)),
    expected = as.vector(normalizationFactors(ods)),
    aberrant = as.vector(aberrant(ods))
)
ggplot(cnts, aes(log10(expected + 1), log10(observed + 1), color = aberrant)) +
    geom_point() +
    scale_color_brewer(palette = 'Set1') +
    theme_classic() +
    theme(legend.position = 'bottom')
ggsave(snakemake@output$exp_obs_counts)

# cnts_long <- melt(
#     cnts,
#     measure.vars = c('expected', 'observed'),
#     variable.name = 'count_type',
#     value.name = 'counts'
# )
# 
# ggplot(cnts_long, aes(log10(counts + 1), after_stat(count), color = count_type)) +
#     geom_density(position = 'identity') +
#     scale_color_brewer(palette = 'Set1') +
#     theme_classic() +
#     theme(legend.position = 'bottom')
# 

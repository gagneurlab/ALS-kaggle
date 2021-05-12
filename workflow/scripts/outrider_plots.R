suppressPackageStartupMessages({
    library(data.table)
    library(OUTRIDER)
    library(ggplot2)
})

anno_groups <- snakemake@params$anno_groups
ods <- readRDS(snakemake@input$ods)
sample_annotation <- fread(snakemake@input$sample_annotation, na.strings = c('NA', ''))
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


# BCV
estimateThetaWithoutAutoCorrect <- function(ods){
    ods1 <- OutriderDataSet(countData=counts(ods), colData=colData(ods))
    # use rowMeans as expected means
    normalizationFactors(ods1) <- matrix(
        rowMeans(counts(ods1)),
        ncol=ncol(ods1),
        nrow=nrow(ods1)
    )
    ods1 <- fit(ods1)
    theta(ods1)
}

before <- data.table(
    autocorrect = "No",
    BCV = 1/sqrt(estimateThetaWithoutAutoCorrect(ods))
)
after <- data.table(
    autocorrect = "Yes",
    BCV = 1/sqrt(theta(ods))
)
bcv_dt <- cbind(geneID = rownames(ods), rbind(before, after))

ggplot(bcv_dt, aes(autocorrect, BCV)) +
    geom_boxplot() +
    theme_bw(base_size = 14) +
    labs(
        x = "Autoencoder corrected",
        y = "BCV",
        title = "Biological coefficient of variation"
    )
ggsave(snakemake@output$bcv)

fwrite(
    bcv_dt[autocorrect == 'Yes'][order(BCV, decreasing = TRUE)],
    snakemake@output$bcv_dt
)


# Heatmaps
dev.off()

png(snakemake@output$gene_vs_sample, width = 1000, height = 1000)
plotCountGeneSampleHeatmap(
    ods,
    normalized = FALSE,
    nGenes = 50,
    colGroups = anno_groups,
    main = 'Gene counts vs Sample Heatmap'
)
dev.off()

png(snakemake@output$heatmap_raw, width = 1000, height = 1000)
plotCountCorHeatmap(
    ods,
    normalized = FALSE,
    colGroups = anno_groups,
    main = 'Raw Count Correlation'
)
dev.off()

png(snakemake@output$heatmap_fit, width = 1000, height = 1000)
plotCountCorHeatmap(
    ods,
    normalized = TRUE,
    colGroups = anno_groups,
    main = 'Normalized Count Correlation'
)
dev.off()

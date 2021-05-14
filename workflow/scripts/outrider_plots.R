suppressPackageStartupMessages({
    library(data.table)
    library(OUTRIDER)
    library(ggplot2)
    library(viridis)
    library(patchwork)
})

ods <- readRDS(snakemake@input$ods)

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

obs_cnt_mtx <- log10(counts(ods, normalized = FALSE) + 1)
exp_cnt_mtx <- log10(normalizationFactors(ods) + 1)

cnts <- data.table(
    geneID=rep(rownames(ods), ncol(ods)), 
    sampleID=rep(colnames(ods), each=nrow(ods)),
    observed = as.vector(obs_cnt_mtx),
    mean_counts = rep(rowMeans(obs_cnt_mtx), ncol(ods)),
    expected = as.vector(exp_cnt_mtx),
    observed_centered = as.vector(obs_cnt_mtx) - rowMeans(obs_cnt_mtx),
    expected_centered = as.vector(exp_cnt_mtx) - rowMeans(exp_cnt_mtx),
    aberrant = as.vector(aberrant(ods))
)

r_squared <- c(
    exp_obs = round(summary(lm(observed ~ expected, data = cnts))$r.squared, 2),
    exp_obs_centered = round(summary(lm(observed_centered ~ expected_centered, data = cnts))$r.squared, 2),
    mean_obs = round(summary(lm(observed ~ mean_counts, data = cnts))$r.squared, 2)
)

p1 <- ggplot(cnts, aes(x=expected_centered, y=observed_centered)) +
    geom_hex(
        aes(fill=stat(log10(count))),
        binwidth = 0.05,
        show.legend = FALSE
    ) +
    scale_fill_gradientn(colours = viridis(30)) +
    geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
    annotate(
        'text',
        x = -1.5,
        y = 2.5,
        label = paste0("r^2 == ", r_squared['exp_obs_centered']),
        parse = TRUE,
        size = 4
    ) +
    labs(
        title = 'Expected vs observed gene counts',
        subtitle = 'Centered by gene',
        x = "Expected log10(counts+1)",
        y = "Observed log10(counts+1)"
    ) +
    theme_classic()

p2 <- ggplot(cnts, aes(x=expected, y=observed)) +
    geom_hex(
        aes(fill=stat(log10(count))),
        binwidth = 0.05,
        show.legend = FALSE
    ) +
    scale_fill_gradientn(colours = viridis(30)) +
    geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
    annotate(
        'text',
        x = 1,
        y = 6,
        label = paste0("r^2 == ", r_squared['exp_obs']),
        parse = TRUE,
        size = 4
    ) +
    labs(
        title = 'Expected vs observed gene counts',
        x = "Expected log10(counts+1)",
        y = "Observed log10(counts+1)"
    ) +
    theme_classic()

p3 <- ggplot(cnts, aes(x=mean_counts, y=observed)) +
    geom_hex(
        aes(fill=stat(log10(count))),
        binwidth = 0.05,
        show.legend = FALSE
    ) +
    scale_fill_gradientn(colours = viridis(30)) +
    annotate(
        'text',
        x = 1,
        y = 6,
        label = paste0("r^2 == ", r_squared['mean_obs']),
        parse = TRUE,
        size = 4
    ) +
    geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
    labs(
        title = 'Mean vs observed gene counts',
        x = "Observed mean log10(counts+1)",
        y = "Observed log10(counts+1)"
    ) +
    theme_classic()

p1 | p2 | p3

ggsave(snakemake@output$exp_obs_counts, width = 12, height = 5)

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

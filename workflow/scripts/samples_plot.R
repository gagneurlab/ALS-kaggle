# Creates a plot containing the expression fold change of candidate samples

suppressPackageStartupMessages({
  library(data.table)
  library(OUTRIDER)
  library(ggplot2)
  library(ggbeeswarm)
  library(cowplot)
})


ods <- readRDS(snakemake@input$ods)

# Specify the candidate samples
candidate_samples <- data.table(
  sample = c('CASE.NEUEK191WYC', 'CASE.NEUBK117YXL', 'CASE.NEUZT557DHF', 'CASE.NEUVX902YNL',
             'CASE.NEULD354RZB', 'CASE.NEUTA689LN5', 'CASE.NEUGW326BRV', 'CASE.NEUME498PCJ', 
             'CASE.NEURR881FKY'),
  gene_name = c('NEK1', 'OPTN', 'OPTN', 'SPG11', 'NOP56', 'TUFM', 'CLPB', 'PIKFYVE', 'ERCC2')
)

# add gene_name
gene_table <- fread(snakemake@config$genes_gtf)
candidate_samples <- merge(candidate_samples, gene_table[, .(gene_name, geneID)])
candidate_samples[, aux := paste(sample, geneID, sep = '-')]

ens_ids <- candidate_samples$geneID %>% unique()

# Obtain the expression fold changes from the ods object
fc <- 2^assays(ods[ens_ids,])$l2fc
fc_dt <- reshape2::melt(fc) %>% as.data.table()
colnames(fc_dt) <- c('geneID', 'sample', 'FC')
fc_dt[, aux := paste(sample, geneID, sep = '-')]
fc_dt[, causal := aux %in% candidate_samples[, aux]]
fc_dt <- merge(fc_dt, gene_table[, .(gene_name, geneID)])

# Plot the fold changes
gf <- ggplot(fc_dt, aes(reorder(gene_name, FC, min), FC)) + 
  geom_beeswarm(size = .15, color = 'gray60', cex = .7) +
  geom_point(data = fc_dt[causal == T], color = 'firebrick') +
  scale_y_log10() + labs(x = '', y = 'Expression fold change') + 
  theme_cowplot() + background_grid() + coord_flip()
save_plot(snakemake@output$samples_plot, gf, base_height = 9, base_width = 8)


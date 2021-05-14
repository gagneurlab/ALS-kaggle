suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

sample_annotation <- fread('data/processed/sample_annotation.tsv')
setnames(sample_annotation, c('ALSFRS-R_Baseline', 'ALSFRS-R_Latest'), c('Baseline', "Latest"))

alsfrs_r <- melt(
    sample_annotation,
    measure.vars = c("Baseline", "Latest"),
    variable.name = 'timepoint',
    value.name = 'ALSFRS_R'
)
setnames(alsfrs_r, 'ALSFRS-R_Progression_Slope', 'progression_slope')
alsfrs_r[, timepoint_num := as.numeric(timepoint)]
alsfrs_r[, progressing := progression_slope < 0]

ggplot(alsfrs_r) +
    geom_histogram(
        aes(x = ALSFRS_R, y = ..density.., color = timepoint, fill = timepoint),
        position = 'identity',
        bins = 25,
        alpha = 0.3
    ) +
    geom_density(
        aes(x = ALSFRS_R, color = timepoint),
        position = 'identity',
        size = 1
    ) +
    scale_fill_brewer(palette = 'Set1') +
    scale_color_brewer(palette = 'Set1') +
    theme_classic() +
    theme(legend.position = 'bottom')


alsfrs_r[, age_group := cut(Age_At_Symptom_Onset, 3)]

ggplot(remove_missing(alsfrs_r, na.rm = TRUE, vars = 'progression_slope'), 
       aes(timepoint, ALSFRS_R, group = Participant_ID, col = progressing),
    ) +
    geom_line() +
    geom_point() +
    facet_wrap(~age_group) +
    scale_color_brewer(palette = 'Dark2') +
    labs(title = 'Progression', x = '', y = 'ALSFRS-R') +
    theme_classic() +
    theme(legend.position = 'bottom')

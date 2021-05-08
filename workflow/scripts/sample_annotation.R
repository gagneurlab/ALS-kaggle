suppressPackageStartupMessages({
    library(data.table)
})

DIR <- snakemake@config$metadata$dir
#DIR <- '/s/raw/als/kaggle/end-als/clinical-data/filtered-metadata/metadata'

# Main sample annotations
sa_participants <- fread(snakemake@input$participants)
sa_dataportal <- fread(snakemake@input$dataportal)

sa_als <- merge(sa_participants, sa_dataportal)
colnames(sa_als) <- gsub(' ', '_', colnames(sa_als))

# Relevant genetic/clinical information
gene_mut <- fread(file.path(DIR, 'clinical/ALS_Gene_Mutations.csv'))
gene_mut <- gene_mut[, .(Participant_ID, ang, c9orf72, fus, setx, sod1, tau, tdp43, vapb, vcp)]
colnames(gene_mut)[2:length(colnames(gene_mut))] <- toupper(colnames(gene_mut)[2:length(colnames(gene_mut))])

# med_log <- fread('/s/raw/als/kaggle/end-als/clinical-data/filtered-metadata/metadata/clinical/ANSWER_ALS_Medications_Log.csv')
# fam_hist <- fread('/s/raw/als/kaggle/end-als/clinical-data/filtered-metadata/metadata/clinical/Family_History_Log.csv')  # column fhpsysp, but very few

med_hist <- fread(file.path(DIR, 'clinical/Medical_History.csv'))
med_hist_abbrev <- med_hist[, .(med_hist = paste(medhxdsc, collapse = ", ")), by = Participant_ID]

vit_signs <- fread(file.path(DIR, 'clinical/Vital_Signs.csv'))
vit_signs <- vit_signs[, .(Participant_ID, height, weight, bmi, hr, bpsys)][, .SD[1], by = Participant_ID]

mt <- merge(vit_signs, gene_mut, all = TRUE)
mt <- merge(mt, med_hist_abbrev, all = TRUE)

sa_als <- merge(sa_als, mt, by = 'Participant_ID', all = TRUE)
fwrite(sa_als, snakemake@output$merged, sep = '\t')

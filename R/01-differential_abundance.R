# -----------------------------------------------------------------------------
# Script 01: Human Plasma Differential Abundance Analysis
# Author: Rani Powers
# Last updated: July 31, 2018
#
# Reads in the human metabolomics eset and performs differential abundance 
# analysis for Nexus cohort alone, HTP cohort alone, and both cohorts combined
# (linear model fit for batch, age, and sex from 00-model_fitting.R).
#
# Saves volcano plots in the Results/Figure_1/ folder.
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')

# Read in the raw data (no model fitting)
raw_eset = readRDS('Data/Human_plasma_metabolomics_filtered.rds')
met_data = as.data.frame(exprs(raw_eset))
sample_data = pData(raw_eset)
feature_data = fData(raw_eset)

# ------------------------------- NEXUS ONLY ----------------------------------
NEXUS = row.names(sample_data[sample_data$Batch_Name == 'Nexus',])
nexus_data = met_data[,NEXUS]
nexus_sample_data = sample_data[NEXUS,]

# Nexus differential abundance - no fitting
nexus_design = model.matrix(~0 + as.factor(Karyotype),
                      data = nexus_sample_data)
colnames(nexus_design) = c('D21', 'T21')
nexus_contrast = makeContrasts(T21vsD21 = T21-D21, levels = nexus_design)
nexus_fit = eBayes(contrasts.fit(lmFit(nexus_data, nexus_design), nexus_contrast))
nexus_sigtable = topTable(nexus_fit, adjust = 'fdr', number = nrow(nexus_data)+1)
nexus_sigtable$Compound_ID = row.names(nexus_sigtable)
nexus_sigtable = join(nexus_sigtable, feature_data)

# Nexus differential abundance - with age and sex fit
nexus_design2 = model.matrix(~0 + Karyotype+Age+Sex, 
                             nexus_sample_data[,c('Karyotype', 'Age', 'Sex')])
colnames(nexus_design2)[1:2] = c('D21', 'T21')
nexus_contrast2 = makeContrasts(T21vsD21 = T21-D21, levels = nexus_design2)
nexus_fit2 = eBayes(contrasts.fit(lmFit(nexus_data, nexus_design2), nexus_contrast2))
nexus_sigtable2 = topTable(nexus_fit2, adjust = 'fdr', number = nrow(nexus_data)+1)
nexus_sigtable2$Compound_ID = row.names(nexus_sigtable2)
nexus_sigtable2 = join(nexus_sigtable2, feature_data)

# Plot volcano with no fitting
pdf('Results/Figure_1/Fig1_Supp1B_nofit.pdf',
    height = 7, width = 6)
plot_volcano(x = nexus_sigtable$logFC, y = -log10(nexus_sigtable$adj.P.Val), plot_fc_lines = F,
             x_cutoffs = c(-.01, .01), xlim = c(-3,3),
             main = 'Nexus\nno fitting', pt_labels = nexus_sigtable$Compound_Name)
dev.off()

# Plot volcano with fitting
pdf('Results/Figure_1/Fig1_Supp1B_fit.pdf',
    height = 7, width = 6)
plot_volcano(x = nexus_sigtable2$logFC, y = -log10(nexus_sigtable2$adj.P.Val), plot_fc_lines = F,
             x_cutoffs = c(-.01, .01), xlim = c(-3,3),
             main = 'Nexus\nmodel fit with batch, age and sex', pt_labels = nexus_sigtable2$Compound_Name)
dev.off()

# ------------------------------- HTP ONLY ------------------------------------ 
HTP = row.names(sample_data[sample_data$Batch_Name == 'HTP',])
htp_data = met_data[,HTP]
htp_sample_data = sample_data[HTP,]

# HTP differential abundance
htp_design = model.matrix(~0 + as.factor(Karyotype),
                            data = htp_sample_data)
colnames(htp_design) = c('D21', 'T21')
htp_contrast = makeContrasts(T21vsD21 = T21-D21, levels = htp_design)
htp_fit = eBayes(contrasts.fit(lmFit(htp_data, htp_design), htp_contrast))
htp_sigtable = topTable(htp_fit, adjust = 'fdr', number = nrow(htp_data)+1)
htp_sigtable$Compound_ID = row.names(htp_sigtable)
htp_sigtable = join(htp_sigtable, feature_data)

# HTP differential abundance - with age and sex fit
htp_design2 = model.matrix(~0 + Karyotype+Age+Sex, 
                           htp_sample_data[,c('Karyotype', 'Age', 'Sex')])
colnames(htp_design2)[1:2] = c('D21', 'T21')
htp_contrast2 = makeContrasts(T21vsD21 = T21-D21, levels = htp_design2)
htp_fit2 = eBayes(contrasts.fit(lmFit(htp_data, htp_design2), htp_contrast2))
htp_sigtable2 = topTable(htp_fit2, adjust = 'fdr', number = nrow(htp_data)+1)
htp_sigtable2$Compound_ID = row.names(htp_sigtable2)
htp_sigtable2 = join(htp_sigtable2, feature_data)

# Plot volcano with no fitting
pdf('Results/Figure_1/Fig1_Supp1C_nofit.pdf',
    height = 7, width = 6)
plot_volcano(x = htp_sigtable$logFC, y = -log10(htp_sigtable$adj.P.Val), plot_fc_lines = F,
             x_cutoffs = c(-.01, .01), xlim = c(-3,3),
             main = 'HTP\nno fitting', pt_labels = htp_sigtable$Compound_Name)
dev.off()

# Plot volcano with fitting
pdf('Results/Figure_1/Fig1_Supp1C_fit.pdf',
    height = 7, width = 6)
plot_volcano(x = htp_sigtable2$logFC, y = -log10(htp_sigtable2$adj.P.Val), plot_fc_lines = F,
             x_cutoffs = c(-.01, .01), xlim = c(-3,3),
             main = 'HTP\nmodel fit with batch, age and sex', pt_labels = htp_sigtable2$Compound_Name)
dev.off()

# ------------------------------- NEXUS & HTP --------------------------------- 
# BOTH differential abundance
eset = readRDS('Data/Human_plasma_metabolomics_filtered.rds')
met_data = as.data.frame(exprs(eset))
sample_data = pData(eset)
feature_data = fData(eset)
both_design = model.matrix(~0 + Karyotype+Age+Sex+Batch, 
                           sample_data[,c('Karyotype', 'Age', 'Sex', 'Batch')])
colnames(both_design)[1:2] = c('D21', 'T21')
both_contrast = makeContrasts(T21vsD21 = T21-D21, levels = both_design)
both_fit = eBayes(contrasts.fit(lmFit(met_data, both_design), both_contrast))
both_sigtable = topTable(both_fit, adjust = 'fdr', number = nrow(met_data)+1)
both_sigtable$Compound_ID = row.names(both_sigtable)
both_sigtable = join(both_sigtable, feature_data)

# Plot volcano with fitting
pdf('Results/Figure_1/Fig1A.pdf',
    height = 7, width = 5)
plot_volcano(x = both_sigtable$logFC, y = -log10(both_sigtable$adj.P.Val), plot_fc_lines = F,
             x_cutoffs = c(-.01, .01), xlim = c(-2,2),
             main = 'Nexus & HTP\nmodel fit with batch, age and sex', pt_labels = both_sigtable$Compound_Name)
dev.off()

# Save table
write.table(both_sigtable, 'Results/Figure_1/Fig1_SuppTable.csv', 
            sep = ',', row.names = F)

# -----------------------------------------------------------------------------
# Script 06: Human Plasma Metabolites & MSD Cytokines
# Author: Rani Powers
# Last updated: August 20, 2018
#
# Correlates Cohort 2 human plasma metabolomics data with MSD cytokine measurements.
# Saves heatmaps and scatterplots in the Results/Figure_6/ folder.
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
source('R/heatmap.3.R')
if (!dir.exists('Results/Figure_6')) { dir.create('Results/Figure_6', recursive = T) }

# Get HTP metabolomics data
eset = readRDS('Data/Human_plasma_metabolomics_filtered.rds')
sample_data = pData(eset)
HTP = sample_data[sample_data$Batch_Name == 'HTP', 'Barcode']
sample_data = sample_data[HTP,]
met_data = exprs(eset)[,HTP]
feature_data = fData(eset)

# Match with MSD data
msd_data = read.csv('Data/Human_plasma_MSD.csv', stringsAsFactors = F)
samples = intersect(unique(msd_data$Sample), colnames(met_data))
met_data = met_data[,samples]
msd_data = msd_data[msd_data$Sample %in% samples,]
sample_data = sample_data[samples,]

msd_data = msd_data[msd_data$Assay != 'IL-8',]
msd_data = as.data.frame(spread(msd_data, key = 'Assay', value = 'Conc_Mean'))
row.names(msd_data) = msd_data$Sample
msd_data$Sample = NULL
msd_data = as.data.frame(t(msd_data[samples,]))

# Fit model on metabolomics data
design = model.matrix(~0+Age+as.factor(Sex), sample_data[,c('Batch', 'Age', 'Sex')])
fit_met = lmFit(met_data, design = design)
met_data = residuals(fit_met, met_data)

# Fit model on MSD data
fit_msd = lmFit(msd_data, design = design)
msd_data = residuals(fit_msd, msd_data)

# Correlate all metabolites and cytokines
ALL_DATA = rbind(met_data[,samples], msd_data[,samples])
data_for_correlation = as.data.frame(t(ALL_DATA))
all_corr_s = cor(data_for_correlation, method = 'spearman', use = 'pairwise.complete.obs')
t21 = sample_data[sample_data$Karyotype == 'T21', 'Barcode']
d21 = sample_data[sample_data$Karyotype == 'D21', 'Barcode']
t21_corr_s = cor(data_for_correlation[t21,], 
                 method = 'spearman', use = 'pairwise.complete.obs')
d21_corr_s = cor(data_for_correlation[d21,], 
                 method = 'spearman', use = 'pairwise.complete.obs')

# Plot just kynurenine vs all cytokines, ranked
kyn_cor_t21 = sort(t21_corr_s['C00328', 92:ncol(t21_corr_s)], 
                   decreasing = T)
pdf('Results/Figure_6/Fig6_kyn_heatmap.pdf',
    height = 10, width = 6)
plot_single_row_heatmap(data_vector = kyn_cor_t21, karyotype = 'T21', 
                        column_label = 'kynurenine', n = length(t21))
dev.off()

# Plot top 2 scatterplots with kynurenine
top2_kyn = names(kyn_cor_t21)[1:2]
data_for_correlation = as.data.frame(data_for_correlation)
for (cyto in top2_kyn){
  pdf(paste0('Results/Figure_6/Fig6_kyn_', cyto, '.pdf'), 
      height = 6, width = 6)
  par(mar = c(8,4,3,2))
  plot_msd_met_scatterplot(met1 = 'C00328', met2 = cyto, plot_fit = F, 
                           df = data_for_correlation, t21 = t21, d21 = d21)
  dev.off()
}

# Plot just quinolinic acid vs all cytokines, ranked
quin_cor_t21 = sort(t21_corr_s['C03722', 92:ncol(t21_corr_s)], 
                   decreasing = T)
pdf('Results/Figure_6/Fig6_quin_heatmap.pdf',
    height = 10, width = 6)
plot_single_row_heatmap(data_vector = quin_cor_t21, karyotype = 'T21', 
                        column_label = 'quinolinic acid', n = length(t21))
dev.off()

# Plot top 2 scatterplots with quinolinic acid
top2_quin = names(quin_cor_t21)[1:2]
for (cyto in top2_quin){
  pdf(paste0('Results/Figure_6/Fig6_quin_', cyto, '.pdf'), 
      height = 6, width = 6)
  par(mar = c(8,4,3,2))
  plot_msd_met_scatterplot(met1 = 'C03722', met2 = cyto, plot_fit = F, 
                           df = data_for_correlation, t21 = t21, d21 = d21)
  dev.off()
}

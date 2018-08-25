# -----------------------------------------------------------------------------
# Script 00: Human Plasma Data Cleaning & Linear Model Fitting
# Author: Rani Powers
# Last updated: July 31, 2018
#
# Reads in the human metabolomics eset, filters low abundance metabolites and 
# metabolites in < 90% of the T21 and/or D21 samples. 
#
# Identifies and filters two outlier samples.
#
# Saves density and PCA plots from pre- and post-adjustment (batch, age, sex)
# in the Results/Figure_1/ folder.
# 
# Fits linear model to data with batch, age and sex as covariates. Saves the
# model residuals in a new eset (Human_plasma_residuals.rds).
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
if (!dir.exists('Results/Figure_1')) { dir.create('Results/Figure_1', recursive = T) }

# Read in human metabolomics data from Nexus and HTP cohorts
eset = readRDS('Data/Human_plasma_metabolomics.rds')
met_data = as.data.frame(exprs(eset))
sample_data = pData(eset)
feature_data = fData(eset)

# Filter low abundance metabolites and metabolites in < 90% of T21 or D21
met_data[met_data < 1000] = NA
met_t21 = sample_data[sample_data$Karyotype == 'T21', 'Barcode']
met_d21 = sample_data[sample_data$Karyotype == 'D21', 'Barcode']
in_t21 = apply(met_data[,met_t21], 1, function(x) sum(!is.na(x))) >= length(met_t21) * .9
in_d21 = apply(met_data[,met_d21], 1, function(x) sum(!is.na(x))) >= length(met_d21) * .9
to_keep = in_t21 | in_d21
feature_data[names(to_keep[to_keep == F]), 'Compound_Name']
to_keep = names(to_keep[to_keep == T])
met_data = met_data[to_keep,]
# 91 metabolites meet the criteria

# Use log2 data
met_data_log2 = log2(met_data[to_keep,])

# Plot density plot using the log2 data (no model fitting)
outliers = c('222762382', '222763227')  # We know these are > 2.5 sd from mean - see line 116
d = list()
pdf('Results/Figure_1/Fig1_Supp2A_withOutliers.pdf',
    height = 5, width = 5)
plot(density(na.omit(as.numeric(met_data_log2[,1]))), 
     col = ifelse(sample_data$Batch[1] == 1, cols[11], cols[8]),
     main = 'Density of log2 data\nno model fitting',
     ylim = c(0, .2), las = 1, xlab = 'log2 intensity (91 metabolites)')
for (i in 2:ncol(met_data_log2)){
  lines(density(na.omit(as.numeric(met_data_log2[,i]))),
        col = ifelse(sample_data$Barcode[i] %in% outliers, 'red',
                     ifelse(sample_data$Batch[i] == 1, cols[11], cols[8])))
  d = c(d, list(density(na.omit(as.numeric(met_data_log2[,i])))))
}
legend('topright', 
       lty = 1, col = c(cols[11], cols[8]),
       legend = c('Nexus', 'HTP'))
dev.off()

# Plot PCA using the log2 data (no model fitting)
met_data_log2_NAexcluded = na.omit(met_data_log2)
log2_pca = summary(prcomp(t(met_data_log2_NAexcluded), center = T, scale. = T))
pdf('Results/Figure_1/Fig1_Supp2C_withOutliers.pdf',
    height = 5, width = 5)
plot_batch_pca(log2_pca, bad = outliers,
               main = 'PCA of log2 data\nno model fitting',
               legend_side = 'topleft')
dev.off()

# Fit linear model with batch, age and sex as covariates
adjusted_met_data = removeBatchEffect(met_data_log2, 
                            batch = sample_data$Batch, 
                            design = model.matrix(~0+Age+Sex, data = sample_data))

# Plot density plot using adjusted data (after model fitting)
pdf('Results/Figure_1/Fig1_Supp2B_withOutliers.pdf',
    height = 5, width = 5)
plot(density(adjusted_met_data[,1]), 
     col = ifelse(sample_data$Batch[1] == 1, cols[11], cols[8]),
     main = 'Density of log2 data after model fitting\n(covariates = batch, age, sex)',
     ylim = c(0, .2), las = 1)
for (i in sample(2:ncol(adjusted_met_data))){
  lines(density(na.omit(adjusted_met_data[,i])),
        col = ifelse(sample_data$Barcode[i] %in% outliers, 'red',
                     ifelse(sample_data$Batch[i] == 1, cols[11], cols[8])))
}
legend('topleft', 
       lty = 1, col = c(cols[11], cols[8]),
       legend = c('Nexus', 'HTP'))
dev.off()

# Plot PCA using adjusted data (after model fitting)
adjusted_NAexcluded = na.omit(adjusted_met_data)
adjusted_pca = summary(prcomp(t(adjusted_NAexcluded), center = T, scale. = T))
# Draw lines showing > 3 SD from the mean on PCA
pca_data = as.data.frame(adjusted_pca$x)
pc1 = pca_data$PC1
pc2 = pca_data$PC2
pdf('Results/Figure_1/Fig1_Supp2D_withOutliers.pdf',
    height = 5, width = 5)
plot_batch_pca(adjusted_pca, bad = outliers,
               main = 'PCA of log2 data after model fitting\n(covariates = batch, age, sex)',
               legend_side = 'bottomright')
abline(v = sd(pc1)*3, col = 'grey')
abline(v = sd(pc1)*-3, col = 'grey')
abline(h = sd(pc2)*3, col = 'grey')
abline(h = sd(pc2)*-3, col = 'grey')
dev.off()

# Remove outliers
sd_cutoff = 2.5
outside_of_3sd_pc1 = row.names(pca_data[pc1 > sd(pc1) * sd_cutoff | pc1 < sd(pc1) * -1*sd_cutoff,])
outside_of_3sd_pc2 = row.names(pca_data[pc2 > sd(pc2) * sd_cutoff | pc2 < sd(pc2) * -1*sd_cutoff,])
outliers = intersect(outside_of_3sd_pc1, outside_of_3sd_pc2)
met_data_log2 = met_data_log2[,!names(met_data_log2) %in% outliers]
sample_data = sample_data[!row.names(sample_data) %in% outliers,]

# Plot density plot using the log2 data (no model fitting)
d = list()
pdf('Results/Figure_1/Fig1_Supp2A.pdf',
    height = 5, width = 5)
plot(density(na.omit(as.numeric(met_data_log2[,1]))), 
     col = ifelse(sample_data$Batch[1] == 1, cols[11], cols[8]),
     main = 'Density of log2 data\nno model fitting',
     ylim = c(0, .2), las = 1, xlab = 'log2 intensity (91 metabolites)')
for (i in 2:ncol(met_data_log2)){
  lines(density(na.omit(as.numeric(met_data_log2[,i]))),
        col = ifelse(sample_data$Barcode[i] %in% outliers, 'red',
                     ifelse(sample_data$Batch[i] == 1, cols[11], cols[8])))
  d = c(d, list(density(na.omit(as.numeric(met_data_log2[,i])))))
}
legend('topright', 
       lty = 1, col = c(cols[11], cols[8]),
       legend = c('Nexus', 'HTP'))
dev.off()

# Plot PCA using the log2 data (no model fitting)
met_data_log2_NAexcluded = na.omit(met_data_log2)
log2_pca = summary(prcomp(t(met_data_log2_NAexcluded), center = T, scale. = T))
pdf('Results/Figure_1/Fig1_Supp2C.pdf',
    height = 5, width = 5)
plot_batch_pca(log2_pca, bad = outliers,
               main = 'PCA of log2 data\nno model fitting',
               legend_side = 'topleft')
dev.off()

# Fit linear model with batch, age and sex as covariates
adjusted_met_data = removeBatchEffect(met_data_log2, 
                                      batch = sample_data$Batch, 
                                      design = model.matrix(~0+Age+Sex, data = sample_data))

# Plot density plot using adjusted data (after model fitting)
pdf('Results/Figure_1/Fig1_Supp2B.pdf',
    height = 5, width = 5)
plot(density(adjusted_met_data[,1]), 
     col = ifelse(sample_data$Batch[1] == 1, cols[11], cols[8]),
     main = 'Density of log2 data after model fitting\n(covariates = batch, age, sex)',
     ylim = c(0, .2), las = 1)
for (i in sample(2:ncol(adjusted_met_data))){
  lines(density(na.omit(adjusted_met_data[,i])),
        col = ifelse(sample_data$Barcode[i] %in% outliers, 'red',
                     ifelse(sample_data$Batch[i] == 1, cols[11], cols[8])))
}
legend('topleft', 
       lty = 1, col = c(cols[11], cols[8]),
       legend = c('Nexus', 'HTP'))
dev.off()

# Plot PCA using adjusted data (after model fitting)
adjusted_NAexcluded = na.omit(adjusted_met_data)
adjusted_pca = summary(prcomp(t(adjusted_NAexcluded), center = T, scale. = T))
pdf('Results/Figure_1/Fig1_Supp2D.pdf',
    height = 5, width = 5)
plot_batch_pca(adjusted_pca, bad = outliers,
               main = 'PCA of log2 data after model fitting\n(covariates = batch, age, sex)',
               legend_side = 'topleft')
dev.off()

# Fit linear model and save residuals
design_covars = model.matrix(~0+as.factor(Batch)+Age+as.factor(Sex), 
                             sample_data[,c('Batch', 'Age', 'Sex')])
fit = lmFit(met_data_log2, design = design_covars)
model_residuals = residuals(fit, met_data_log2)
saveRDS(model_residuals, 'Data/Human_plasma_residuals.rds')

# Save final, filtered metabolomics data
EXPRS_DATA = met_data_log2
PHENO_DATA = AnnotatedDataFrame(data = sample_data[names(EXPRS_DATA),],
                                dimLabels = c('Barcode', 'Info'))
FEATURE_DATA = AnnotatedDataFrame(feature_data[row.names(EXPRS_DATA),])
eset = new("ExpressionSet", 
            exprs = EXPRS_DATA, 
            phenoData = PHENO_DATA, 
            featureData = FEATURE_DATA)
saveRDS(eset, 'Data/Human_plasma_metabolomics_filtered.rds')

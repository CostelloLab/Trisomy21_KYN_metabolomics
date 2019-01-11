# -----------------------------------------------------------------------------
# Script 00: Human Plasma Data Cleaning & Linear Model Fitting
# Author: Rani Powers
# Last updated: December 11, 2018
#
# Reads in Data/human_plasma_metabolomics_eset.rds, filters low abundance 
# metabolites and metabolites in < 90% of the T21 and/or D21 samples. 
#
# Outputs:
# 1) Saves Data/human_plasma_metabolomics_eset_filtered.rds - eset with the raw
# data after the filtering steps described above.
#
# 2) Saves density plots from pre- and post-adjustment (batch, age, sex) in the
# Results/Supplementary_Figures/ folder.
# 
# 3) Fits linear model to data with batch, age and gender as covariates. Saves 
# the model residuals in a new eset (Data/human_plasma_residuals.rds).
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
if (!dir.exists('Results')){
  dir.create('Results/Supplementary_Figures/', recursive = T)
}

# Read in human metabolomics data from Nexus and HTP cohorts
eset = readRDS('Data/human_plasma_metabolomics_eset.rds')
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
d = list()
pdf('Results/Supplementary_Figures/Supp_Fig_1b_density_plot_before_fitting.pdf',
    height = 5, width = 5)
plot(density(na.omit(as.numeric(met_data_log2[,1]))), 
     col = ifelse(sample_data$Batch[1] == 1, cols[11], cols[8]),
     main = 'Density of log2 data\nno model fitting',
     ylim = c(0, .2), las = 1, xlab = 'log2 intensity (91 metabolites)')
for (i in 2:ncol(met_data_log2)){
  lines(density(na.omit(as.numeric(met_data_log2[,i]))),
        col = ifelse(sample_data$Batch[i] == 1, cols[11], cols[8]))
  d = c(d, list(density(na.omit(as.numeric(met_data_log2[,i])))))
}
legend('topright', 
       lty = 1, col = c(cols[11], cols[8]),
       legend = c('Cohort 1', 'Cohort 2'))
dev.off()

# Fit linear model with batch, age and sex as covariates
# (Use removeBatchEffect to get adjusted intensity values for plotting)
adjusted_met_data = removeBatchEffect(met_data_log2, 
                                      batch = sample_data$Batch, 
                                      design = model.matrix(~0+Age+Sex, data = sample_data))

# Plot density plot using adjusted data (after model fitting)
pdf('Results/Supplementary_Figures/Supp_Fig_1c_density_plot_after_fitting.pdf',
    height = 5, width = 5)
plot(density(adjusted_met_data[,1]), 
     col = ifelse(sample_data$Batch[1] == 1, cols[11], cols[8]),
     main = 'Density of log2 data after model fitting\n(covariates = batch, age, sex)',
     ylim = c(0, .2), las = 1)
for (i in sample(2:ncol(adjusted_met_data))){
  lines(density(na.omit(adjusted_met_data[,i])),
        col = ifelse(sample_data$Batch[i] == 1, cols[11], cols[8]))
}
legend('topleft', 
       lty = 1, col = c(cols[11], cols[8]),
       legend = c('Cohort 1', 'Cohort 2'))
dev.off()

# Fit linear model and save residuals
design_covars = model.matrix(~0+as.factor(Batch)+Age+as.factor(Sex), 
                             sample_data[,c('Batch', 'Age', 'Sex')])
fit = lmFit(met_data_log2, design = design_covars)
model_residuals = residuals(fit, met_data_log2)
saveRDS(model_residuals, 'Data/human_plasma_residuals.rds')

# Save final, filtered metabolomics data
EXPRS_DATA = met_data_log2
PHENO_DATA = AnnotatedDataFrame(data = sample_data[names(EXPRS_DATA),],
                                dimLabels = c('Barcode', 'Info'))
FEATURE_DATA = AnnotatedDataFrame(feature_data[row.names(EXPRS_DATA),])
eset = new("ExpressionSet", 
           exprs = EXPRS_DATA, 
           phenoData = PHENO_DATA, 
           featureData = FEATURE_DATA)
saveRDS(eset, 'Data/human_plasma_metabolomics_eset_cleaned.rds')
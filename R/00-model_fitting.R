# -----------------------------------------------------------------------------
# Script 00: Human Plasma Data Cleaning & Linear Model Fitting
# Author: Rani Powers
#
# Reads in the human metabolomics eset, filters low abundance metabolites and 
# metabolites in < 90% of the T21 and/or D21 samples. Performs age-, sex- and
# batch-adjustment using linear model.
#
# Outputs:
# 1) Saves Data/Human_plasma_metabolomics_eset_filtered.rds - eset with the raw
# data after the filtering steps described above.
#
# 2) Saves Data/Human_plasma_metabolomics_eset_filtered_adjusted.rds - eset with 
# the adjusted data after age/sex/batch adjustment.
#
# 3) Saves supplementary density plots from pre- and post-adjustment (batch, age, sex)
# in the Results/Supplementary_Figures/ folder.
# 
# 4) Saves Results/Supplementary_Files/Supplementary_Data_1.csv with cohort 
# details and Results/Supplementary_Files/Supplementary_Data_2.csv with raw 
# metabolomics data
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')

# Read in human metabolomics data from Nexus and HTP cohorts
eset = readRDS('Data/Human_plasma_metabolomics_eset.rds')
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
to_keep = names(to_keep[to_keep == T])
met_data = met_data[to_keep,]
feature_data = feature_data[to_keep,]
# 91 metabolites meet the criteria

# Use log2 data
met_data_log2 = log2(met_data[to_keep,])

## PLOT 
# Density plot using the log2 data (no model fitting)
d = list()
pdf('Results/Supplementary_Figures/SuppFig_1B_densityPlot_beforeFitting.pdf',
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

## PLOT
# Density plot using adjusted data (after model fitting)
pdf('Results/Supplementary_Figures/SuppFig_1C_densityPlot_afterFitting.pdf',
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

# Save final, filtered metabolomics data
EXPRS_DATA = met_data_log2
PHENO_DATA = AnnotatedDataFrame(data = sample_data[names(EXPRS_DATA),],
                                dimLabels = c('Barcode', 'Info'))
FEATURE_DATA = AnnotatedDataFrame(feature_data[row.names(EXPRS_DATA),])
eset = new("ExpressionSet", 
           exprs = EXPRS_DATA, 
           phenoData = PHENO_DATA, 
           featureData = FEATURE_DATA)
saveRDS(eset, 'Data/Human_plasma_metabolomics_eset_filtered.rds')

# Save final, filtered metabolomics data with age/batch/sex adj data
EXPRS_DATA = as.data.frame(adjusted_met_data)
PHENO_DATA = AnnotatedDataFrame(data = sample_data[names(EXPRS_DATA),],
                                dimLabels = c('Barcode', 'Info'))
FEATURE_DATA = AnnotatedDataFrame(feature_data[row.names(EXPRS_DATA),])
eset = new("ExpressionSet", 
           exprs = EXPRS_DATA, 
           phenoData = PHENO_DATA, 
           featureData = FEATURE_DATA)
saveRDS(eset, 'Data/Human_plasma_metabolomics_eset_filtered_adjusted.rds')

# Save supplementary files
supp_cohort_data = sample_data[,c('Barcode', 'Sex', 'Age', 'Karyotype')]
names(supp_cohort_data)[1] = 'Sample_ID'
write.table(supp_cohort_data, 
            'Results/Supplementary_Files/Supplementary_Data_1AB.csv',
            sep = ',', row.names = F)

supp_met_data2 = cbind(feature_data[,1:5],
                       met_data[,sample_data$Batch == 1])
write.table(supp_met_data2,
            'Results/Supplementary_Files/Supplementary_Data_2.csv',
            sep = ',', row.names = F)

supp_met_data3 = cbind(feature_data[,c(1:2,6:8)],
                       met_data[,sample_data$Batch == 2])
write.table(supp_met_data3,
            'Results/Supplementary_Files/Supplementary_Data_3.csv',
            sep = ',', row.names = F)
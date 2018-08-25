# -----------------------------------------------------------------------------
# Script 04: Human Plasma Lasso Biomarker Analysis
# Author: Rani Powers & Jim Costello
# Last updated: August 1, 2018
#
# Removes batch effects from Nexus & HTP metabolomics data and plots ROC with 
# each metabolite as an individual predictor in the Results/Figure_2/ folder.
#
# For multivariate analysis, this script then imputes missing data within each 
# karyotype, z-scores the whole data set, and uses NROUNDS (10,000) CV lasso
# models to plot distribution of coefficients of the most frequently selected 
# features.
# 
# A metabolic signature of trisomy 21 is learned and the ROC plot is saved to
# the Results/Figure_2/ folder.
#
# NOTE: Setting NROUNDS to 10,000 (the number of models used in our paper) will
# result in this script running for a few hours. To decrease run time, set
# NROUNDS to lower number, e.g. 1000.
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
if (!dir.exists('Results/Figure_2')) { dir.create('Results/Figure_2', recursive = T) }
source('R/helpers.R')
source('R/lasso.R')
NROUNDS = 1000

# Use the removeBatchEffect data instead of residuals
eset = readRDS('Data/Human_plasma_metabolomics_filtered.rds')
unadj_met_data = exprs(eset)
sample_data = pData(eset)
feature_data = fData(eset)
adj_met_data = removeBatchEffect(unadj_met_data, 
                                 batch = sample_data$Batch, 
                                 covariates = cbind(sample_data$Age, as.factor(sample_data$Sex)))

# Get the 16 significantly differentially abundant metabolites
sig_table = read.csv('Results/Figure_1/Fig1_SuppTable.csv', stringsAsFactors = F)
de_mets = sig_table[sig_table$adj.P.Val < .01, 'Compound_ID']

# ---------------------------- UNIVARIATE ROC ---------------------------------
# ROC plot for individual metabolite predictors
sig_met_data = adj_met_data[de_mets,]
karyo_flip = ifelse(sample_data[colnames(sig_met_data), 'Karyotype'] == 'T21', 'D21', 'T21')
pdf('Results/Figure_2/Fig2A.pdf',
    height = 10, width = 10)
met_cols = c(brewer.pal(11, 'Spectral'), brewer.pal(5, 'Set1'))
plot(0, 0, col = 'white', ylim = c(0, 1), xlim = c(0, 1), 
     ylab = 'True positive rate', xlab = 'False positive rate', 
     main = 'Performance of individual metabolites as predictors')
met_names = c()
aucs = c()
for(i in 1:nrow(sig_met_data)) {
  pred = prediction(sig_met_data[i,], sample_data[colnames(sig_met_data), 'Karyotype'])
  perf_auc = performance(pred, 'auc')
  if (perf_auc@y.values[[1]] < 0.5) {
    pred = prediction(sig_met_data[i,], karyo_flip)
    perf_auc = performance(pred, 'auc')
  }
  aucs = c(aucs, perf_auc@y.values[[1]])
  perf = performance(pred, 'tpr', 'fpr')
  lines(perf@x.values[[1]], perf@y.values[[1]], col = met_cols[i])
  segments(0, 0, 1, 1, col = 'grey')
  met_names = c(met_names, 
                paste0(feature_data$Compound_Name[match(row.names(sig_met_data)[i],feature_data$Compound_ID)], 
                      ', AUC = ', round(perf_auc@y.values[[1]], 2)))
}
legend('bottomright', legend = met_names, lty = 1, col = met_cols)
dev.off() 

# -------------------------- MULTIVARIATE LASSO -------------------------------
# Impute missing values within karyotype
t21 = sample_data[sample_data$Karyotype == 'T21', 'Barcode']
d21 = sample_data[sample_data$Karyotype == 'D21', 'Barcode']
t21_data = impute.knn(adj_met_data[,t21], k = 5)$data
d21_data = impute.knn(adj_met_data[,d21], k = 5)$data
met_data = cbind(t21_data, d21_data)

# Z-score by row (metabolite) across whole data set
met_data_z = as.data.frame(scale(t(met_data)))
met_data_z$Karyotype = sample_data[row.names(met_data_z), 'Karyotype']

# Run cross-validated lasso on whole data set
lasso_result = LassoTally(train.dataset = met_data_z, good.acc = 0, 
                           pos.class.label = 'T21', n.rounds = NROUNDS)
saveRDS(lasso_result, 'Data/lasso_result.rds')

# Find the top 80% of metabolites selected in models with > 90% accuracy
lasso_result = readRDS('Data/lasso_result.rds')
mean_auc =  sort(lasso_result$AUC)
quantile(mean_auc, probs = 0.025)
quantile(mean_auc, probs = 0.975)
hist(mean_auc)
mean(mean_auc)
high_acc = lasso_result$Accuracy > 0.9
coeffs_df = lasso_result$AllCoeffsDF[-1,]

# Extract just the coefficient distributions for the most accurate models
good_df = cbind(Original_ID = coeffs_df$Original_ID, coeffs_df[,-1][,high_acc])
good_df$Original_ID = feature_data[as.character(good_df$Original_ID), 'Compound_Name']
feature_tally = apply(good_df[,-1], 1, function(x) sum(!is.na(x)))
names(feature_tally) = good_df$Original_ID
features = as.data.frame(feature_tally)
features$Compound_Name = row.names(features)
features$Compound_ID = sapply(features$Compound_Name, function(x){
  feature_data[feature_data$Compound_Name == x, 'Compound_ID']
})
features = features[order(features$feature_tally, decreasing = T),]

# 58 models had accuracy > 90%. Get the top 80% of features (11) in these models.
top_80percent = sum(lasso_result$Accuracy >= .9) * .8
top_80percent_features = features[features$feature_tally > top_80percent,]
top_coeffs_df = good_df[good_df$Original_ID %in% top_80percent_features$Compound_Name,]

pdf('Results/Figure_2/Fig2B.pdf', height = 6, width = 7)
par(mar = c(6,15,5,2))
CoeffsDFBoxplot(top_coeffs_df)
dev.off()

# Use just these top 11 features in a lasso model
biomarkers = top_80percent_features$Compound_ID
biomarker_lasso_result = LassoTally(train.dataset = met_data_z[,c(biomarkers, 'Karyotype')], 
                                    good.acc = 0, pos.class.label = 'T21', 
                                    n.rounds = NROUNDS)
saveRDS(biomarker_lasso_result, 'Data/biomarker_lasso_result.rds')

biomarker_lasso_result = readRDS('Data/biomarker_lasso_result.rds')
mean_auc =  sort(biomarker_lasso_result$AUC)
quantile(mean_auc, probs = 0.025)
quantile(mean_auc, probs = 0.975)
hist(mean_auc)
mean(mean_auc)

# --------------------------- MULTIVARIATE ROC --------------------------------

# ROC plot for biomarker lasso models
pdf('Results/Figure_2/Fig2C.pdf',
    height = 10, width = 10)
PlotROCWithCI(biomarker_lasso_result)
dev.off() 
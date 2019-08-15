# -----------------------------------------------------------------------------
# Script 01: Human Plasma Differential Abundance Analysis
# Author: Rani Powers
#
# Reads in the human metabolomics eset and performs differential abundance 
# analysis for Cohort 1 alone, Cohort 2 alone, and both cohorts combined
# (linear model fit for batch, age, and sex from 00-model_fitting.R).
#
# Outputs:
# 1) Saves volcano plot in the Results/Figure1/ folder for Cohort 1 and
# Cohort 2 combined
#
# 2) Saves the differential abundance results for Cohort 1 alone, Cohort 2 alone
# and combined, with and without covariates, in Results/Supplementary_Files/
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
FDR_CUTOFF = .05

# Read in the raw data (no model fitting has been performed yet)
raw_eset = readRDS('Data/Human_plasma_metabolomics_eset_filtered.rds')
met_data = as.data.frame(exprs(raw_eset))
sample_data = pData(raw_eset)
feature_data = fData(raw_eset)

# ----------------------- COHORT 1 NO FITTING ----------------------------------

NEXUS = row.names(sample_data[sample_data$Batch_Name == 'Nexus',])
nexus_data = met_data[,NEXUS]
nexus_sample_data = sample_data[NEXUS,]

# Cohort 1 differential abundance - no fitting
nexus_design = model.matrix(~0 + as.factor(Karyotype),
                            data = nexus_sample_data)
colnames(nexus_design) = c('D21', 'T21')
nexus_contrast = makeContrasts(T21vsD21 = T21-D21, levels = nexus_design)

# Make table for diff abundance results
nexus_fit = eBayes(contrasts.fit(lmFit(nexus_data, nexus_design), nexus_contrast))
nexus_sigtable = topTable(nexus_fit, adjust = 'fdr', number = nrow(nexus_data)+1,
                          confint = T)
nexus_sigtable$Compound_ID = row.names(nexus_sigtable)
nexus_sigtable = join(nexus_sigtable, feature_data)

# Calculate unmoderated t-test
nexus_fit_unmod = contrasts.fit(lmFit(nexus_data, nexus_design), nexus_contrast)
unmod_t = (nexus_fit_unmod$coefficients / nexus_fit_unmod$stdev.unscaled / nexus_fit_unmod$sigma)[,1]
unmod_df = data.frame(Unmodified_t = unmod_t,
                      Unmodified_pvalue = 2*pt(abs(unmod_t), length(unmod_t)-1, lower=F),
                      Unmodified_DF = nexus_fit_unmod$df.residual)
unmod_df$Compound_ID = row.names(unmod_df)
unmod_df$Unmodified_fdr = p.adjust(unmod_df$Unmodified_pvalue, method = 'fdr')
nexus_sigtable = join(nexus_sigtable, unmod_df)
nexus_sigtable = nexus_sigtable[,c('Compound_ID', 'Compound_Name', 'logFC',
                                   'CI.L', 'CI.R', 'AveExpr', 
                                   'P.Value', 'adj.P.Val', 'Unmodified_t',
                                   'Unmodified_pvalue', 'Unmodified_fdr', 'Unmodified_DF')]
names(nexus_sigtable)[7:8] = c('Modified_pvalue', 'Modified_fdr')

# ---------------------- COHORT 1 AGE & SEX FIT --------------------------------

nexus_design2 = model.matrix(~0 + Karyotype+Age+Sex, 
                             nexus_sample_data[,c('Karyotype', 'Age', 'Sex')])
colnames(nexus_design2)[1:2] = c('D21', 'T21')
nexus_contrast2 = makeContrasts(T21vsD21 = T21-D21, levels = nexus_design2)

# Make table for diff abundance results
nexus_fit2 = eBayes(contrasts.fit(lmFit(nexus_data, nexus_design2), nexus_contrast2))
nexus_sigtable2 = topTable(nexus_fit2, adjust = 'fdr', number = nrow(nexus_data)+1,
                           confint = T)
nexus_sigtable2$Compound_ID = row.names(nexus_sigtable2)
nexus_sigtable2 = join(nexus_sigtable2, feature_data)

# Calculate unmoderated t-test
nexus_fit2_unmod = contrasts.fit(lmFit(nexus_data, nexus_design2), nexus_contrast2)
unmod_t2 = (nexus_fit2_unmod$coefficients / nexus_fit2_unmod$stdev.unscaled / nexus_fit2_unmod$sigma)[,1]
unmod_df2 = data.frame(Unmodified_t = unmod_t2,
                       Unmodified_pvalue = 2*pt(abs(unmod_t2), length(unmod_t2)-1, lower=F),
                       Unmodified_DF = nexus_fit2_unmod$df.residual)
unmod_df2$Compound_ID = row.names(unmod_df2)
unmod_df2$Unmodified_fdr = p.adjust(unmod_df2$Unmodified_pvalue, method = 'fdr')
nexus_sigtable2 = join(nexus_sigtable2, unmod_df2)
nexus_sigtable2 = nexus_sigtable2[,c('Compound_ID', 'Compound_Name', 'logFC',
                                     'CI.L', 'CI.R', 'AveExpr', 
                                     'P.Value', 'adj.P.Val', 'Unmodified_t',
                                     'Unmodified_pvalue', 'Unmodified_fdr', 'Unmodified_DF')]
names(nexus_sigtable2)[7:8] = c('Modified_pvalue', 'Modified_fdr')

# ------------------------ COHORT 1 SUPP TABLE --------------------------------- 

# Save all Cohort 1 only data
nexus_sigtable2 = nexus_sigtable2[,-2]
names(nexus_sigtable2)[-1] = sapply(names(nexus_sigtable2)[-1], function(s){
  paste0('Fit_', s)
})
all_cohort1 = join(nexus_sigtable, nexus_sigtable2)

# Format file
all_cohort1 = all_cohort1[c(1:6,9,12,10,11,
                            13:16,19,22,20,21)]
names(all_cohort1) = c('Compound_ID', 'Compound_Name',
                       'Log2_FC', 'Lower_CI', 'Upper_CI',
                       'Average_Intensity', 't', 'DF',
                       'p_value', 'Adjusted_p_value',
                       'WithCovars_Log2_FC', 'WithCovars_Lower_CI', 'WithCovars_Upper_CI',
                       'WithCovars_Average_Intensity', 'WithCovars_t', 'WithCovars_DF',
                       'WithCovars_p_value', 'WithCovars_Adjusted_p_value')

write.table(all_cohort1, 
            'Results/Supplementary_Files/Supplementary_Data_4B.csv',
            sep = ',', row.names = F)

# ------------------------ COHORT 2 NO FITTING --------------------------------- 

HTP = row.names(sample_data[sample_data$Batch_Name == 'HTP',])
htp_data = met_data[,HTP]
htp_sample_data = sample_data[HTP,]

# HTP differential abundance
htp_design = model.matrix(~0 + as.factor(Karyotype),
                          data = htp_sample_data)
colnames(htp_design) = c('D21', 'T21')
htp_contrast = makeContrasts(T21vsD21 = T21-D21, levels = htp_design)

# Make table for diff abundance results
htp_fit = eBayes(contrasts.fit(lmFit(htp_data, htp_design), htp_contrast))
htp_sigtable = topTable(htp_fit, adjust = 'fdr', number = nrow(htp_data)+1,
                        confint = T)
htp_sigtable$Compound_ID = row.names(htp_sigtable)
htp_sigtable = join(htp_sigtable, feature_data)

# Calculate unmoderated t-test
htp_fit_unmod = contrasts.fit(lmFit(htp_data, htp_design), htp_contrast)
unmod_t = (htp_fit_unmod$coefficients / htp_fit_unmod$stdev.unscaled / htp_fit_unmod$sigma)[,1]
unmod_df = data.frame(Unmodified_t = unmod_t,
                      Unmodified_pvalue = 2*pt(abs(unmod_t), length(unmod_t)-1, lower=F),
                      Unmodified_DF = htp_fit_unmod$df.residual)
unmod_df$Compound_ID = row.names(unmod_df)
unmod_df$Unmodified_fdr = p.adjust(unmod_df$Unmodified_pvalue, method = 'fdr')
htp_sigtable = join(htp_sigtable, unmod_df)
htp_sigtable = htp_sigtable[,c('Compound_ID', 'Compound_Name', 'logFC',
                               'CI.L', 'CI.R', 'AveExpr', 
                               'P.Value', 'adj.P.Val', 'Unmodified_t',
                               'Unmodified_pvalue', 'Unmodified_fdr', 'Unmodified_DF')]
names(htp_sigtable)[7:8] = c('Modified_pvalue', 'Modified_fdr')

# ------------------------ COHORT 2 AGE & SEX FIT ------------------------------ 

# Make table for diff abundance results
htp_design2 = model.matrix(~0 + Karyotype+Age+Sex, 
                           htp_sample_data[,c('Karyotype', 'Age', 'Sex')])
colnames(htp_design2)[1:2] = c('D21', 'T21')
htp_contrast2 = makeContrasts(T21vsD21 = T21-D21, levels = htp_design2)
htp_fit2 = eBayes(contrasts.fit(lmFit(htp_data, htp_design2), htp_contrast2))
htp_sigtable2 = topTable(htp_fit2, adjust = 'fdr', number = nrow(htp_data)+1,
                         confint = T)
htp_sigtable2$Compound_ID = row.names(htp_sigtable2)
htp_sigtable2 = join(htp_sigtable2, feature_data)

# Calculate unmoderated t-test
htp_fit2_unmod = contrasts.fit(lmFit(htp_data, htp_design2), htp_contrast2)
unmod_t2 = (htp_fit2_unmod$coefficients / htp_fit2_unmod$stdev.unscaled / htp_fit2_unmod$sigma)[,1]
unmod_df2 = data.frame(Unmodified_t = unmod_t2,
                       Unmodified_pvalue = 2*pt(abs(unmod_t2), length(unmod_t2)-1, lower=F),
                       Unmodified_DF = htp_fit2_unmod$df.residual)
unmod_df2$Compound_ID = row.names(unmod_df2)
unmod_df2$Unmodified_fdr = p.adjust(unmod_df2$Unmodified_pvalue, method = 'fdr')
htp_sigtable2 = join(htp_sigtable2, unmod_df2)
htp_sigtable2 = htp_sigtable2[,c('Compound_ID', 'Compound_Name', 'logFC',
                                 'CI.L', 'CI.R', 'AveExpr', 
                                 'P.Value', 'adj.P.Val', 'Unmodified_t',
                                 'Unmodified_pvalue', 'Unmodified_fdr', 'Unmodified_DF')]
names(htp_sigtable2)[7:8] = c('Modified_pvalue', 'Modified_fdr')

# ------------------------ COHORT 2 SUPP TABLE --------------------------------- 

# Save all Cohort 2 only data
htp_sigtable2 = htp_sigtable2[,-2]
names(htp_sigtable2)[-1] = sapply(names(htp_sigtable2)[-1], function(s){
  paste0('Fit_', s)
})
all_cohort2 = join(htp_sigtable, htp_sigtable2)

# Format file
all_cohort2 = all_cohort2[c(1:6,9,12,10,11,
                            13:16,19,22,20,21)]
names(all_cohort2) = c('Compound_ID', 'Compound_Name',
                       'Log2_FC', 'Lower_CI', 'Upper_CI',
                       'Average_Intensity', 't', 'DF',
                       'p_value', 'Adjusted_p_value',
                       'WithCovars_Log2_FC', 'WithCovars_Lower_CI', 'WithCovars_Upper_CI',
                       'WithCovars_Average_Intensity', 'WithCovars_t', 'WithCovars_DF',
                       'WithCovars_p_value', 'WithCovars_Adjusted_p_value')

write.table(all_cohort2, 
            'Results/Supplementary_Files/Supplementary_Data_4C.csv',
            sep = ',', row.names = F)

# ----------------------- COHORT 1 & 2 BATCH FIT ONLY -------------------------- 

# BOTH differential abundance
eset = readRDS('Data/Human_plasma_metabolomics_eset_filtered.rds')
met_data = as.data.frame(exprs(eset))
sample_data = pData(eset)
feature_data = fData(eset)
both_design = model.matrix(~0 + Karyotype+Batch, 
                           sample_data[,c('Karyotype', 'Batch')])
colnames(both_design)[1:2] = c('D21', 'T21')

# Make table for diff abundance results
both_contrast = makeContrasts(T21vsD21 = T21-D21, levels = both_design)
both_fit = eBayes(contrasts.fit(lmFit(met_data, both_design), both_contrast))
both_sigtable = topTable(both_fit, adjust = 'fdr', number = nrow(met_data)+1,
                         confint = T)
both_sigtable$Compound_ID = row.names(both_sigtable)
both_sigtable = join(both_sigtable, feature_data)

# Calculate unmoderated t-test
both_fit_unmod = contrasts.fit(lmFit(met_data, both_design), both_contrast)
unmod_t = (both_fit_unmod$coefficients / both_fit_unmod$stdev.unscaled / both_fit_unmod$sigma)[,1]
unmod_df = data.frame(Unmodified_t = unmod_t,
                      Unmodified_pvalue = 2*pt(abs(unmod_t), length(unmod_t)-1, lower=F),
                      Unmodified_DF = both_fit_unmod$df.residual)
unmod_df$Compound_ID = row.names(unmod_df)
unmod_df$Unmodified_fdr = p.adjust(unmod_df$Unmodified_pvalue, method = 'fdr')
both_sigtable = join(both_sigtable, unmod_df)
both_sigtable = both_sigtable[,c('Compound_ID', 'Compound_Name', 'logFC',
                                 'CI.L', 'CI.R', 'AveExpr', 
                                 'P.Value', 'adj.P.Val', 'Unmodified_t',
                                 'Unmodified_pvalue', 'Unmodified_fdr', 'Unmodified_DF')]
names(both_sigtable)[7:8] = c('Modified_pvalue', 'Modified_fdr')

# -------------------- COHORT 1 & 2 BATCH, AGE SEX FIT ------------------------- 

both_design2 = model.matrix(~0 + Karyotype+Age+Sex+Batch, 
                            sample_data[,c('Karyotype', 'Age', 'Sex', 'Batch')])
colnames(both_design2)[1:2] = c('D21', 'T21')

# Make table for diff abundance results
both_contrast2 = makeContrasts(T21vsD21 = T21-D21, levels = both_design2)
both_fit2 = eBayes(contrasts.fit(lmFit(met_data, both_design2), both_contrast2))
both_sigtable2 = topTable(both_fit2, adjust = 'fdr', number = nrow(met_data)+1,
                          confint = T)
both_sigtable2$Compound_ID = row.names(both_sigtable2)
both_sigtable2 = join(both_sigtable2, feature_data)

# Calculate unmoderated t-test
both_fit2_unmod = contrasts.fit(lmFit(met_data, both_design2), both_contrast2)
unmod_t2 = (both_fit2_unmod$coefficients / both_fit2_unmod$stdev.unscaled / both_fit2_unmod$sigma)[,1]
unmod_df2 = data.frame(Unmodified_t = unmod_t2,
                       Unmodified_pvalue = 2*pt(abs(unmod_t2), length(unmod_t2)-1, lower=F),
                       Unmodified_DF = both_fit2_unmod$df.residual)
unmod_df2$Compound_ID = row.names(unmod_df2)
unmod_df2$Unmodified_fdr = p.adjust(unmod_df2$Unmodified_pvalue, method = 'fdr')
both_sigtable2 = join(both_sigtable2, unmod_df2)
both_sigtable2 = both_sigtable2[,c('Compound_ID', 'Compound_Name', 'logFC',
                                   'CI.L', 'CI.R', 'AveExpr', 
                                   'P.Value', 'adj.P.Val', 'Unmodified_t',
                                   'Unmodified_pvalue', 'Unmodified_fdr', 'Unmodified_DF')]
names(both_sigtable2)[7:8] = c('Modified_pvalue', 'Modified_fdr')

# ------------------------- COHORT 1 & 2 PLOTTING ------------------------------ 

# Plot volcano with batch, age, sex fitting
pdf(paste0('Results/Figure_1/Fig_1A_volcano_bothCohorts_withFit.pdf'),
    height = 7, width = 5)
plot_volcano(x = both_sigtable2$logFC, 
             y = -log10(both_sigtable2$Unmodified_fdr), 
             plot_fc_lines = F,
             y_cutoff = -log10(.05),
             x_cutoffs = c(-.01, .01), xlim = c(-2,2),
             main = 'Combined Nexus & HTP\nmodel fit with batch, age and sex', 
             pt_labels = both_sigtable2$Compound_Name)
dev.off()

# ------------------------ COHORT 2 SUPP TABLE --------------------------------- 

# Save all Cohort 1 and 2 data
both_sigtable2 = both_sigtable2[,-2]
names(both_sigtable2)[-1] = sapply(names(both_sigtable2)[-1], function(s){
  paste0('Fit_', s)
})
all_cohorts = join(both_sigtable, both_sigtable2)

# Format file
all_cohorts = all_cohorts[c(1:6,9,12,10,11,
                            13:16,19,22,20,21)]
names(all_cohorts) = c('Compound_ID', 'Compound_Name',
                       'Log2_FC', 'Lower_CI', 'Upper_CI',
                       'Average_Intensity', 't', 'DF',
                       'p_value', 'Adjusted_p_value',
                       'WithCovars_Log2_FC', 'WithCovars_Lower_CI', 'WithCovars_Upper_CI',
                       'WithCovars_Average_Intensity', 'WithCovars_t', 'WithCovars_DF',
                       'WithCovars_p_value', 'WithCovars_Adjusted_p_value')

write.table(all_cohorts, 
            'Results/Supplementary_Files/Supplementary_Data_4A.csv',
            sep = ',', row.names = F)

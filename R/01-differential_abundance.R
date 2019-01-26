# -----------------------------------------------------------------------------
# Script 01: Human Plasma Differential Abundance Analysis
# Author: Rani Powers
# Last updated: December 11, 2018
#
# Reads in Data/human_plasma_metabolomics_eset_filtered.rds and performs 
# differential abundance analysis on all data using a linear model fit with 
# batch, age, and sex as covariates.
#
# Outputs:
# 1) Saves Fig 1a volcano plot for combined Cohort 1 & 2 data in the 
# Results/Figure1/ directory
#
# 2) Saves the differential abundance results for Cohort 1 alone, Cohort 2
# alone, and both combined in Results/Supplementary_Files/Supp_Data_2...txt
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
FDR_CUTOFF = .05

# Read in the raw data (no model fitting has been performed yet)
raw_eset = readRDS('Data/human_plasma_metabolomics_eset_cleaned.rds')
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

# Make table for diff abundance results with unmoderated t-test 
nexus_fit_unmod = contrasts.fit(lmFit(nexus_data, nexus_design), nexus_contrast)
unmod_t = (nexus_fit_unmod$coefficients / nexus_fit_unmod$stdev.unscaled / nexus_fit_unmod$sigma)[,1]
unmod_df = data.frame(Unmod_t = unmod_t,
                      Unmod_Pvalue = 2*pt(abs(unmod_t), length(unmod_t)-1, lower=F))
unmod_df$Compound_ID = row.names(unmod_df)
unmod_df$Unmod_adjP = p.adjust(unmod_df$Unmod_Pvalue, method = 'fdr')
unmod_df$log2FC = sapply(unmod_df$Compound_ID, function(met){
  mean(as.numeric(nexus_data[met, nexus_sample_data$Karyotype == 'T21']), na.rm = T)-
    mean(as.numeric(nexus_data[met, nexus_sample_data$Karyotype == 'D21']), na.rm = T)
})
nexus_sigtable = join(unmod_df, feature_data[,1:2])
nexus_sigtable = nexus_sigtable[,c(3,6,5,2,4)]

# ---------------------- COHORT 1 AGE & SEX FIT --------------------------------

nexus_design2 = model.matrix(~0 + Karyotype+Age+Sex, 
                             nexus_sample_data[,c('Karyotype', 'Age', 'Sex')])
colnames(nexus_design2)[1:2] = c('D21', 'T21')
nexus_contrast2 = makeContrasts(T21vsD21 = T21-D21, levels = nexus_design2)

# Make table for diff abundance results using unmoderated t-test
nexus_fit = topTable(eBayes(contrasts.fit(lmFit(nexus_data, nexus_design2), nexus_contrast2)),
                   number = nrow(nexus_data))
nexus_fit = data.frame(Compound_ID = row.names(nexus_fit), log2FC_withCovars = nexus_fit$logFC)
nexus_fit_unmod2 = contrasts.fit(lmFit(nexus_data, nexus_design2), nexus_contrast2)
unmod_t2 = (nexus_fit_unmod2$coefficients / nexus_fit_unmod2$stdev.unscaled / nexus_fit_unmod2$sigma)[,1]
unmod_df2 = data.frame(Unmod_t_withCovars = unmod_t2,
                       Unmod_Pvalue_withCovars = 2*pt(abs(unmod_t2), length(unmod_t2)-1, lower=F))
unmod_df2$Compound_ID = row.names(unmod_df2)
unmod_df2$Unmod_adjP_withCovars = p.adjust(unmod_df2$Unmod_Pvalue_withCovars, method = 'fdr')
unmod_df2 = join(unmod_df2, nexus_fit)
nexus_sigtable = join(nexus_sigtable, unmod_df2[,c(3,5,2,4)])

# Save Cohort 1 only data
nexus_sigtable = nexus_sigtable[order(nexus_sigtable$Unmod_adjP_withCovars, decreasing = F),]
write.table(nexus_sigtable, 'Results/Supplementary_Files/Supp_Data_2b_cohort_1_only.csv',
            sep = ',', row.names = F)

# ------------------------ COHORT 2 NO FITTING --------------------------------- 

HTP = row.names(sample_data[sample_data$Batch_Name == 'HTP',])
htp_data = met_data[,HTP]
htp_sample_data = sample_data[HTP,]

# Cohort 2 differential abundance - no fitting
htp_design = model.matrix(~0 + as.factor(Karyotype),
                            data = htp_sample_data)
colnames(htp_design) = c('D21', 'T21')
htp_contrast = makeContrasts(T21vsD21 = T21-D21, levels = htp_design)

# Make table for diff abundance results with unmoderated t-test 
htp_fit_unmod = contrasts.fit(lmFit(htp_data, htp_design), htp_contrast)
unmod_t = (htp_fit_unmod$coefficients / htp_fit_unmod$stdev.unscaled / htp_fit_unmod$sigma)[,1]
unmod_df = data.frame(Unmod_t = unmod_t,
                      Unmod_Pvalue = 2*pt(abs(unmod_t), length(unmod_t)-1, lower=F))
unmod_df$Compound_ID = row.names(unmod_df)
unmod_df$Unmod_adjP = p.adjust(unmod_df$Unmod_Pvalue, method = 'fdr')
unmod_df$log2FC = sapply(unmod_df$Compound_ID, function(met){
  mean(as.numeric(htp_data[met, htp_sample_data$Karyotype == 'T21']), na.rm = T)-
    mean(as.numeric(htp_data[met, htp_sample_data$Karyotype == 'D21']), na.rm = T)
})
htp_sigtable = join(unmod_df, feature_data[,1:2])
htp_sigtable = htp_sigtable[,c(3,6,5,2,4)]

# ---------------------- COHORT 1 AGE & SEX FIT --------------------------------

htp_design2 = model.matrix(~0 + Karyotype+Age+Sex, 
                             htp_sample_data[,c('Karyotype', 'Age', 'Sex')])
colnames(htp_design2)[1:2] = c('D21', 'T21')
htp_contrast2 = makeContrasts(T21vsD21 = T21-D21, levels = htp_design2)

# Make table for diff abundance results using unmoderated t-test
htp_fit = topTable(eBayes(contrasts.fit(lmFit(htp_data, htp_design2), htp_contrast2)),
                   number = nrow(htp_data))
htp_fit = data.frame(Compound_ID = row.names(htp_fit), log2FC_withCovars = htp_fit$logFC)
htp_fit_unmod2 = contrasts.fit(lmFit(htp_data, htp_design2), htp_contrast2)
unmod_t2 = (htp_fit_unmod2$coefficients / htp_fit_unmod2$stdev.unscaled / htp_fit_unmod2$sigma)[,1]
unmod_df2 = data.frame(Unmod_t_withCovars = unmod_t2,
                       Unmod_Pvalue_withCovars = 2*pt(abs(unmod_t2), length(unmod_t2)-1, lower=F))
unmod_df2$Compound_ID = row.names(unmod_df2)
unmod_df2$Unmod_adjP_withCovars = p.adjust(unmod_df2$Unmod_Pvalue_withCovars, method = 'fdr')
unmod_df2 = join(unmod_df2, htp_fit)
htp_sigtable = join(htp_sigtable, unmod_df2[,c(3,5,2,4)])

# Save Cohort 1 only data
htp_sigtable = htp_sigtable[order(htp_sigtable$Unmod_adjP_withCovars, decreasing = F),]
write.table(htp_sigtable, 'Results/Supplementary_Files/Supp_Data_2c_cohort_2_only.csv',
            sep = ',', row.names = F)

# ----------------------- COHORT 1 & 2 BATCH FIT ONLY -------------------------- 

# BOTH differential abundance
eset = readRDS('Data/human_plasma_metabolomics_eset_cleaned.rds')
both_data = as.data.frame(exprs(eset))
sample_data = pData(eset)
feature_data = fData(eset)
both_design = model.matrix(~0 + Karyotype+Batch, 
                           sample_data[,c('Karyotype', 'Batch')])
colnames(both_design)[1:2] = c('D21', 'T21')
both_contrast = makeContrasts(T21vsD21 = T21-D21, levels = both_design)

# Make table for diff abundance results using unmoderated t-test
both_fit = topTable(eBayes(contrasts.fit(lmFit(both_data, both_design), both_contrast)),
                    number = nrow(both_data))
both_fit = data.frame(Compound_ID = row.names(both_fit), log2FC_withCovars = both_fit$logFC)
both_fit_unmod = contrasts.fit(lmFit(both_data, both_design), both_contrast)
unmod_t = (both_fit_unmod$coefficients / both_fit_unmod$stdev.unscaled / both_fit_unmod$sigma)[,1]
unmod_df = data.frame(Unmod_t = unmod_t,
                      Unmod_Pvalue = 2*pt(abs(unmod_t), length(unmod_t)-1, lower=F))
unmod_df$Compound_ID = row.names(unmod_df)
unmod_df$Unmod_adjP = p.adjust(unmod_df$Unmod_Pvalue, method = 'fdr')
unmod_df$log2FC = sapply(unmod_df$Compound_ID, function(met){
  mean(as.numeric(both_data[met, sample_data$Karyotype == 'T21']), na.rm = T)-
    mean(as.numeric(both_data[met, sample_data$Karyotype == 'D21']), na.rm = T)
})
both_sigtable = join(unmod_df, feature_data[,1:2])
both_sigtable = both_sigtable[,c(3,6,5,2,4)]

# -------------------- COHORT 1 & 2 BATCH, AGE SEX FIT ------------------------- 

both_design2 = model.matrix(~0 + Karyotype+Age+Sex+Batch, 
                            sample_data[,c('Karyotype', 'Age', 'Sex', 'Batch')])
colnames(both_design2)[1:2] = c('D21', 'T21')
both_contrast2 = makeContrasts(T21vsD21 = T21-D21, levels = both_design2)

# Make table for diff abundance results using unmoderated t-test
both_fit2 = topTable(eBayes(contrasts.fit(lmFit(both_data, both_design2), both_contrast2)),
                   number = nrow(both_data))
both_fit2 = data.frame(Compound_ID = row.names(both_fit2), log2FC_withCovars = both_fit2$logFC)
both_fit_unmod2 = contrasts.fit(lmFit(both_data, both_design2), both_contrast2)
unmod_t2 = (both_fit_unmod2$coefficients / both_fit_unmod2$stdev.unscaled / both_fit_unmod2$sigma)[,1]
unmod_df2 = data.frame(Unmod_t_withCovars = unmod_t2,
                       Unmod_Pvalue_withCovars = 2*pt(abs(unmod_t2), length(unmod_t2)-1, lower=F))
unmod_df2$Compound_ID = row.names(unmod_df2)
unmod_df2$Unmod_adjP_withCovars = p.adjust(unmod_df2$Unmod_Pvalue_withCovars, method = 'fdr')
unmod_df2 = join(unmod_df2, both_fit)
both_sigtable = join(both_sigtable, unmod_df2[,c(3,5,2,4)])

# Save Cohort 1 and 2 only data
both_sigtable = both_sigtable[order(both_sigtable$Unmod_adjP_withCovars, decreasing = F),]
write.table(both_sigtable, 'Results/Supplementary_Files/Supp_Data_2a_cohort_1_and_2.csv',
            sep = ',', row.names = F)

# ------------------------- COHORT 1 & 2 PLOTTING ------------------------------ 

# Plot volcano with batch, age, sex fitting
pdf(paste0('Results/Figure_1/Fig_1a_volcano_both_cohorts_with_fit.pdf'),
    height = 7, width = 5)
plot_volcano(x = both_sigtable$log2FC_withCovars, 
             y = -log10(both_sigtable$Unmod_adjP_withCovars), 
             plot_fc_lines = F,
             y_cutoff = -log10(.05),
             x_cutoffs = c(-.01, .01), xlim = c(-2,2),
             main = 'Combined Cohort 1 & Cohort 2\nmodel fit with batch, age and sex', 
             pt_labels = both_sigtable$Compound_Name)
dev.off()

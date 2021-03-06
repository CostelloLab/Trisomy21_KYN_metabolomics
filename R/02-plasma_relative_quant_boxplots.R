# -----------------------------------------------------------------------------
# Script 02: Human Plasma T21 vs D21 Boxplots (relative quant)
# Author: Rani Powers
#
# Saves boxplots for kynurenine, quinolinic acid, KYN/TRP, and QA/PA ratios in
# Results/Figure_1/ folder. Additional boxplots for supp figure 1 are saved in 
# the Results/Supplementary_Figures/ folder.
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
FDR_CUTOFF = .05

# Get the raw data for the 91 metabolites
eset = readRDS('Data/Human_plasma_metabolomics_eset_filtered.rds')
met_raw = as.data.frame(t(exprs(eset)))
sample_data = pData(eset)
feature_data = fData(eset)
feature_data = rbind.fill(feature_data,
                          data.frame(Compound_ID = 'Kyn_Tryp',
                                     Compound_Name = 'Kynurenine/tryptophan ratio'),
                          data.frame(Compound_ID = 'QA_PA',
                                     Compound_Name = 'Quinolinic acid/picolinic acid ratio'))
row.names(feature_data) = feature_data$Compound_ID

# Add ratios of interest to the log2 data
met_raw$Kyn_Tryp = met_raw$C00328 - met_raw$C00078
met_raw$QA_PA = met_raw$C03722 - met_raw$C10164
met_raw = as.data.frame(t(met_raw))

# Read in the differential results for all (non-ratio) metabolites
met_sigtable = read.csv('Results/Supplementary_Files/Supplementary_Data_4A.csv',
                        stringsAsFactors = F)
row.names(met_sigtable) = met_sigtable$Compound_ID

# Add nominal p-value to met_sigtable for ratios
qa_pa_test = t.test(met_raw['QA_PA', sample_data$Karyotype == 'T21'],
                    met_raw['QA_PA', sample_data$Karyotype == 'D21'])
kyn_tryp_test = t.test(met_raw['Kyn_Tryp', sample_data$Karyotype == 'T21'],
                       met_raw['Kyn_Tryp', sample_data$Karyotype == 'D21'])
met_sigtable = rbind.fill(met_sigtable,
                          data.frame(Compound_ID = c('QA_PA', 'Kyn_Tryp'),
                                     Compound_Name = c('QA to PA ratio', 'Kyn to Trp ratio'),
                                     Log2_FC = c(qa_pa_test$estimate[1]-qa_pa_test$estimate[2], 
                                                 kyn_tryp_test$estimate[1]-kyn_tryp_test$estimate[2]),
                                     p_value = c(qa_pa_test$p.value, kyn_tryp_test$p.value),
                                     Adjusted_p_value = c(qa_pa_test$p.value, kyn_tryp_test$p.value)))
row.names(met_sigtable) = met_sigtable$Compound_ID

# Use the removeBatchEffects() data for plotting to get adjusted intensity
met_adj = removeBatchEffect(met_raw, batch = sample_data$Batch, 
                            design = model.matrix(~0+Age+Sex, data = sample_data))
data_for_boxplots = as.data.frame(t(met_adj))
data_for_boxplots$Karyotype = sample_data[row.names(data_for_boxplots),'Karyotype']

# Plot boxplots for the metabolites in Figure 1c
fig1b_metabolites = c('C03722', 'C00328')
for (met in fig1b_metabolites){
  met_name = gsub(' |/', '_', met_sigtable[met, 'Compound_Name'])
  pdf(paste0('Results/Figure_1/Fig_1B_boxplot_', met_name, '.pdf'),
      height = 5, width = 5)
  if (met_sigtable[met, 'WithCovars_Log2_FC'] < 0 & met_sigtable[met, 'WithCovars_Adjusted_p_value'] <= FDR_CUTOFF){
    cols_to_use = c(dark_blue, light_blue)
  } else if (met_sigtable[met, 'WithCovars_Log2_FC'] > 0 & met_sigtable[met, 'WithCovars_Adjusted_p_value'] <= FDR_CUTOFF){
    cols_to_use = c(dark_red, light_red)
  } else{
    cols_to_use = c('#515151', '#A0A0A0')
  }
  plot_boxplot(met, data_for_boxplots, feature_data, 
               point_borders = 'black',
               ylab = 'batch, age and sex adjusted intensity',
               col1 = cols_to_use[1], col2 = cols_to_use[2],
               pval_column = 'WithCovars_p_value',
               qval_column = 'WithCovars_Adjusted_p_value')
  dev.off()
}

# Plot boxplots for the metabolites in Figure 1c
fig1c_metabolites = c('Kyn_Tryp', 'QA_PA')
for (met in fig1c_metabolites){
  met_name = gsub(' |/', '_', met_sigtable[met, 'Compound_Name'])
  pdf(paste0('Results/Figure_1/Fig_1C_boxplot_', met_name, '.pdf'),
      height = 5, width = 5)
  if (met_sigtable[met, 'Log2_FC'] < 0 & met_sigtable[met, 'Adjusted_p_value'] <= FDR_CUTOFF){
    cols_to_use = c(dark_blue, light_blue)
  } else if (met_sigtable[met, 'Log2_FC'] > 0 & met_sigtable[met, 'Adjusted_p_value'] <= FDR_CUTOFF){
    cols_to_use = c(dark_red, light_red)
  } else{
    cols_to_use = c('#515151', '#A0A0A0')
  }
  plot_boxplot(met, data_for_boxplots, feature_data, 
               point_borders = 'black',
               ylab = 'batch, age and sex adjusted intensity',
               col1 = cols_to_use[1], col2 = cols_to_use[2],
               pval_column = 'p_value', qval_column = 'Adjusted_p_value')
  dev.off()
}

# Plot boxplots for all Supp Fig 1 metabolites
supp_fig1d_metabolites = c('C00065', 'C00082', 'C00155', 
                           'C00042', 'C02630', 'C05422')
for (met in supp_fig1d_metabolites){
  met_name = gsub(' |/', '_', met_sigtable[met, 'Compound_Name'])
  pdf(paste0('Results/Supplementary_Figures/SuppFig_1D_boxplot_', met_name, '.pdf'),
      height = 5, width = 5)
  if (met_sigtable[met, 'WithCovars_Log2_FC'] < 0 & met_sigtable[met, 'WithCovars_Adjusted_p_value'] <= FDR_CUTOFF){
    cols_to_use = c(dark_blue, light_blue)
  } else if (met_sigtable[met, 'WithCovars_Log2_FC'] > 0 & met_sigtable[met, 'WithCovars_Adjusted_p_value'] <= FDR_CUTOFF){
    cols_to_use = c(dark_red, light_red)
  } else{
    cols_to_use = c('#515151', '#A0A0A0')
  }
  plot_boxplot(met, data_for_boxplots, feature_data, 
               point_borders = 'black',
               ylab = 'batch, age and sex adjusted intensity',
               col1 = cols_to_use[1], col2 = cols_to_use[2],
               pval_column = 'WithCovars_p_value',
               qval_column = 'WithCovars_Adjusted_p_value')
  dev.off()
}

# Plot boxplots for all Supp Fig 2 metabolites
supp_fig2_metabolites = c('C00780', 'C00078', 'C00463', 
                          'C05635', 'C00954', 'C00637',
                          'C00108', 'C10164', 'C03453')
for (met in supp_fig2_metabolites){
  met_name = gsub(' |/', '_', met_sigtable[met, 'Compound_Name'])
  pdf(paste0('Results/Supplementary_Figures/Supp_Fig_2_boxplot_', met_name, '.pdf'),
      height = 5, width = 5)
  if (met_sigtable[met, 'WithCovars_Log2_FC'] < 0 & met_sigtable[met, 'WithCovars_Adjusted_p_value'] <= FDR_CUTOFF){
    cols_to_use = c(dark_blue, light_blue)
  } else if (met_sigtable[met, 'WithCovars_Log2_FC'] > 0 & met_sigtable[met, 'WithCovars_Adjusted_p_value'] <= FDR_CUTOFF){
    cols_to_use = c(dark_red, light_red)
  } else{
    cols_to_use = c('#515151', '#A0A0A0')
  }
  plot_boxplot(met, data_for_boxplots, feature_data, 
               point_borders = 'black',
               ylab = 'batch, age and sex adjusted intensity',
               col1 = cols_to_use[1], col2 = cols_to_use[2],
               pval_column = 'WithCovars_p_value',
               qval_column = 'WithCovars_Adjusted_p_value')
  dev.off()
}

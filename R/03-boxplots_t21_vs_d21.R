# -----------------------------------------------------------------------------
# Script 03: Human Plasma T21 vs D21 Boxplots
# Author: Rani Powers
# Last updated: July 31, 2018
#
# Saves boxplots for most differentially abundant metabolites in the 
# Results/Figure_1/ folder. Additional boxplots for tryptophan biosynthesis and
# catabolism pathways are saved in the Results/Figure_3/ folder.
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
if (!dir.exists('Results/Figure_3')) { dir.create('Results/Figure_3', recursive = T) }

# Get the model residuals for the 91 metabolites
met_data = readRDS('Data/Human_plasma_residuals.rds')
eset = readRDS('Data/Human_plasma_metabolomics_filtered.rds')
sample_data = pData(eset)
feature_data = fData(eset)
data_for_boxplots = as.data.frame(t(met_data))
data_for_boxplots$Karyotype = sample_data[row.names(data_for_boxplots),'Karyotype']

# P-values etc from the differential abundance analysis
met_sigtable = read.csv('Results/Figure_1/Fig1_SuppTable.csv',
                     stringsAsFactors = F)
met_sigtable = met_sigtable[,c(7,8,1,4,5)]

# Add ratios
data_for_boxplots$Kyn_Tryp = data_for_boxplots$C00328 - data_for_boxplots$C00078
data_for_boxplots$QA_PA = data_for_boxplots$C03722 - data_for_boxplots$C10164
kyn_tryp_t = t.test(Kyn_Tryp ~ Karyotype, data = data_for_boxplots)
qa_pa_t = t.test(QA_PA ~ Karyotype, data = data_for_boxplots)
met_sigtable = rbind(met_sigtable,
                     data.frame(Compound_ID = 'Kyn_Tryp', 
                                Compound_Name = 'Kynurenine to tryptophan ratio',
                                logFC = kyn_tryp_t$estimate[2] - kyn_tryp_t$estimate[1],
                                P.Value = kyn_tryp_t$p.value,
                                adj.P.Val = p.adjust(c(kyn_tryp_t$p.value, met_sigtable$P.Value)[1])),
                     data.frame(Compound_ID = 'QA_PA', 
                                Compound_Name = 'Quinolinic acid to picolinic acid ratio',
                                logFC = qa_pa_t$estimate[2] - qa_pa_t$estimate[1],
                                P.Value = qa_pa_t$p.value,
                                adj.P.Val = p.adjust(c(qa_pa_t$p.value, met_sigtable$P.Value))[1]))
row.names(met_sigtable) = met_sigtable$Compound_ID

# Plot boxplots for the 6 metabolites in Figure 1
fig1_metabolites = c('C00065', 'C05422', 'C00155', 'C00042', 'C00082', 'C02630')
for (met in fig1_metabolites){
  met_name = gsub(' |/', '_', met_sigtable[met, 'Compound_Name'])
  pdf(paste0('Results/Figure_1/Fig1B_', met_name, '.pdf'),
      height = 5, width = 5)
  if (met_sigtable[met, 'logFC'] < 0 & met_sigtable[met, 'adj.P.Val'] <= .01){
    cols_to_use = c(dark_blue, light_blue)
  } else if (met_sigtable[met, 'logFC'] > 0 & met_sigtable[met, 'adj.P.Val'] <= .01){
    cols_to_use = c(dark_red, light_red)
  } else{
    cols_to_use = c('#515151', '#A0A0A0')
  }
  plot_boxplot(met, data_for_boxplots, feature_data, 
               ylab = 'model residuals',
               col1 = cols_to_use[1], col2 = cols_to_use[2])
  dev.off()
}

# Plot boxplots for the tryp pathway metabolites in Figure 3
fig3_metabolites = c('C00078', 'C00780', 'C00463', 'C05635', 'C00328',
                     'C00637', 'C00954', 'C00108', 'C03722', 'C03453', 'C10164')
for (met in fig3_metabolites){
  met_name = gsub(' |/', '_', met_sigtable[met, 'Compound_Name'])
  pdf(paste0('Results/Figure_3/Fig3_', met_name, '.pdf'),
      height = 5, width = 5)
  if (met_sigtable[met, 'logFC'] < 0 & met_sigtable[met, 'adj.P.Val'] <= .01){
    cols_to_use = c(dark_blue, light_blue)
  } else if (met_sigtable[met, 'logFC'] > 0 & met_sigtable[met, 'adj.P.Val'] <= .01){
    cols_to_use = c(dark_red, light_red)
  } else{
    cols_to_use = c('#515151', '#A0A0A0')
  }
  plot_boxplot(met, data_for_boxplots, feature_data, 
               ylab = 'model residuals',
               col1 = cols_to_use[1], col2 = cols_to_use[2])
  dev.off()
}

# Plot ratio boxplots for Figure 3 Supp 1
fig3_metabolites1 = c('Kyn_Tryp', 'QA_PA')
ratio_features = data.frame(Compound_ID = fig3_metabolites1,
                            Compound_Name = c('Kynurenine / tryptophan ratio',
                                              'Quinolinic acid / picolinic acid ratio'))
feature_data = rbind.fill(feature_data, ratio_features)
row.names(feature_data) = feature_data$Compound_ID
for (met in fig3_metabolites1){
  met_name = gsub(' |/', '_', met_sigtable[met, 'Compound_Name'])
  pdf(paste0('Results/Figure_3/Fig3_Supp1_', met_name, '.pdf'),
      height = 5, width = 5)
  if (met_sigtable[met, 'logFC'] < 0 & met_sigtable[met, 'adj.P.Val'] <= .01){
    cols_to_use = c(dark_blue, light_blue)
  } else if (met_sigtable[met, 'logFC'] > 0 & met_sigtable[met, 'adj.P.Val'] <= .01){
    cols_to_use = c(dark_red, light_red)
  } else{
    cols_to_use = c('#515151', '#A0A0A0')
  }
  plot_boxplot(met, data_for_boxplots, feature_data, 
               ylab = 'model residuals',
               col1 = cols_to_use[1], col2 = cols_to_use[2])
  dev.off()
}

# Plot boxplots for Figure 3 Supp 2
fig3_metabolites2 = c('C00463', 'C00064', 'C00025', 'C00022', 'C00108', 
                      'C00078', 'C00065')
for (met in fig3_metabolites2){
  met_name = gsub(' |/', '_', met_sigtable[met, 'Compound_Name'])
  pdf(paste0('Results/Figure_3/Fig3_Supp2_', met_name, '.pdf'),
      height = 5, width = 5)
  if (met_sigtable[met, 'logFC'] < 0 & met_sigtable[met, 'adj.P.Val'] <= .01){
    cols_to_use = c(dark_blue, light_blue)
  } else if (met_sigtable[met, 'logFC'] > 0 & met_sigtable[met, 'adj.P.Val'] <= .01){
    cols_to_use = c(dark_red, light_red)
  } else{
    cols_to_use = c('#515151', '#A0A0A0')
  }
  plot_boxplot(met, data_for_boxplots, feature_data, 
               ylab = 'model residuals',
               col1 = cols_to_use[1], col2 = cols_to_use[2])
  dev.off()
}

# Plot boxplots for Figure 3 Supp 1 that are only in HTP
eset = readRDS('Data/Human_plasma_metabolomics.rds')
sample_data = pData(eset)
HTP = sample_data[sample_data$Batch_Name == 'HTP', 'Barcode']
htp_data = as.data.frame(exprs(eset))[,HTP]
htp_feature_data = fData(eset)
htp_sample_data = sample_data[HTP,]

# Add ratio
htp_data = as.data.frame(log2(t(htp_data)+1))
htp_data$QA_KA = htp_data$C03722 - htp_data$C01717
htp_data = t(htp_data)

# Differential abundance to get p-values
htp_design = model.matrix(~0 + Karyotype+Age+Sex, 
                           htp_sample_data[,c('Karyotype', 'Age', 'Sex')])
colnames(htp_design)[1:2] = c('D21', 'T21')
htp_contrast = makeContrasts(T21vsD21 = T21-D21, levels = htp_design)
htp_fit = eBayes(contrasts.fit(lmFit(htp_data, htp_design), htp_contrast))
htp_sigtable = topTable(htp_fit, adjust = 'fdr', number = nrow(htp_data)+1)
htp_sigtable$Compound_ID = row.names(htp_sigtable)
htp_sigtable = join(htp_sigtable, feature_data)
row.names(htp_sigtable) = htp_sigtable$Compound_ID
htp_sigtable[htp_sigtable$Compound_ID == 'C01717', 'Compound_Name'] = 'Kynurenic acid'
htp_sigtable[htp_sigtable$Compound_ID == 'QA_KA', 'Compound_Name'] = 'Quinolinic acid to kynurenic acid ratio'

# Get model residuals for plotting
fit = lmFit(htp_data, design = htp_design)
model_residuals = residuals(fit, htp_data)
htp_residuals = as.data.frame(t(model_residuals))
htp_residuals$Karyotype = htp_sample_data[row.names(htp_residuals), 'Karyotype']

# Plot
fig3_metabolites1 = c('C01717', 'QA_KA')
feature_data = rbind.fill(feature_data,
                          data.frame(Compound_ID = fig3_metabolites1,
                                     Compound_Name = c('Kynurenic acid',
                                                       'Quinolinic acid to kynurenic acid ratio')))
row.names(feature_data) = feature_data$Compound_ID
met_sigtable = rbind(met_sigtable, htp_sigtable[fig3_metabolites1,c(7,8,1,4,5)])
for (met in fig3_metabolites1){
  met_name = gsub(' |/', '_', htp_sigtable[met, 'Compound_Name'])
  pdf(paste0('Results/Figure_3/Fig3_Supp1_', met_name, '.pdf'),
      height = 5, width = 5)
  if (htp_sigtable[met, 'logFC'] < 0 & htp_sigtable[met, 'adj.P.Val'] <= .01){
    cols_to_use = c(dark_blue, light_blue)
  } else if (htp_sigtable[met, 'logFC'] > 0 & htp_sigtable[met, 'adj.P.Val'] <= .01){
    cols_to_use = c(dark_red, light_red)
  } else{
    cols_to_use = c('#515151', '#A0A0A0')
  }
  plot_boxplot(met, df = htp_residuals, feature_data = feature_data, 
               ylab = 'model residuals',
               col1 = cols_to_use[1], col2 = cols_to_use[2])
  dev.off()
}

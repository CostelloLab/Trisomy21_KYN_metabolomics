# -----------------------------------------------------------------------------
# Script 04: Human CSF KYN, TRP, and QUIN (absolute quant)
# Author: Rani Powers
#
# Adjusts absolute quant metabolomics data for age and sex. Saves boxplots in 
# Results/Figure_1/ folder. Outputs supplementary data to the
# Results/Supplementary_Files/ folder.
# -----------------------------------------------------------------------------

source('R/helpers.R')

# Read human sample annotations
sample_data = read.table('Data/All_sample_annotations.txt',
                         sep = '\t', stringsAsFactors = F,
                         header = T, row.names = 1)

# Read human plasma absolute quant data
dat = read.csv('Data/Human_CSF_metabolomics_absolute.csv',
                   stringsAsFactors = F, row.names = 1)

# Fit model on log-transformed data (removeBatchEffect expects log values)
dat_log = log2(dat)
sample_annots = sample_data[row.names(dat_log),]

# Adjust
sample_annots$Sex = as.factor(sample_annots$Sex)
dat_adj = as.data.frame(t(removeBatchEffect(t(dat_log), 
                                            covariates = model.matrix(~0+Sex+Age, data = sample_annots),
                                            design = model.matrix(~Karyotype, data = sample_annots))))

# Undo log
dat_adj = as.data.frame(2^as.matrix(dat_adj))
sum(dat == 0) == sum(dat_adj == 0) # all zero values stayed 0 and weren't adjusted - good

# Pre-compute all p-values so that they are all corrected together
dat$Karyotype = sample_annots[row.names(dat), 'Karyotype']
dat_adj$Karyotype = sample_annots[row.names(dat_adj), 'Karyotype']
t21 = sample_annots$Karyotype == 'T21'
pvals_table = data.frame(Compound = names(dat)[-ncol(dat)],
                         RawData_Pval = sapply(names(dat)[-ncol(dat)], function(s){
                           t.test(get(s) ~ Karyotype, data = dat)$p.value
                         }),
                         RawData_T21direction = sapply(names(dat)[-ncol(dat)], function(s){
                           res = t.test(get(s) ~ Karyotype, data = dat)
                           ifelse(res$estimate[[2]] > res$estimate[[1]], 'Up', 'Down') 
                         }),
                         AdjustedData_Pval = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           t.test(get(s) ~ Karyotype, data = dat_adj)$p.value
                         }),
                         AdjustedData_T21direction = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           res = t.test(get(s) ~ Karyotype, data = dat_adj)
                           ifelse(res$estimate[[2]] > res$estimate[[1]], 'Up', 'Down')
                         }),
                         AdjustedData_D21_mean = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           mean(dat_adj[!t21, s])
                         }),
                         AdjustedData_D21_med = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           median(dat_adj[!t21, s])
                         }),
                         AdjustedData_D21_IQR = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           IQR(dat_adj[!t21, s])
                         }),
                         AdjustedData_T21_mean = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           mean(dat_adj[t21, s])
                         }),
                         AdjustedData_T21_med = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           median(dat_adj[t21, s])
                         }),
                         AdjustedData_T21_IQR = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           IQR(dat_adj[t21, s])
                         }),
                         AdjustedData_All_mean = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           mean(dat_adj[, s])
                         }),
                         AdjustedData_All_med = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           median(dat_adj[, s])
                         }),
                         AdjustedData_All_IQR = sapply(names(dat_adj)[-ncol(dat_adj)], function(s){
                           IQR(dat_adj[, s])
                         }))

pvals_table$RawData_FDR_Qval = NA
pvals_table$AdjustedData_FDR_Qval = NA
i = which(!pvals_table$Compound %in% c('KYNTRP_ratio', 'QAPA_ratio'))
pvals_table$RawData_FDR_Qval[i] = p.adjust(pvals_table$RawData_Pval[i], 'fdr')
pvals_table$AdjustedData_FDR_Qval[i] = p.adjust(pvals_table$AdjustedData_Pval[i], 'fdr')
dat$Karyotype = dat_adj$Karyotype = NULL
pvals_table2 = pvals_table[order(pvals_table$AdjustedData_Pval, decreasing = F),
                          c(1,15,16,2:14)]

write.table(pvals_table2, 'Results/Supplementary_Files/Supplementary_Data_10.csv',
              sep = ',', row.names = F)


# Make all plots using adjusted data
fitting = 'withAgeSexAdjustment'
dat_to_plot = dat_adj
fitting_label = 'with model fit for age and sex'
qval_column = 'AdjustedData_FDR_Qval'
direction_column = 'AdjustedData_T21direction'

# Plot boxplots of absolute levels and log2 absolute levels
compounds = c('KYN', 'L_formyl_kynurenine', 'KYNTRP', 'TRP')
for (compound in compounds){
  compound_column = grep(paste0('^', compound, '_ratio'), names(dat_to_plot), value = T)
  if (compound == 'KYNTRP'){
    qval = pvals_table[pvals_table$Compound == compound_column, 'AdjustedData_Pval']
  } else{
    qval = pvals_table[pvals_table$Compound == compound_column, qval_column]
  }
  units = ifelse(compound == 'TRP', 'uM\n(calculated against heavy TRP standard)',
                 'uM\n(calculated against heavy KYN standard)')
  direction = as.character(pvals_table[pvals_table$Compound == compound_column, direction_column])
  qval_sub = paste0('FDR-adjusted pval = ', round(qval, 10))
  plot_main = paste0(compound, ' absolute quant in plasma\n', fitting_label)
    
  # Y axis labels
  if (compound == 'KYNTRP'){
    plot_ylab = 'KYN (uM) / TRP (uM)'
  } else{
    plot_ylab = paste0(compound, ' ', units)
  }
    
  # Plot the absolute levels on normal scale
  df = data.frame(Karyotype = sample_annots[row.names(dat_to_plot), 'Karyotype'], 
                  dat_to_plot[,compound_column])
  names(df)[2] = compound
  pdf(paste0('Results/Figure_1/Fig_1F_boxplot_', compound, '.pdf'), 
      height = 6, width = 6)
  par(mar = c(5,5,5,2))
  plot_t21_d21_boxplot(df, compound, qval, direction, log = F,
                       plot_main = plot_main,
                       plot_ylab = plot_ylab,
                       plot_sub = qval_sub)
  dev.off()
}

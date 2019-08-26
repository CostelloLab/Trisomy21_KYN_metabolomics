# -----------------------------------------------------------------------------
# Script 05: Human Plasma Metabolites & MSD Cytokines
# Author: Rani Powers
#
# Performs differential analysis on human MSD data. Correlates human plasma 
# metabolomics data with MSD cytokine measurements. Saves heatmaps and 
# scatterplots in the Results/Figure_2/ folder. Outputs supplementary data to 
# Results/Supplementary_Files/.
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
source('R/heatmap.3.R')
FDR_CUTOFF = .05

# Read human MSD data
msd_data = read.csv('Data/Human_MSD.csv',
                   stringsAsFactors = F,
                   row.names = 1)

# Read human sample annotations
sample_data = read.table('Data/All_sample_annotations.txt',
                         stringsAsFactors = F, sep = '\t', header = T,
                         row.names = 1)
sample_data = sample_data[names(msd_data),]

# Read human plasma data and subset to the individuals who also had MSD
met_data = exprs(readRDS('Data/Human_plasma_metabolomics_eset_filtered.rds'))[,names(msd_data)]

# Fit model on metabolomics data
design = model.matrix(~0+Age+as.factor(Sex), sample_data[,c('Age', 'Sex')])
met_data2 = removeBatchEffect(met_data, covariates = design)
  
# Fit model on MSD data
msd_data2 = removeBatchEffect(msd_data, covariates = design)

# -----------------------

# Correlate all metabolites and cytokines
ALL_DATA = rbind(met_data, msd_data)
data_for_correlation = as.data.frame(t(ALL_DATA))
data_for_correlation$Kyn_Tryp_ratio = data_for_correlation$C00328 - data_for_correlation$C00078
t21 = row.names(sample_data[sample_data$Karyotype == 'T21',])
d21 = row.names(sample_data[sample_data$Karyotype == 'D21',])

# Correlation results
all_correlation_results = data.frame(Measurement1 = 'A', Measurement2 = 'A',
                                     All_Spearman_Corr = 0, All_Spearman_Pval = 0,
                                     T21_Spearman_Corr = 0, T21_Spearman_Pval = 0,
                                     D21_Spearman_Corr = 0, D21_Spearman_Pval = 0)
for (measurement1 in c('Kyn_Tryp_ratio')){
  for (measurement2 in names(data_for_correlation)[92:ncol(data_for_correlation)]){
    cat(measurement2, '\n')
    all_corr = cor.test(data_for_correlation[,measurement1],
                        data_for_correlation[,measurement2],
                        method = 'spearman', exact = F,
                        use = 'pairwise.complete.obs')
    t21_corr = cor.test(data_for_correlation[t21,measurement1],
                        data_for_correlation[t21,measurement2],
                        method = 'spearman', exact = F,
                        use = 'pairwise.complete.obs')
    d21_corr = cor.test(data_for_correlation[d21,measurement1],
                        data_for_correlation[d21,measurement2],
                        method = 'spearman', exact = F,
                        use = 'pairwise.complete.obs')
    
    all_correlation_results = rbind(all_correlation_results,
                                    data.frame(Measurement1 = measurement1, 
                                               Measurement2 = measurement2,
                                               All_Spearman_Corr = all_corr$estimate, 
                                               All_Spearman_Pval = all_corr$p.value,
                                               T21_Spearman_Corr = t21_corr$estimate, 
                                               T21_Spearman_Pval = t21_corr$p.value,
                                               D21_Spearman_Corr = d21_corr$estimate, 
                                               D21_Spearman_Pval = d21_corr$p.value))
  }
}

all_correlation_results = all_correlation_results[-1,]
all_correlation_results$All_Spearman_Qval = p.adjust(all_correlation_results$All_Spearman_Pval, method = 'fdr')
all_correlation_results$T21_Spearman_Qval = p.adjust(all_correlation_results$T21_Spearman_Pval, method = 'fdr')
all_correlation_results$D21_Spearman_Qval = p.adjust(all_correlation_results$D21_Spearman_Pval, method = 'fdr')

all_correlation_results = all_correlation_results[,c(1:4,9,5:6,10,7:8,11)]
all_correlation_results = all_correlation_results[order(all_correlation_results$All_Spearman_Qval, decreasing = F),]

write.table(all_correlation_results, 
            paste0('Results/Supplementary_Files/Supp_Data_5_MSD_correlation_pvalues.csv'), 
            sep = ',', row.names = F)

sig_cytokines = as.character(all_correlation_results[all_correlation_results$All_Spearman_Qval < .1 |
                                                       all_correlation_results$T21_Spearman_Qval < .1 |
                                                       all_correlation_results$D21_Spearman_Qval < .1, 'Measurement2'])

for (cytokine in c('TNF-A_Crnic_Proinflammatory', 'IP-10_Crnic_Chemokine',
                   'IL-10_Crnic_Proinflammatory', 'IL-29_Crnic_Uplex')){
  
  corr_s = cor.test(data_for_correlation$C00328, data_for_correlation[,cytokine],
                    method = 'spearman', exact = F,
                    use = 'pairwise.complete.obs')
  qval = all_correlation_results[all_correlation_results$Measurement1 == 'C00328' &
                                   all_correlation_results$Measurement2 == cytokine, 
                                 'All_Spearman_Qval']
  
  pdf(paste0('Results/Figure_2/Fig_2c_scatterplot_',
             strsplit(cytokine, '_')[[1]][1], 
             '_', units, '.pdf'), height = 5, width = 5)
  plot(data_for_correlation$C00328, data_for_correlation[,cytokine], 
       pch = 21, 
       xlab = paste0('kynurenine ', units), 
       ylab = paste0(cytokine, ' ', units),
       main = paste0('Kynurenine and ', cytokine, '\ncor = ',
                     round(corr_s$estimate, 3), ', all samples adj p-val = ', 
                     round(qval, 2)),
       bg = ifelse(row.names(data_for_correlation) %in% t21, cols[1], cols[10]))
  dev.off()
}

# Plot heatmaps of MSD data (z-scored by D21 only data)
msd_data = msd_data2 # use adjusted?
d21_msd = msd_data[,sample_data$Karyotype == 'D21']
mu = apply(d21_msd, 1, mean, na.rm = T)
sigma = apply(d21_msd, 1, sd, na.rm = T)
for (i in 1:nrow(d21_msd)){
  d21_msd[i,] = as.numeric(d21_msd[i,] - mu[i])/sigma[i]
}

t21_msd = msd_data[,sample_data$Karyotype == 'T21']
for (i in 1:nrow(t21_msd)){
  t21_msd[i,] = as.numeric(t21_msd[i,] - mu[i])/sigma[i]
}

# KS test for all cytokines, D21 vs T21
d21_msd = as.matrix(d21_msd)
t21_msd = as.matrix(t21_msd)
ks_test_results = data.frame(Cytokine = row.names(d21_msd),
                             D = sapply(1:nrow(d21_msd), function(i){
                               ks.test(d21_msd[i,], t21_msd[i,], exact = F)$statistic
                             }),
                             P_value = sapply(1:nrow(d21_msd), function(i){
                               ks.test(d21_msd[i,], t21_msd[i,], exact = F)$p.value
                             }),
                             Log2MedianFC = sapply(1:nrow(d21_msd), function(i){
                               median(t21_msd[i,]) - median(d21_msd[i,])
                             }))
ks_test_results$Adj_P_value = p.adjust(ks_test_results$P_value, method = 'fdr')
ks_test_results = ks_test_results[order(ks_test_results$Adj_P_value, decreasing = F),]

write.table(ks_test_results, 
            paste0('Results/Supplementary_Files/Supplementary_Data_13.csv'), 
            sep = ',', row.names = F)




sig_cytokines2 = as.character(ks_test_results[ks_test_results$Adj_P_value < FDR_CUTOFF, 'Cytokine'])

# Plot D21 z-scores heatmap for significantly different metabolites
heatmap_title = 'Sig cytokines T21 vs D21, using age- and sex-adjusted log2 intensity MSD values'

pdf('Results/Figure_2/Fig_2b_heatmap_D21.pdf',
    height = 6, width = 10)
heatmap.3(d21_msd[sig_cytokines2,],
          symkey = T,
          breaks = zscore_breaks,
          Rowv = F, Colv = T,
          dendrogram = 'column',
          main = paste0('D21 ', heatmap_title),
          col = heatmap_palette)
dev.off()
pdf('Results/Figure_2/Fig_2b_heatmap_T21.pdf',
    height = 6, width = 10)
heatmap.3(t21_msd[sig_cytokines2,],
          symkey = T,
          breaks = zscore_breaks,
          Rowv = F, Colv = T,
          dendrogram = 'column',
          main = paste0('T21 ', heatmap_title),
          col = heatmap_palette)
dev.off()

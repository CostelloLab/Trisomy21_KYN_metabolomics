# -----------------------------------------------------------------------------
# Script 05: T21 and D21 Fibroblast Metabolomics Timecourse
# Author: Rani Powers
# Last updated: August 1, 2018
#
# Saves line plots in the Results/Figure_4/ folder.
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')

# Data for 3 experiments (called 4, 5 and 6) cells
dat = read.csv('Data/Fibroblast_flux_metabolomics.csv',
               stringsAsFactors = F)
names(dat) = gsub('\\.', '_', names(dat))
names(dat) = gsub('X5', '5', names(dat))
names(dat) = gsub('X3', '3', names(dat))

# Remove batch effects between the 3 replicate experiments
dat[,7:ncol(dat)] = log2(dat[,7:ncol(dat)]+1)
dat_rbe = cbind(dat[,1:7],
                t(removeBatchEffect(t(dat[8:ncol(dat)]), 
                                    batch = factor(dat$Experiment))))

# Add raw ratios
dat_rbe$Kyn_Tryp_Ratio = dat_rbe$L_Kynurenine - dat_rbe$Tryptophan
dat_rbe$Nic_Tryp_Ratio = dat_rbe$Nicotinamide - dat_rbe$Tryptophan

# Plot colors
t21_IFNminus_col = '#4A5059'
d21_IFNminus_col = '#4A5059'
t21_IFNplus_col = cols[1]
d21_IFNplus_col = cols[10]
IFNminus_lty = 1
IFNplus_lty = 1
cell_lines_shapes = c(8, 9, 15, 17, 7, 19)
names(cell_lines_shapes) = unique(dat_rbe$Cell_Line)
sample_group_labels = c('0hr', '1hr', '1hr', '6hr',
                        '6hr', '24hr', '24hr')
treatment_labels = c('IFN-', 'IFN-', 'IFN+', 'IFN-', 'IFN+', 'IFN-', 'IFN+')

# Functions
sd_all = function(x, df){
  out = c()
  timepoints = unique(df$Timepoint)[order(unique(df$Timepoint), decreasing = F)]
  for (timepoint in timepoints){
    out = c(out, sd(df[df$Timepoint == timepoint, x], na.rm = T))
  }
  return(out)
}
agg = function(x, df){
  # Collapse the Expt 4, 5, 6 replicates for each cell line
  # We end up with n = 6 means (from the 6 cell lines) per treatment and timepoint
  out = aggregate(get(x) ~ Timepoint + Treatment + Cell_Line + Karyotype, data = df, FUN = mean)
  names(out)[ncol(out)] = x
  # Duplicate the 0hr value for IFN+
  zero_hr = out[out$Timepoint == 0,]
  zero_hr$Treatment = 'IFN+'
  return(rbind(out, zero_hr))
}
pval_all_treatment = function(x, df){
  # P-value comparing IFN- to IFN+ within the same karyotype group
  pvals = c()
  timepoints = timepoints = unique(df$Timepoint)[order(unique(df$Timepoint), decreasing = F)]
  for (timepoint in timepoints){
    pvals = c(pvals, 
              t.test(df[df$Timepoint == timepoint & df$Treatment == 'IFN-', x],
                     df[df$Timepoint == timepoint & df$Treatment == 'IFN+', x])$p.value)
  }
  return(pvals)
}
convert_pvals = function(p){
  new_p = ifelse(p < .005, '***',
                 ifelse(p < .01, '**',
                        ifelse(p < .05, '*', '')))
  return(new_p)
}

# ----------------------------- PLOT CELLS ------------------------------------
cells = dat_rbe[dat_rbe$Sample_Type == 'Cells',]
cells = cells[, !apply(cells, 2, function(x) all(is.na(x)))]
cells = cbind(ID = paste(cells$Sample_Type, cells$Timepoint, cells$Treatment, sep = '_'), 
              cells)
cells$ID = factor(cells$ID, 
                  levels = c('Cells_0_IFN-', 'Cells_1_IFN-', 'Cells_1_IFN+', 'Cells_6_IFN-',
                             'Cells_6_IFN+', 'Cells_24_IFN-', 'Cells_24_IFN+'))

ALL_PVALS <- data.frame(Compound = 'A', Sample_Type = 'A', Comparison = 'A', 
                        Karyotype = 'A', Timepoint = 0, FC = 0, Pval = 0)

# Cells - lineplot
for (compound in c('Tryptophan', 'L_Kynurenine', 'Kyn_Tryp_Ratio')){
  
  # CELLS - For each experiment, collapse the 3 experiments at each timepoint by mean
  temp = cells
  zero_hr = temp[temp$Timepoint == 0,]
  zero_hr$Treatment = 'IFN+'
  temp = rbind(temp, zero_hr)
  temp = temp[order(temp$Timepoint, decreasing = F),]
  
  T21_sds_IFNminus = sd_all(compound, temp[temp$Treatment == 'IFN-' & temp$Karyotype == 'T21',]) # 4 SDs = 1 per timepoint
  T21_sds_IFNplus = sd_all(compound, temp[temp$Treatment == 'IFN+' & temp$Karyotype == 'T21',])
  T21_fc = median(temp[temp$Treatment == 'IFN-' & temp$Karyotype == 'T21', compound], na.rm = T) -
    median(temp[temp$Treatment == 'IFN+' & temp$Karyotype == 'T21', compound], na.rm = T)
  
  D21_sds_IFNminus = sd_all(compound, temp[temp$Treatment == 'IFN-' & temp$Karyotype == 'D21',])
  D21_sds_IFNplus = sd_all(compound, temp[temp$Treatment == 'IFN+' & temp$Karyotype == 'D21',])
  D21_fc = median(temp[temp$Treatment == 'IFN-' & temp$Karyotype == 'D21', compound], na.rm = T) -
    median(temp[temp$Treatment == 'IFN+' & temp$Karyotype == 'D21', compound], na.rm = T)
  
  # P-value for IFN+ vs IFN-
  T21_pvals = convert_pvals(pval_all_treatment(compound, temp[temp$Karyotype == 'T21',])) # 4 pvals = 1 per timepoint
  D21_pvals = convert_pvals(pval_all_treatment(compound, temp[temp$Karyotype == 'D21',]))
  
  # Track all p-values
  ALL_PVALS = rbind(ALL_PVALS,
                    data.frame(Compound = compound, Sample_Type = 'Cells', Comparison = 'IFN+ vs IFN-',
                               Karyotype = 'T21', Timepoint = c(0,1,6,24),
                               FC = T21_fc, 
                               Pval = pval_all_treatment(compound, temp[temp$Karyotype == 'T21',])),
                    data.frame(Compound = compound, Sample_Type = 'Cells', Comparison = 'IFN+ vs IFN-',
                               Karyotype = 'D21', Timepoint = c(0,1,6,24),
                               FC= D21_fc, 
                               Pval = pval_all_treatment(compound, temp[temp$Karyotype == 'D21',])))
  
  # Plot the mean of the 3 cell lines per timepoint
  means = aggregate(get(compound) ~ Timepoint + Treatment + Karyotype, data = temp, FUN = mean)
  names(means)[ncol(means)] = compound
  means = means[order(means$Timepoint, decreasing = F),]
  
  # Configure axis limits
  max_sd = max(T21_sds_IFNminus, T21_sds_IFNplus, D21_sds_IFNminus, D21_sds_IFNplus, na.rm = T)
  min_y = min(means[,ncol(means)], na.rm = T) - max_sd
  max_y = max(means[,ncol(means)], na.rm = T) + max_sd
  
  pdf(paste0('Results/Figure_4/Fig4_', compound, '_cells.pdf'), 
      height = 5, width = 10)
  par(mar = c(5, 6, 4, 2), mfrow = c(1,2))
  
  # D21 plot
  plot(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'D21', compound],
       main = paste0(compound, ' in D21 cells'),
       pch = 21, bg = d21_IFNminus_col, las = 1, 
       xlab = 'Timepoint', xaxt = 'n',
       ylab = '', lwd = 2, lty = IFNminus_lty,
       type = 'b', col = d21_IFNminus_col,
       ylim = c(min_y, max_y))
  if (grepl('[Rr]atio', compound)){
    mtext(side = 2, line = 4, paste0(compound, ' (3 cell lines, 3 reps per cell line +/- SD)'))
    abline(h = 0, col = 'lightgrey')
  } else {
    mtext(side = 2, line = 4, paste0('Log2 ', compound, ' intensity (3 cell lines, 3 reps per cell line +/- SD)'))
  }
  axis(side = 1, at = 1:4, labels = c('0 hr', '1 hr', '6 hrs', '24 hrs'))
  mtext(side = 1, at = 1:4, D21_pvals, line = 2)
  arrows(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'D21', compound] - D21_sds_IFNminus, 
         1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'D21', compound] + D21_sds_IFNminus, 
         length=0.05, angle=90, code=3, col = d21_IFNminus_col, lwd = 2)
  points(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'D21', compound],
         pch = 21, bg = d21_IFNplus_col, lwd = 2,
         type = 'b', col = d21_IFNplus_col, lty = IFNplus_lty)
  arrows(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'D21', compound] - D21_sds_IFNplus, 
         1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'D21', compound] + D21_sds_IFNplus, 
         length=0.05, angle=90, code=3, col = d21_IFNplus_col, lwd = 2)
  legend('bottomleft', pch = 19, col = c(d21_IFNminus_col, d21_IFNplus_col),
         legend = c('D21, IFN-', 'D21, IFN+'))
  
  # T21 plot
  plot(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'T21', compound],
       main = paste0(compound, ' in T21 cells'),
       pch = 21, bg = t21_IFNminus_col, las = 1, 
       xlab = 'Timepoint', xaxt = 'n',
       ylab = '', lwd = 2, lty = IFNminus_lty,
       type = 'b', col = t21_IFNminus_col,
       ylim = c(min_y, max_y))
  if (grepl('[Rr]atio', compound)){
    mtext(side = 2, line = 4, paste0(compound, ' (3 cell lines, 3 reps per cell line +/- SD)'))
    abline(h = 0, col = 'lightgrey')
  } else {
    mtext(side = 2, line = 4, paste0('Log2 ', compound, ' intensity (3 cell lines, 3 reps per cell line +/- SD)'))
  }
  axis(side = 1, at = 1:4, labels = c('0 hr', '1 hr', '6 hrs', '24 hrs'))
  mtext(side = 1, at = 1:4, T21_pvals, line = 2)
  arrows(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'T21', compound] - T21_sds_IFNminus, 
         1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'T21', compound] + T21_sds_IFNminus, 
         length=0.05, angle=90, code=3, col = t21_IFNminus_col, lwd = 2)
  points(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'T21', compound],
         pch = 21, bg = t21_IFNplus_col, lwd = 2,
         type = 'b', col = t21_IFNplus_col, lty = IFNplus_lty)
  arrows(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'T21', compound] - T21_sds_IFNplus, 
         1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'T21', compound] + T21_sds_IFNplus, 
         length=0.05, angle=90, code=3, col = t21_IFNplus_col, lwd = 2)
  legend('bottomleft', pch = 19, col = c(t21_IFNminus_col, t21_IFNplus_col),
         legend = c('T21, IFN-', 'T21, IFN+'))
  dev.off()
}

# ---------------------------------- PLOT SUPS --------------------------------

sups = dat_rbe[dat_rbe$Sample_Type == 'Supernatant',]
sups = sups[, !apply(sups, 2, function(x) all(is.na(x)))]
sups = cbind(ID = paste(sups$Sample_Type, sups$Timepoint, sups$Treatment, sep = '_'), 
             sups)
sups$ID = factor(sups$ID, 
                 levels = c('Supernatant_0_IFN-', 'Supernatant_1_IFN-', 
                            'Supernatant_1_IFN+', 'Supernatant_6_IFN-', 
                            'Supernatant_6_IFN+', 'Supernatant_24_IFN-', 'Supernatant_24_IFN+'))

# Sups - lineplot
for (compound in c('Tryptophan', 'L_Kynurenine', 'Kyn_Tryp_Ratio')){
  
  # SUPS - For each experiment, collapse the 3 experiments at each timepoint by mean
  temp = sups
  zero_hr = temp[temp$Timepoint == 0,]
  zero_hr$Treatment = 'IFN+'
  temp = rbind(temp, zero_hr)
  temp = temp[order(temp$Timepoint, decreasing = F),]
  
  T21_sds_IFNminus = sd_all(compound, temp[temp$Treatment == 'IFN-' & temp$Karyotype == 'T21',]) # 4 SDs = 1 per timepoint
  T21_sds_IFNplus = sd_all(compound, temp[temp$Treatment == 'IFN+' & temp$Karyotype == 'T21',])
  
  D21_sds_IFNminus = sd_all(compound, temp[temp$Treatment == 'IFN-' & temp$Karyotype == 'D21',])
  D21_sds_IFNplus = sd_all(compound, temp[temp$Treatment == 'IFN+' & temp$Karyotype == 'D21',])
  
  # P-value for IFN+ vs IFN-
  T21_pvals = convert_pvals(pval_all_treatment(compound, temp[temp$Karyotype == 'T21',])) # 4 pvals = 1 per timepoint
  D21_pvals = convert_pvals(pval_all_treatment(compound, temp[temp$Karyotype == 'D21',]))
  
  # Track all p-values
  ALL_PVALS = rbind(ALL_PVALS,
                    data.frame(Compound = compound, Sample_Type = 'Supernatant', Comparison = 'IFN+ vs IFN-',
                               Karyotype = 'T21', Timepoint = c(0,1,6,24),
                               FC = 0, 
                               Pval = pval_all_treatment(compound, temp[temp$Karyotype == 'T21',])),
                    data.frame(Compound = compound, Sample_Type = 'Supernatant', Comparison = 'IFN+ vs IFN-',
                               Karyotype = 'D21', Timepoint = c(0,1,6,24),
                               FC = 0, 
                               Pval = pval_all_treatment(compound, temp[temp$Karyotype == 'D21',])))
  
  # Plot the mean of the 3 cell lines per timepoint
  means = aggregate(get(compound) ~ Timepoint + Treatment + Karyotype, data = temp, FUN = mean)
  names(means)[ncol(means)] = compound
  means = means[order(means$Timepoint, decreasing = F),]
  
  # Configure axis limits
  max_sd = max(T21_sds_IFNminus, T21_sds_IFNplus, D21_sds_IFNminus, D21_sds_IFNplus, na.rm = T)
  min_y = min(means[,ncol(means)], na.rm = T) - max_sd
  max_y = max(means[,ncol(means)], na.rm = T) + max_sd
  
  pdf(paste0('Results/Figure_4/Fig4_', compound, '_sups.pdf'), 
      height = 5, width = 10)
  par(mar = c(5, 6, 4, 2), mfrow = c(1,2))
  
  # D21 plot
  plot(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'D21', compound],
       main = paste0(compound, ' in D21 supernatant'),
       pch = 21, bg = d21_IFNminus_col, las = 1, 
       xlab = 'Timepoint', xaxt = 'n',
       ylab = '', lwd = 2, lty = IFNminus_lty,
       type = 'b', col = d21_IFNminus_col,
       ylim = c(min_y, max_y))
  if (grepl('[Rr]atio', compound)){
    mtext(side = 2, line = 4, paste0(compound, ' (3 cell lines, 3 reps per cell line +/- SD)'))
    abline(h = 0, col = 'lightgrey')
  } else {
    mtext(side = 2, line = 4, paste0('Log2 ', compound, ' intensity (3 cell lines, 3 reps per cell line +/- SD)'))
  }
  axis(side = 1, at = 1:4, labels = c('0 hr', '1 hr', '6 hrs', '24 hrs'))
  mtext(side = 1, at = 1:4, D21_pvals, line = 2)
  arrows(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'D21', compound] - D21_sds_IFNminus, 
         1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'D21', compound] + D21_sds_IFNminus, 
         length=0.05, angle=90, code=3, col = d21_IFNminus_col, lwd = 2)
  points(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'D21', compound],
         pch = 21, bg = d21_IFNplus_col, lwd = 2,
         type = 'b', col = d21_IFNplus_col, lty = IFNplus_lty)
  arrows(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'D21', compound] - D21_sds_IFNplus, 
         1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'D21', compound] + D21_sds_IFNplus, 
         length=0.05, angle=90, code=3, col = d21_IFNplus_col, lwd = 2)
  legend('bottomleft', pch = 19, col = c(d21_IFNminus_col, d21_IFNplus_col),
         legend = c('D21, IFN-', 'D21, IFN+'))
  
  # T21 plot
  plot(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'T21', compound],
       main = paste0(compound, ' in T21 supernatant'),
       pch = 21, bg = t21_IFNminus_col, las = 1, 
       xlab = 'Timepoint', xaxt = 'n',
       ylab = '', lwd = 2, lty = IFNminus_lty,
       type = 'b', col = t21_IFNminus_col,
       ylim = c(min_y, max_y))
  if (grepl('[Rr]atio', compound)){
    mtext(side = 2, line = 4, paste0(compound, ' (3 cell lines, 3 reps per cell line +/- SD)'))
    abline(h = 0, col = 'lightgrey')
  } else {
    mtext(side = 2, line = 4, paste0('Log2 ', compound, ' intensity (3 cell lines, 3 reps per cell line +/- SD)'))
  }
  axis(side = 1, at = 1:4, labels = c('0 hr', '1 hr', '6 hrs', '24 hrs'))
  mtext(side = 1, at = 1:4, T21_pvals, line = 2)
  arrows(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'T21', compound] - T21_sds_IFNminus, 
         1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'T21', compound] + T21_sds_IFNminus, 
         length=0.05, angle=90, code=3, col = t21_IFNminus_col, lwd = 2)
  points(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'T21', compound],
         pch = 21, bg = t21_IFNplus_col, lwd = 2,
         type = 'b', col = t21_IFNplus_col, lty = IFNplus_lty)
  arrows(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'T21', compound] - T21_sds_IFNplus, 
         1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'T21', compound] + T21_sds_IFNplus, 
         length=0.05, angle=90, code=3, col = t21_IFNplus_col, lwd = 2)
  legend('bottomleft', pch = 19, col = c(t21_IFNminus_col, t21_IFNplus_col),
         legend = c('T21, IFN-', 'T21, IFN+'))
  dev.off()
}

# ---------------------- PLOT KYN CELLS / TRYP SUPS ------------------------ #
all_dat = dat_rbe
all_dat = all_dat[, !apply(all_dat, 2, function(x) all(is.na(x)))]
all_dat = cbind(ID = paste(all_dat$Sample_Type, all_dat$Timepoint, all_dat$Treatment, sep = '_'), 
                all_dat)
all_dat$ID = factor(all_dat$ID, 
                  levels = c('Cells_0_IFN-', 'Cells_1_IFN-', 'Cells_1_IFN+', 'Cells_6_IFN-',
                             'Cells_6_IFN+', 'Cells_24_IFN-', 'Cells_24_IFN+',
                             'Supernatant_0_IFN-', 'Supernatant_1_IFN-', 'Supernatant_1_IFN+', 'Supernatant_6_IFN-',
                             'Supernatant_6_IFN+', 'Supernatant_24_IFN-', 'Supernatant_24_IFN+'))

temp_cells = all_dat[all_dat$Sample_Type == 'Cells',]
temp_sups = all_dat[all_dat$Sample_Type == 'Supernatant',]
kyn_sup = temp_sups$L_Kynurenine
tryp_cells = temp_cells$Tryptophan
to_plot = cbind(temp_cells, Kyn_sup_tryp_cells = kyn_sup - tryp_cells)

for (compound in c('Kyn_sup_tryp_cells')){
  
  # SUPS - For each experiment, collapse the 3 experiments at each timepoint by mean
  temp = to_plot
  zero_hr = temp[temp$Timepoint == 0,]
  zero_hr$Treatment = 'IFN+'
  temp = rbind(temp, zero_hr)
  temp = temp[order(temp$Timepoint, decreasing = F),]
  
  T21_sds_IFNminus = sd_all(compound, temp[temp$Treatment == 'IFN-' & temp$Karyotype == 'T21',]) # 4 SDs = 1 per timepoint
  T21_sds_IFNplus = sd_all(compound, temp[temp$Treatment == 'IFN+' & temp$Karyotype == 'T21',])
  
  D21_sds_IFNminus = sd_all(compound, temp[temp$Treatment == 'IFN-' & temp$Karyotype == 'D21',])
  D21_sds_IFNplus = sd_all(compound, temp[temp$Treatment == 'IFN+' & temp$Karyotype == 'D21',])
  
  # P-value for IFN+ vs IFN-
  T21_pvals = convert_pvals(pval_all_treatment(compound, temp[temp$Karyotype == 'T21',])) # 4 pvals = 1 per timepoint
  D21_pvals = convert_pvals(pval_all_treatment(compound, temp[temp$Karyotype == 'D21',]))
  
  # Track all p-values
  ALL_PVALS = rbind(ALL_PVALS,
                    data.frame(Compound = compound, Sample_Type = 'Cells/Sups', Comparison = 'IFN+ vs IFN-',
                               Karyotype = 'T21', Timepoint = c(0,1,6,24),
                               FC = 0,
                               Pval = pval_all_treatment(compound, temp[temp$Karyotype == 'T21',])),
                 data.frame(Compound = compound, Sample_Type = 'Cells/Sups', Comparison = 'IFN+ vs IFN-',
                               Karyotype = 'D21', Timepoint = c(0,1,6,24),
                              FC = 0,
                               Pval = pval_all_treatment(compound, temp[temp$Karyotype == 'D21',])))
  
  # Plot the mean of the 3 cell lines per timepoint
  means = aggregate(get(compound) ~ Timepoint + Treatment + Karyotype, data = temp, FUN = mean)
  names(means)[ncol(means)] = compound
  means = means[order(means$Timepoint, decreasing = F),]
  
  # Configure axis limits
  max_sd = max(T21_sds_IFNminus, T21_sds_IFNplus, D21_sds_IFNminus, D21_sds_IFNplus, na.rm = T)
  min_y = min(means[,ncol(means)], na.rm = T) - max_sd
  max_y = max(means[,ncol(means)], na.rm = T) + max_sd
  
  pdf(paste0('Results/Figure_4/Fig4_', compound, '_sups.pdf'), 
      height = 5, width = 10)
  par(mar = c(5, 6, 4, 2), mfrow = c(1,2))
  
  # D21 plot
  plot(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'D21', compound],
       main = paste0(compound, ' in D21 supernatant'),
       pch = 21, bg = d21_IFNminus_col, las = 1, 
       xlab = 'Timepoint', xaxt = 'n',
       ylab = '', lwd = 2, lty = IFNminus_lty,
       type = 'b', col = d21_IFNminus_col,
       ylim = c(min_y, max_y))
  if (grepl('[Rr]atio', compound)){
    mtext(side = 2, line = 4, paste0(compound, ' (3 cell lines, 3 reps per cell line +/- SD)'))
    abline(h = 0, col = 'lightgrey')
  } else {
    mtext(side = 2, line = 4, paste0('Log2 ', compound, ' intensity (3 cell lines, 3 reps per cell line +/- SD)'))
  }
  axis(side = 1, at = 1:4, labels = c('0 hr', '1 hr', '6 hrs', '24 hrs'))
  mtext(side = 1, at = 1:4, D21_pvals, line = 2)
  arrows(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'D21', compound] - D21_sds_IFNminus, 
         1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'D21', compound] + D21_sds_IFNminus, 
         length=0.05, angle=90, code=3, col = d21_IFNminus_col, lwd = 2)
  points(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'D21', compound],
         pch = 21, bg = d21_IFNplus_col, lwd = 2,
         type = 'b', col = d21_IFNplus_col, lty = IFNplus_lty)
  arrows(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'D21', compound] - D21_sds_IFNplus, 
         1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'D21', compound] + D21_sds_IFNplus, 
         length=0.05, angle=90, code=3, col = d21_IFNplus_col, lwd = 2)
  legend('bottomleft', pch = 19, col = c(d21_IFNminus_col, d21_IFNplus_col),
         legend = c('D21, IFN-', 'D21, IFN+'))
  
  # T21 plot
  plot(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'T21', compound],
       main = paste0(compound, ' in T21 supernatant'),
       pch = 21, bg = t21_IFNminus_col, las = 1, 
       xlab = 'Timepoint', xaxt = 'n',
       ylab = '', lwd = 2, lty = IFNminus_lty,
       type = 'b', col = t21_IFNminus_col,
       ylim = c(min_y, max_y))
  if (grepl('[Rr]atio', compound)){
    mtext(side = 2, line = 4, paste0(compound, ' (3 cell lines, 3 reps per cell line +/- SD)'))
    abline(h = 0, col = 'lightgrey')
  } else {
    mtext(side = 2, line = 4, paste0('Log2 ', compound, ' intensity (3 cell lines, 3 reps per cell line +/- SD)'))
  }
  axis(side = 1, at = 1:4, labels = c('0 hr', '1 hr', '6 hrs', '24 hrs'))
  mtext(side = 1, at = 1:4, T21_pvals, line = 2)
  arrows(1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'T21', compound] - T21_sds_IFNminus, 
         1:4, means[means$Treatment == 'IFN-' & means$Karyotype == 'T21', compound] + T21_sds_IFNminus, 
         length=0.05, angle=90, code=3, col = t21_IFNminus_col, lwd = 2)
  points(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'T21', compound],
         pch = 21, bg = t21_IFNplus_col, lwd = 2,
         type = 'b', col = t21_IFNplus_col, lty = IFNplus_lty)
  arrows(1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'T21', compound] - T21_sds_IFNplus, 
         1:4, means[means$Treatment == 'IFN+' & means$Karyotype == 'T21', compound] + T21_sds_IFNplus, 
         length=0.05, angle=90, code=3, col = t21_IFNplus_col, lwd = 2)
  legend('bottomleft', pch = 19, col = c(t21_IFNminus_col, t21_IFNplus_col),
         legend = c('T21, IFN-', 'T21, IFN+'))
  dev.off()
}

ALL_PVALS = ALL_PVALS[-1,]
ALL_PVALS$Adj_Pval = p.adjust(ALL_PVALS$Pval, method = 'fdr')
write.table(ALL_PVALS, 'Results/Figure_4/Fig4_pvals_table.csv',
            sep = ',', row.names = F)

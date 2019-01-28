# -----------------------------------------------------------------------------
# Script 05: Mouse Plasma Kynurenine Levels
# Author: Rani Powers
# Last updated: December 23, 2018
#
# Plots basal kynurenine levels in WT, Dp10, Dp16 and Dp17 mouse strains and
# calculates statistical significance.
#
# Saves plots in the Results/Figure_4/ folder.
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
if (!dir.exists('Results/Figure_4')) { dir.create('Results/Figure_4', recursive = T) }

# Read mouse data
mouse_data = readRDS('Data/mouse_metabolomics_data.rds')

# Log2 transform and remove batch + sex effect
mouse_data_log2 = t(log2(mouse_data[,7] + 1))
mouse_rbe = as.data.frame(t(removeBatchEffect(mouse_data_log2, batch = factor(mouse_data$Batch),
                                          batch2 = factor(mouse_data$Sex))))
names(mouse_rbe)[1] = 'Kynurenine'
row.names(mouse_rbe) = row.names(mouse_data) = mouse_data$Mouse_ID
mouse_rbe$Strain = mouse_data$Strain

# Plot
pdf('Results/Figure_4/Fig_4b.pdf',
    height = 6, width = 6)
boxplot(Kynurenine ~ Strain, data = mouse_rbe, log = 'y', las = 1,
        ylab = 'log2 kynurenine intensity',
        col = 'lightgrey', main = 'Kynurenine intensity\n(after batch + sex adjustment)')
beeswarm(Kynurenine ~ Strain, data = mouse_rbe, pch = 21,
         bg = c(light_blue, cols[9], light_red, cols[4]), add = T)
t.test(Kynurenine ~ Strain, data = mouse_rbe[mouse_rbe$Strain %in% c('WT', 'Dp16'),])
t.test(Kynurenine ~ Strain, data = mouse_rbe[mouse_rbe$Strain %in% c('WT', 'Dp10'),])
t.test(Kynurenine ~ Strain, data = mouse_rbe[mouse_rbe$Strain %in% c('WT', 'Dp17'),])
dev.off()

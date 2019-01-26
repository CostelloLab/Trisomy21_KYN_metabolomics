# -----------------------------------------------------------------------------
# Script 07: Mouse Plasma Kynurenine Levels
# Author: Rani Powers
# Last updated: August 23, 2018
#
# Plots basal kynurenine levels in WT, Dp10, Dp16 and Dp17 mouse strains and
# calculates statistical significance. Plots kynurenine levels in Dp16 D21 and
# T21 mice with and without poly(I:C) treatment.
#
# Saves plots in the Results/Figure_5/ folder.
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
if (!dir.exists('Results/Figure_5')) { dir.create('Results/Figure_5', recursive = T) }

# Read mouse data
mouse_data = read.csv('Data/Mouse_plasma_metabolomics.csv',
                      stringsAsFactors = F, header = T)
row.names(mouse_data) = apply(mouse_data, 1, function(x){
  paste(x[3], x[7], sep = '_')
})

# ------------------- Panel A - kynurenine in all strains ------------------- #
A_data = mouse_data[mouse_data$Plot == 'A', ]
A_data = cbind(Group_Label = ifelse(A_data$Genotype == 'D21', 'D21', A_data$Strain),
               A_data)

# Remove batch + sex effects
A_data_log2 = t(log2(A_data[,9] + 1))
A_rbe = as.data.frame(t(removeBatchEffect(A_data_log2, batch = factor(A_data$Batch),
                                          batch2 = factor(A_data$Sex))))
names(A_rbe)[1] = 'Kynurenine'
row.names(A_rbe) = row.names(A_data)
A_rbe$Group_Label = A_data[row.names(A_rbe), 'Group_Label']

# Plot
pdf('Results/Figure_5/Fig5A.pdf',
    height = 6, width = 6)
boxplot(Kynurenine ~ Group_Label, data = A_rbe, log = 'y', las = 1,
        ylab = 'log2 kynurenine intensity',
        col = 'lightgrey', main = 'Kynurenine intensity\n(after batch + sex adjustment)')
beeswarm(Kynurenine ~ Group_Label, data = A_rbe, pch = 21,
         bg = c(light_blue, cols[9], light_red, cols[4]), add = T)
t.test(Kynurenine ~ Group_Label, data = A_rbe[A_rbe$Group_Label %in% c('D21', 'Dp16'),])
t.test(Kynurenine ~ Group_Label, data = A_rbe[A_rbe$Group_Label %in% c('D21', 'Dp10'),])
t.test(Kynurenine ~ Group_Label, data = A_rbe[A_rbe$Group_Label %in% c('D21', 'Dp17'),])
dev.off()

# -------------------- Panel B - PolyI:C in Dp16 and WT --------------------- #
B_data = mouse_data[mouse_data$Plot == 'B', ]
B_data$Strain = factor(B_data$Strain)

# Remove batch effects with batch + sex (sex = batch so no explicit adjustment)
table(B_data$Sex, B_data$Batch)
B_data_log2 = t(log2(B_data[,8] + 1))
B_rbe = as.data.frame(t(removeBatchEffect(B_data_log2,
                                          batch = B_data$Batch)))
names(B_rbe)[1] = 'Kynurenine'
row.names(B_rbe) = row.names(B_data)

# Plot
B_rbe$Group_Label = sapply(row.names(B_rbe), function(s){
  paste(B_data[s, 'Genotype'], B_data[s, 'Treatment'], sep = '_')})
B_rbe$Batch = B_data[row.names(B_rbe), 'Batch']
B_rbe$Sex = B_data[row.names(B_rbe), 'Sex']

pdf('Results/Figure_5/Fig5B.pdf',
    height = 6, width = 6)
boxplot(list(D21_Sham = B_rbe[B_rbe$Group_Label == 'D21_Sham', 'Kynurenine'],
             D21_PolyIC = B_rbe[B_rbe$Group_Label == 'D21_Poly I:C', 'Kynurenine'],
             T21_Sham = B_rbe[B_rbe$Group_Label == 'T21_Sham', 'Kynurenine'],
             T21_PolyIC = B_rbe[B_rbe$Group_Label == 'T21_Poly I:C', 'Kynurenine']), 
        log = 'y', las = 1,
        ylab = 'log2 kynurenine intensity',
        col = 'lightgrey', main = 'Kynurenine intensity\n(Dp16 only, after batch/sex adjustment)')
beeswarm(list(D21_Sham = B_rbe[B_rbe$Group_Label == 'D21_Sham', 'Kynurenine'],
              D21_PolyIC = B_rbe[B_rbe$Group_Label == 'D21_Poly I:C', 'Kynurenine'],
              T21_Sham = B_rbe[B_rbe$Group_Label == 'T21_Sham', 'Kynurenine'],
              T21_PolyIC = B_rbe[B_rbe$Group_Label == 'T21_Poly I:C', 'Kynurenine']),
         pch = 21, 
         bg = c(light_blue, dark_blue, light_red, dark_red), add = T)
d21_p = t.test(B_rbe[B_rbe$Group_Label == 'D21_Sham', 'Kynurenine'],
               B_rbe[B_rbe$Group_Label == 'D21_Poly I:C', 'Kynurenine'])$p.value
t21_p = t.test(B_rbe[B_rbe$Group_Label == 'T21_Sham', 'Kynurenine'],
               B_rbe[B_rbe$Group_Label == 'T21_Poly I:C', 'Kynurenine'])$p.value
poly_p = t.test(B_rbe[B_rbe$Group_Label == 'D21_Poly I:C', 'Kynurenine'],
                B_rbe[B_rbe$Group_Label == 'T21_Poly I:C', 'Kynurenine'])$p.value
dev.off()
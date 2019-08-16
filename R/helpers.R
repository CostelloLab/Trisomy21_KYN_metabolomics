# -------------------------------- PACKAGES --------------------------------- #

library(plyr)
library(affy)
library(limma)
library(Biobase)
library(preprocessCore)
library(impute)
library(RColorBrewer)
library(beeswarm)
library(tidyr)
library(dplyr)

# --------------------------------- COLORS ---------------------------------- #

cols = brewer.pal(11, 'Spectral')

# T21 and D21 colors
d21_col = cols[10]
t21_col = cols[1]

# Colors for metabolites that increased or decreased in expression
up_col = cols[1]
down_col = cols[10]

# Boxplot colors
dark_blue = cols[10]
light_blue = '#87c7ee'
dark_red = cols[1]
light_red = '#d48daa'

# Heatmap palette
heatmap_palette = colorRampPalette(c(rep(cols[10], 3), 'black', rep(cols[1], 3)))(n = 100)
spearman_breaks = seq(-1, 1, length=101)
zscore_breaks = seq(-4, 4, length=101)

# --------------------------------- PLOTS ----------------------------------- #

# PCA plotting function
plot_batch_pca = function(pca_summary, main = '', legend_side = 'topleft',
                          bad = c()){
  plot(pca_summary$x[,1], pca_summary$x[,2],
       main = main, las = 1,
       xlab = paste0('PC1 (',
                     round(pca_summary$importance[2,1],2)*100, 
                     '%)'),
       ylab = paste0('PC2 (',
                     round(pca_summary$importance[2,2],2)*100, 
                     '%)'),
       bg = ifelse(sample_data$Barcode %in% bad, 'red', 
                   ifelse(sample_data$Batch == 1, cols[11], cols[8])),
       pch = ifelse(sample_data$Karyotype == 'T21', 24, 
                    ifelse(sample_data$Karyotype == 'D21', 21, 1)))
  legend(legend_side, 
         legend = c('Nexus T21',
                    'Nexus D21',
                    'HTP T21',
                    'HTP D21'),
         cex = .5,
         pt.cex = 1,
         pch = c(24, 21, 24, 21),
         pt.bg = c(rep(cols[11], 2), rep(cols[8], 2)))
}

# Volcano plotting function
plot_volcano = function(x, y, main = '', xlab = 'log2 FC', ylab = '-log10 adj p',
                        x_cutoffs = log2(c(.8, 1.2)), las = 1,
                        y_cutoff = -log10(.01), plot_fc_lines = T,
                        up_col = '#9E0142', down_col = '#3288BD',
                        to_color = NULL,    # indexes of points to color - will be automatically colored up/down
                        pt_labels = NULL,   # all labels - will automatically only write the sig ones
                        xlim = NULL, ylim = NULL){
  sig_up = which(x > x_cutoffs[2] & y > y_cutoff) 
  sig_down = which(x < x_cutoffs[1] & y > y_cutoff)
  if (is.null(to_color)){
    to_color = c(sig_up, sig_down)
  }
  pt_colors = rep('grey', length(x))
  pt_colors[to_color] = ifelse(x[to_color] > 0, up_col,
                               ifelse(x[to_color] < 0, down_col, 'grey'))
  plot(x, y,
       main = main,
       xlab = xlab,
       ylab = ylab,
       xlim = xlim,
       ylim = ylim,
       las = las,
       pch = 21,
       bg = pt_colors)
  if (plot_fc_lines){
    abline(v = c(x_cutoffs[1], x_cutoffs[2]), 
           col = c(down_col, up_col))
  }
  abline(v = 0, 
         lty = 2,
         col = 'grey')
  abline(h = y_cutoff, col = 'grey')
  if (!is.null(pt_labels)){
    if (any(sig_up)){
      text(x = x[sig_up],
           y = y[sig_up],
           labels = pt_labels[sig_up],
           pos = 4,
           cex = .5)
    }
    if (any(sig_down)){
      text(x = x[sig_down],
           y = y[sig_down],
           labels = pt_labels[sig_down],
           pos = 2,
           cex = .5)
    }
  }
}

# Boxplot function
plot_boxplot = function(metabolite, df, feature_data, ylim = NULL, ylab = '',
                        point_borders = 'none',  # c('black', 'color', 'none')
                        col1 = '#013c7c', col2 = '#9cb7d9',
                        pval_column = 'P.Value', qval_column = 'adj.P.Val'){
  met_name = feature_data[metabolite, 'Compound_Name']
  y_max = max(as.numeric(df[,metabolite]), na.rm = T)
  y_min = min(as.numeric(df[,metabolite]), na.rm = T)
  if (is.null(ylim)){
    ylim = c(y_min, y_max)
  }
  boxplot(get(metabolite)~Karyotype, data = df,
          outline = F,
          main = met_name,
          las = 1,
          whisklty = 0, staplelty = 0,
          ylab = ylab,
          ylim = ylim,
          col = 'lightgrey',
          sub = paste0('unadjusted p = ',
                       round(met_sigtable[met_sigtable$Compound_ID == metabolite, pval_column], 10),
                       '\nFDR-adjusted p = ',
                       round(met_sigtable[met_sigtable$Compound_ID == metabolite, qval_column], 10)))
  if (point_borders == 'black'){
    beeswarm(get(metabolite)~Karyotype, data = df,
             pch = 21,
             method = 'swarm',
             spacing = 1,
             pwbg = ifelse(df$Karyotype == 'T21', col1, col2),
             add = T)
  }
  if (point_borders == 'color') {
    border_col1 = ifelse(col1 == dark_blue, '#25617F',
                         ifelse(col1 == dark_red, '#680130',
                                ifelse(col1 == '#515151', '#353535', 'black')))
    border_col2 = ifelse(col2 == light_blue, '#689bb5',
                         ifelse(col2 == light_red, '#93657a',
                                ifelse(col2 == '#A0A0A0', '#515151', 'black')))
    beeswarm(get(metabolite)~Karyotype, data = df,
             pch = 21,
             method = 'swarm',
             spacing = 1,
             pwbg = ifelse(df$Karyotype == 'T21', col1, col2),
             pwcol = ifelse(df$Karyotype == 'T21', border_col1, border_col2),
             add = T)
  }
  if (point_borders == 'none'){
    beeswarm(get(metabolite)~Karyotype, data = df,
             pch = 19,
             method = 'swarm',
             spacing = 1,
             pwcol = ifelse(df$Karyotype == 'T21', col1, col2),
             add = T)
  }
}

# Single row heatmap
plot_single_row_heatmap = function(data_vector, karyotype, column_label, n,
                                   key_label = 'Spearman Corr'){
  par(cex.main = .7)
  suppressWarnings(
    heatmap.3(cbind(rep(0, length(data_vector)), data_vector),
              main = paste0('Spearman correlation across all\nn = ',
                            n, ' ', karyotype, ' samples, adjusted for age + gender'),
              Rowv = F, Colv = F,
              breaks = spearman_breaks,
              symkey = T,
              symbreaks = T,
              KeyValueName = 'Spearman Corr',
              dendrogram = 'none',
              labCol = c('', column_label),
              cexCol = .7,
              cellnote = cbind(rep(NA, length(data_vector)), 
                               sapply(data_vector, round, 2)),
              notecol = 'black',
              margins = c(5,8),
              col = heatmap_palette)
  )
}

plot_msd_met_scatterplot = function(met1, met2, df, t21, d21, plot_fit = F){
  xmin = min(c(df[t21, met1], df[d21, met1]), na.rm = T)
  xmax = max(c(df[t21, met1], df[d21, met1]), na.rm = T)
  ymin = min(c(df[t21, met2], df[d21, met2]), na.rm = T)
  ymax = max(c(df[t21, met2], df[d21, met2]), na.rm = T)
  t21_corr_s = cor(df[t21, met1], df[t21, met2], method = 'spearman',
                   use = 'pairwise.complete.obs')
  plot(df[t21, met1], 
       df[t21, met2],
       pch = 19, col = cols[2], 
       main = paste0(met1, ' and ', met2, 
                     ' in T21 and D21 samples\nresiduals after adjustment for age + gender'),
       las = 1,
       xlab = paste0(met1, ' (residuals)'),
       ylab = paste0(met2, ' (residuals)'),
       xlim = c(xmin, xmax), 
       ylim = c(ymin, ymax))
  points(df[d21, met1], 
         df[d21, met2],
         pch = 19, col = 'black')
  mtext(text = paste0('T21 Spearman: ', round(t21_corr_s, 2)),
        side = 1, line = 5)
  if (plot_fit){
    #abline(lm(get(met2) ~ get(met1), df[t21,]), col = cols[2])
    #abline(lm(get(met2) ~ get(met1), df[d21,]), col = 'black')
    #abline(lm(get(met2) ~ get(met1), df), col = 'grey', lty = 2)
    lw_t21 = loess(get(met2) ~ get(met1), df[t21,], col = cols[2], span = 5)
    i = order(df[t21,met1])
    lines(df[t21,met1][i], lw_t21$fitted[i], col = cols[2])
    lw_d21 = loess(get(met2) ~ get(met1), df[d21,], col = 'black', span = 5)
    i = order(df[d21,met1])
    lines(df[d21,met1][i], lw_d21$fitted[i], col = 'black')
    legend('topleft', c('T21', 'D21'),
           pch = c(19, 19), col = c(cols[2], 'black'))
  }
}

# More colors
up_down_cols = list(light_red = '#D48DAA', dark_red = '#9E0142',
                    light_blue = '#72AECC', dark_blue = '#3288BD')
get_pval_cols = function(qval, direction){
  if (qval < .05){
    if (direction == 'Up'){
      plot_cols = c(up_down_cols$light_red, up_down_cols$dark_red)
    } else{
      plot_cols = c(up_down_cols$light_blue, up_down_cols$dark_blue)
    }
  } else {
    plot_cols = 'darkgrey'
  }
  return(plot_cols)
}

# Plots
plot_t21_d21_boxplot <- function(df, compound, qval, direction, log = F, 
                                 zero_log_threshold = 1e-4, exclude_zeros = F,
                                 plot_main = NA, plot_ylab = NA, plot_sub = NA){
  # df is a 2-column df with column names: Karyotype, <compound>
  if (is.na(plot_main)){
    plot_main = compound
  }
  if (is.na(plot_ylab)){
    plot_ylab = compound
  }
  if (is.na(plot_sub)){
    plot_sub = ''
  }
  
  # Determine whether there are any samples with zero values
  tab = as.data.frame(table(df$Karyotype, df[,compound] == 0))
  d21_zeros = tab[tab$Var1 == 'D21' & tab$Var2 == 'TRUE', 'Freq']
  t21_zeros = tab[tab$Var1 == 'T21' & tab$Var2 == 'TRUE', 'Freq']
  
  # Optionally exclude or threshold the zeros when showing data on log scale
  phrase_end = ')'
  if (log){
    plot_log = 'y'
    
    # Calculate y axis limits
    ymin = ifelse((min(df[,compound], na.rm = T)-(0.5*sd(df[,compound], na.rm = T))) <= 0,
                  zero_log_threshold,
                  (min(df[,compound], na.rm = T)-(0.5*sd(df[,compound], na.rm = T))))
    ylimits = c(ymin, (max(df[,compound], na.rm = T)+sd(df[,compound], na.rm = T)))
    
    if (exclude_zeros){
      df[,compound][df[,compound] == 0] = NA
      phrase_end = ' excluded)'
    } else{
      df[,compound][df[,compound] == 0] = zero_log_threshold
    }
    
  } else{
    plot_log = ''
    ylimits = c((min(df[,compound], na.rm = T)-sd(df[,compound], na.rm = T)), 
                (max(df[,compound], na.rm = T)+sd(df[,compound], na.rm = T)))
  }
  
  boxplot(get(compound) ~ Karyotype, data = df, outline = F, ylim = ylimits,
          las = 1, col = 'lightgrey', log = plot_log,
          main = plot_main, ylab = plot_ylab, sub = plot_sub)
  beeswarm(get(compound) ~ Karyotype, data = df, add = T,
           pch = 16, col = get_pval_cols(qval, direction))
  
  # Add text under the x axis if some samples had zero values
  if (length(d21_zeros) > 0){
    mtext(text = paste0('(n = ', d21_zeros,
                        ' zero values', phrase_end), 
          side = 1, at = 1, line = 2, cex = .75)
  }
  if (length(t21_zeros) > 0){
    mtext(text = paste0('(n = ', t21_zeros,
                        ' zero values', phrase_end), 
          side = 1, at = 2, line = 2, cex = .75)
  }
}
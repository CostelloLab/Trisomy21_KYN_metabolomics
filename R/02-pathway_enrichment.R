# -----------------------------------------------------------------------------
# Script 02: Human Plasma Pathway Enrichment Analysis
# Author: Rani Powers
# Last updated: July 31, 2018
#
# Annotates metabolites in the human metabolomics eset to Kegg pathways and
# calculates p-values based on a hypergeometric test.
#
# Saves enrichment plot in the Results/Figure1/ folder.
# -----------------------------------------------------------------------------

# Load libraries, color palettes and plotting functions
source('R/helpers.R')
library(KEGGREST)

# Get the data on the 91 metabolites
eset = readRDS('Data/Human_plasma_metabolomics_filtered.rds')
met_data = exprs(eset)
sample_data = pData(eset)
feature_data = fData(eset)

# -------------------------- PATHWAY ANNOTATION -------------------------------
# Get KEGG IDs, excluding 8 HMDB and other non-KEGG IDs (83 total metabolites)
kegg_ids = feature_data$Compound_ID
kegg_ids = kegg_ids[!grepl('^HMDB|^ac|^CID', kegg_ids)]

# Use KEGG API to look up the pathways for each metabolite
l = list()
p = c()
not_found = c()
for (kegg in kegg_ids){
  cat(paste0('Annotating ', kegg, '\n'))
  query = keggGet(kegg)
  if ('PATHWAY' %in% names(query[[1]])){
    l[kegg] = list(query[[1]]$PATHWAY)
    p = unique(c(p, list(query[[1]]$PATHWAY)))
  }
  else{
    not_found = c(not_found, kegg)
  }
}

# 4 metabolites were not found by keggGet()
not_found

# Total, 12 of 91 do not have pathway annotations
# Of these, 4 are in Kegg but have no pathway annotations, and 8 are not in Kegg
missing = feature_data$Compound_ID[!feature_data$Compound_ID %in% names(l)]
# feature_data[missing,]  # mainly carnitines and gamma glutamyls

# Save pathway names their Kegg pathway IDs
all_annotations = unlist(l)
names(all_annotations) = gsub('C[0-9]+\\.', '', names(all_annotations))
all_annotations = all_annotations[!duplicated(all_annotations)]
pathway_ids = data.frame(Pathway = as.character(all_annotations), 
                         Pathway_ID = names(all_annotations),
                         stringsAsFactors = F)

# Store the total # metabolites in each pathway
pathway_ids$Mets_in_Pathway = sapply(pathway_ids$Pathway, function(x){
  sum(unlist(l) == x)
})
row.names(pathway_ids) = pathway_ids$Pathway_ID

# Make a list for each pathway that contains its metabolites
pathways = pathway_ids$Pathway
mets_in_pathways = list()
for (path in pathways){
  mets_in_pathways[path] = list(NA)
  for (kegg in names(l)){
    temp = l[kegg]
    if (path %in% as.character(unlist(temp))){
      mets_in_pathways[[path]] = na.omit(c(mets_in_pathways[[path]], kegg))
    }
  }
}

# We will only do pathway enrichment based on the metabolites that were in Kegg
all_mets = names(l)         # 79 metabolites total have a pathway annotation
sig_table = read.csv('Results/Figure_1/Fig1_SuppTable.csv',
                     header = T, stringsAsFactors = F)
sig_table = sig_table[sig_table$Compound_ID %in% all_mets,]

# Of the 79, 18 metabolites are differentially expressed
is_met_de = all_mets %in% sig_table[sig_table$adj.P.Val <= .01, 'Compound_ID']
num_sig = sum(is_met_de)
sig_mets = all_mets[is_met_de]

# Filter out pathways that use the exact same set of metabolites
mets_in_pathways = lapply(mets_in_pathways, sort)
sets = list()
for (pathway_name in names(mets_in_pathways)){
  mets = paste0(unlist(mets_in_pathways[pathway_name]), collapse = '')
  if (any(mets == names(sets))){
    sets[which(mets == names(sets))] = list(c(sets[[which(mets == names(sets))]], 
                                              pathway_name))
  } else{
    sets = c(sets, pathway_name)
    names(sets)[length(sets)] = mets
  }
}

# 116 pathways that aren't duplicated, 76 of which have > 4 metabolites in them 
not_duplicated = unlist(sets[lapply(sets, length) == 1])
at_least_5 = names(not_duplicated)[sapply(names(not_duplicated), function(p) { sum(grepl('C', strsplit(p, '')[[1]])) > 4})]

# 20 sets have duplicate pathway annotations, only one of which has > 4 metabolites in it
duplicated_sets = sets[lapply(sets, length) > 1]
dup_at_least_5 = names(duplicated_sets)[sapply(names(duplicated_sets), function(p) { sum(grepl('C', strsplit(p, '')[[1]])) > 4})]
duplicated_sets[dup_at_least_5] # Citrate cycle (TCA cycle)
  
# 52 final pathways will be used for the enrichment calculation
final_pathways = c(not_duplicated[at_least_5], 
                   C00022C00026C00036C00042C00122C00149C00158 = 'Citrate cycle (TCA cycle)')

# -------------------------- PATHWAY ENRICHMENT -------------------------------
ids = c()
names = c()
Q = c()
M = c()
N = c()
pvals = c()
sets_seen = c()
for (path in final_pathways){
  
  q = m = n = 0
  path_id = pathway_ids[pathway_ids$Pathway == path, 'Pathway_ID']
  mets = as.character(unlist(mets_in_pathways[path]))
  
  # m = total number of mets in pathway
  m = length(mets)
  
  # q = number of mets in pathway that are diff expr
  q = length(intersect(mets, sig_mets))
  
  # n = 79 total mets with annotations - m
  n = length(all_mets) - m
  
  # If there were at least 3 mets in the pathway, calculate p-value
  M = c(M, m)
  Q = c(Q, q)
  N = c(N, n)
  ids = c(ids, path_id)
  names = c(names, path)
  pvals = c(pvals, phyper(q, m, n, length(sig_mets), lower.tail = F))
}
hypergeom_results = data.frame(Pathway_ID = ids,
                               Pathway_Name = names,
                               Diff_expr_in_pathway = Q, 
                               Total_mets_in_pathway = M, 
                               Total_mets_not_in_pathway= N, 
                               Total_diff_expr_mets = length(sig_mets),
                               pval = pvals,
                               stringsAsFactors = F)

hypergeom_results$qval = p.adjust(hypergeom_results$pval, method = 'bonferroni')
hypergeom_results = hypergeom_results[order(hypergeom_results$pval, decreasing = F),]
write.table(hypergeom_results,
            'Results/Figure_1/Fig1C.csv', sep = ',', row.names = F)

# Plot pathway enrichment for top 5 pathways
pathways_to_plot = hypergeom_results$Pathway_Name[1:5]

# Rank by fold change
sig_table = sig_table[order(sig_table$logFC, decreasing = F),]

pdf('Results/Figure_1/Fig1C.pdf',
    height = 8, width = 5)
par(mfrow = c(5,1))
for (i in 1:length(pathways_to_plot)){
  id = pathways_to_plot[i]
  mets_in_pathway = mets_in_pathways[[id]]
  ranks = sig_table$Compound_ID %in% mets_in_pathway
  fc = sig_table$logFC[ranks]
  sig = sig_table$adj.P.Val[ranks] <= .01
  plot(which(ranks), rep(1, length(mets_in_pathway)),
       main = id,
       sub = paste0('(unadj p = ', 
                    round(hypergeom_results[hypergeom_results$Pathway_Name == id, 'pval'], 3), 
                    ', adj p = ', 
                    round(hypergeom_results[hypergeom_results$Pathway_Name == id, 'qval'], 3), ')'),
       xlab = 'Metabolite rank by increasing fold change',
       xlim = c(0,92),
       ylab = '', ylim = c(.8,1.2),
       pch = '.', yaxt = 'n',
       col = ifelse(fc < 0 & sig, cols[10], 
                    ifelse(fc > 0 & sig, cols[1], 'grey')))
  for (i in 1:length(mets_in_pathway)){
    abline(v = which(ranks)[i],
           lwd = 4,
           col = ifelse(fc[i] < 0 & sig[i], cols[10], 
                        ifelse(fc[i] > 0 & sig[i], cols[1], 'grey')))
  }
}
dev.off()
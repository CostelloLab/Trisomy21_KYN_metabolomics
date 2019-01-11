---
output:
  html_document: default
  pdf_document: default
---
## Trisomy 21 activates the kynurenine pathway via increased dosage of interferon receptors

Our preprint is available on bioRxiv [here](https://www.biorxiv.org/content/early/2018/08/29/403642).

**Authors:**

Rani K. Powers<sup>1,2,3#</sup>, Kelly D. Sullivan<sup>1,3,4,5#</sup>, Rachel Culp-Hill<sup>6</sup>, Michael P. Ludwig<sup>1,3</sup>, Keith P. Smith<sup>1</sup>, Katherine A. Waugh<sup>1</sup>, Ross Minter<sup>1</sup>, Kathryn D. Tuttle<sup>1</sup>, Hannah C. Lewis<sup>1</sup>, Angela L. Rachubinski<sup>1,5</sup>, Ross E. Granrath<sup>1</sup>, Rebecca B. Wilkerson<sup>6</sup>, Darcy E. Kahn<sup>1</sup>, Molishree Joshi<sup>4</sup>, Angelo D’Alessandro<sup>1,6</sup>, James C. Costello<sup>2,3,&ast;</sup> & Joaquin M. Espinosa<sup>1,3,4,7&ast;</sup>

**Affiliations:**

<sup>1</sup> Linda Crnic Institute for Down Syndrome, University of Colorado Anschutz Medical Campus, Aurora, CO, USA 

<sup>2</sup> Computational Bioscience Program, University of Colorado Anschutz Medical Campus, Aurora, CO, USA 

<sup>3</sup> Department of Pharmacology, University of Colorado Anschutz Medical Campus, Aurora, CO, USA 

<sup>4</sup> Functional Genomics Facility, University of Colorado Anschutz Medical Campus, Aurora, CO, USA

<sup>5</sup> Department of Pediatrics, University of Colorado Anschutz Medical Campus, Aurora, CO, USA

<sup>6</sup> Department of Biochemistry and Molecular Genetics, University of Colorado Anschutz Medical Campus,
Aurora, CO, USA

<sup>7</sup> Department of Molecular, Cellular and Developmental Biology, University of Colorado Boulder, Boulder, CO, USA

<sup>#</sup> These authors contributed equally to this work

**&ast; Corresponding Authors:**

James C. Costello: james.costello@ucdenver.edu

Joaquin M. Espinosa: joaquin.espinosa@ucdenver.edu

---

### Included in this repository

To ensure reproducibility, all data for this project can be found in the Data directory and the R code to generate all figures can be found in the R directory.

#### Data

* **human_plasma_metabolomics.csv** - Cohort 1 and Cohort 2 participant characteristics and raw, LC-MS plasma metabolomics intensity values (relative abundance)
* **human_plasma_metabolomics_eset.rds** - Same data as Human_plasma_metabolomics.csv formatted as an R [eset](https://www.bioconductor.org/packages/3.7/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) object
* **human_plasma_MSD.csv** - Mesoscale discovery assay data from Cohort 2 human plasma samples
* **human_fibroblast_flux.csv** - Targeted metabolic flux data from timecourse experiment on human fibroblasts with IFN-α stimulation
* **human_fibroblast_flux_KO.csv** - Targeted metabolic flux data from timecourse experiment on human fibroblasts with IFNAR knockout and IFN-α stimulation
* **Mouse_plasma_metabolomics.csv** - Targeted metabolomics data on plasma from WT and T21 (Dp10, Dp16 and Dp17 strains) mice

#### Code

To reproduce the figures from our paper, make sure you have all dependencies installed then run each of the numbered R scripts in order. These scripts will create a directory called Results and output plots and/or tables to an appropriate subdirectory named for each figure (e.g. Results/Figure_1). General descriptions of each script are below, and additional information can be found at the beginning of the scripts themselves.

Dependencies:

* plyr (>= 1.8.4), tidyr (>= 0.8.0), dplyr (>= 0.7.4)
* affy (>= 1.54.0), Biobase (>= 2.36.2), preprocessCore (>= 1.38.1)
* impute (>= 1.50.1), limma (>= 3.32.10), glmnet (>= 2.0-13)
* RColorBrewer (>= 1.1-2), beeswarm (>= 0.2.3)
* ROCR (>= 1.0-7), pROC (>= 1.10.0), classInt (>= 0.2-3)

R Scripts:

* **00-model_fitting.R** - pre-processes human metabolomics data and fits linear model with batch, age and sex covariates. Saves density plots for Supplementary Figure 1b-c in the Results/Supplementary_Figures directory.
* **01-differential_abundance.R** - performs differential abundance analysis on human plasma samples. Saves plots for Figure 1A and Figure 1 figure supplement 1B-C in the Results/Figure_1 directory.
* **02-pathway_enrichment.R** - performs hypergeometric testing for pathway enrichment on human plasma samples. Saves plot for Figure 1C and supplemental table of results in the Results/Figure_1 directory.
* **03-boxplots_t21_vs_d21** - plots boxplots of differentially abundant metabolites, tryptophan pathway metabolites and metabolite ratios. Saves plots for Figure 1B, Figure 3, and Figure 3 supplement 1 in the Results/Figure_1 and Results/Figure_3 directories.
* **04-lasso_biomarker.R** - performs univariate and multivariate biomarker analysis to learn a metabolic signature of trisomy 21. Saves plots for Figure 2A-C in the Results/Figure_2 directory.
* **05-fibroblast_timecourse.R** - pre-processes and plots timecourse data from the flux metabolomics experiment on human fibroblasts. Saves plots for Figure 4A-C and Figure 4 figure supplement 1A-D in the Results/Figure_4 directory.
* **06-cytokine_correlations.R** - correlates kynurenine and quinolinic acid in human plasma with Mesoscale Discovery assay data. Saves plots for Figure 6 in the Results/Figure_6 directory.
* **07-mouse_kynurenine.R** - plots kynurenine levels in mouse plasma with and without treatment with poly(I:C). Saves plots for Figure 5 in the Results/Figure_5 directory.


## Trisomy 21 activates the kynurenine pathway via increased dosage of interferon receptors

Our preprint is available on bioRxiv [here](https://www.biorxiv.org/content/early/2018/08/29/403642).

**Authors:**

Rani K. Powers<sup>1,2,3</sup>, Rachel Culp-Hill<sup>4</sup>, Michael P. Ludwig<sup>1,3</sup>, Keith P. Smith<sup>1</sup>, Katherine A. Waugh<sup>1</sup>, Ross Minter<sup>1</sup>, Kathryn D. Tuttle<sup>1</sup>, Hannah C. Lewis<sup>1</sup>, Angela L. Rachubinski<sup>1,5</sup>, Ross E. Granrath<sup>1</sup>, Maria Carmona<sup>6,7</sup>, Rebecca B. Wilkerson<sup>4</sup>, Darcy E. Kahn<sup>1</sup>, Molishree Joshi<sup>8</sup>, Alberto Lleo<sup>6</sup>, Rafael Blesa<sup>6</sup>, Juan Fortea<sup>6,7</sup>, Angelo D’Alessandro<sup>1,4</sup>, James C. Costello<sup>2,3</sup>, Kelly D. Sullivan<sup>1,3,5,8 *</sup>, & Joaquin M. Espinosa<sup>1,3,8,9*</sup>


**Affiliations:**

<sup>1</sup> Linda Crnic Institute for Down Syndrome, University of Colorado Anschutz Medical Campus, Aurora, CO, USA 

<sup>2</sup> Computational Bioscience Program, University of Colorado Anschutz Medical Campus, Aurora, CO, USA 

<sup>3</sup> Department of Pharmacology, University of Colorado Anschutz Medical Campus, Aurora, CO, USA 

<sup>4</sup> Department of Biochemistry and Molecular Genetics, University of Colorado Anschutz Medical Campus,
Aurora, CO, USA

<sup>5</sup> Department of Pediatrics, University of Colorado Anschutz Medical Campus, Aurora, CO, USA

<sup>6</sup> Department of Neurology, Hospital de la Santa Creu i Sant Pau, Biomedical Research Institute Sant Pau, Universitat Autonoma de Barcelona, CIBERNED, Barcelona, Spain

<sup>7</sup> Down Medical Center, Catalan Down Syndrome Foundation, Barcelona, Spain

<sup>8</sup> Functional Genomics Facility, University of Colorado Anschutz Medical Campus, Aurora, CO, USA

<sup>9</sup> Department of Molecular, Cellular and Developmental Biology, University of Colorado Boulder, Boulder, CO, USA

**&ast; Corresponding Authors:**

Kelly D. Sullivan: kelly.d.sullivan@ucdenver.edu

Joaquin M. Espinosa: joaquin.espinosa@ucdenver.edu

---

### Included in this repository

To ensure reproducibility, all data for this project can be found in the Data/ directory and the R code to generate all figures and supplementary data from the manuscript can be found in the R/ directory.

#### Data

* **Human_plasma_metabolomics_eset.rds** - Cohort 1 and Cohort 2 participant characteristics and raw, LC-MS plasma metabolomics intensity values (relative abundance) formatted as an R [eset](https://www.bioconductor.org/packages/3.7/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf)
* **human_MSD_data.rds** - Mesoscale discovery assay data from Cohort 2 human plasma samples
* **mouse_metabolomics_data.rds** - Targeted metabolomics data on plasma from WT and T21 (Dp10, Dp16 and Dp17 strains) mice

#### Code

To reproduce the figures from our paper, make sure you have all dependencies installed then run each of the numbered R scripts in order. These scripts will output plots and/or tables to an appropriate subdirectory named for each figure (e.g. Results/Figure_1). General descriptions of each script are below, and additional information can be found at the beginning of the scripts themselves.

Dependencies:

* plyr (>= 1.8.4), tidyr (>= 0.8.0), dplyr (>= 0.7.4)
* affy (>= 1.54.0), Biobase (>= 2.36.2), preprocessCore (>= 1.38.1)
* impute (>= 1.50.1), limma (>= 3.32.10), glmnet (>= 2.0-13)
* RColorBrewer (>= 1.1-2), beeswarm (>= 0.2.3)
* ROCR (>= 1.0-7), pROC (>= 1.10.0), classInt (>= 0.2-3)

R Scripts:

* **00-model_fitting.R** - Pre-processes human metabolomics data and fits linear model with batch, age and sex covariates. Saves plots and supplementary data.
* **01-differential_abundance.R** - Performs differential abundance analysis on human plasma samples from D21 and T21 individuals. Saves plots and supplementary data.
* **02-boxplots_t21_vs_d21** - Saves boxplots for kynurenine, quinolinic acid, KYN/TRP, and QA/PA ratios in the Results/Figure_1/ folder. Additional boxplots for Supp Figure 1 are saved in the Results/Supplementary_Figures/ folder.
* **03-cytokine_correlations.R** - Correlates Cohort 2 human plasma metabolomics data with MSD cytokine measurements. Saves heatmaps and scatterplots in the Results/Figure_2/ folder. Outputs supplementary data to Results/Supplementary_Files/.
* **04-fibroblast_timecourse.R** - Pre-processes and plots timecourse data from the flux metabolomics experiment on human fibroblasts. Saves plots for Figure 3 in the Results/Figure_3 directory.
* **05-mouse_kynurenine.R** - Plots basal kynurenine levels in WT, Dp10, Dp16 and Dp17 mouse strains and calculates statistical significance. Saves plots in the Results/Figure_4/ folder.


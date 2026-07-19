# Analyses_of_37_composts
R code for pubication Logo et al. 2025 in Applied and Environmental Microbiology 

"Analyses of 37 composts revealed microbial taxa associated with disease suppressiveness"

All analysis were conducted in R (v4.3.1) using RStudio (v2023.06). The R session outputs are seperatly saved as txt files for the different R scripts.

The raw sequencing data are published as an SRA at the NCBI: Bacteria (PRJNA1200637), Fungi (PRJNA1201072)

This respitory contains the following: 

*Figures in publication*
* Figures_Logo_et_al_AEM_2025:  Folder containing all Figures in the main part and supplementary of the publication

*Analysis Bioassays, compost characteristics and bacterial communties of 37 composts*
* Bioassay_normal_distribution.R: R code which helps to automize the statistical analysis of the bioassay data. Output file: Bioassay_normal_distribution_filled.csv
* Figures_analysis_paper1_GTDB: Main R code analysis for the bioassays, compost characteristics and bacterial community analysis, most figures and results were generated with this file
* Figures_analysis_paper_1_GTDB_calc.R:  R code to calculate alpha diversity, beta diversity etc. for bacterial communities
* setup.R:  Set up R enviroment
* data_Figures_analysis_paper1_GTDB: Data needed for analysis

*Analysis of fungal communties*
* Fungi_37_analysis.R: R code for community analysis, sequencing overview, indicator analysis for fungi, results generated for Figure S13, Table S7 additional figures including the fungal communties are also in the combined R script (Fungi_Bacteria_37_analysis), the data needed is stored in data_Fungi_Bacteria_37_analysis
* Fungi_37_calc.R: R code to calculate alpha diverisy, beta diversity ect. for fungal communities
* setup_F24.R: Set up R enviroment*data_Fungi_37_analysis.R
* data_Fungi_37_analysis: Data needed for analysis

*Combined Figures used in publication of bacterial and fungal communities*
* Fungi_Bacteria_37_analysis: R code for community analysis including both fungi and bacteria creating Figure 3, 4, 5C, 6 and supplementary Figure S11, S13, S15, S17
* setup_F_B_paper1.R: Set up R environemnt
* data_Fungi_Bacteria_37_analysis: Data needed for R file

*Additional files*
* R_sessions_info: Sessions infos (including package versions loaded) for the different R scripts

*Additional comments*
* Some figures were manually redesigned in powerpoint and/or inkscape, for example to add icons, readjust legends or readjust axis, labels etc.
* If not spefically mentioned above, Table were manually assembled in Excle based on the different output files of analyses.








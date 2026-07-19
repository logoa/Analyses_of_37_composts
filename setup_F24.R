# Set up working environment for Fungi BLW2
# Author: Anja Logo
# Started: 09.09.24
# Last changes: 09.09.24

# Packages
library(conflicted) # To check for conflicts between packaged
library(data.table)
library(RColorBrewer) # additional color 
library(multcompView) # letters for direct comparisons
library(corrplot) # Correlations plots
library(rstatix) # Wilcoxon text
library(Hmisc)
library(ggpubr)
library(dunn.test) # non-paramteric direct comparison
library(ggrepel) # plot lables without overlap
library(adiv) # to claculate alpha diversity?
library(vegan) # Multivariate analysis
library(ALDEx2) # Differential abundance analysis
library(Maaslin2) # Differential abundance analysis
library(gtools) # data manipulation
library(parallel) #Parallelization of processes
#library(car)
library(plyr)
library(indicspecies) # For indicator species analysis
library(tidyverse) 
library(ggnewscale) # Allows more than one coloring in R! Very useful
library(VennDiagram)
library(ggVennDiagram)
library(gridExtra) # Grid.arrange needed to draw the venn diagramm
library(cowplot) # needed for arrange of alpha-diversity plots
library(seqinr) # To read fasta files into R
library(msa) # Sequence allignment
library(ape) # Conversion of DNA sequence format
library(ggtree)
library(ComplexHeatmap)
library(plotly)
library(scales)


# Conflicts of functions

conflicted::conflicts_prefer(dplyr::mutate)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::arrange)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::desc)
conflicted::conflicts_prefer(dplyr::summarise)
conflicted::conflicts_prefer(dplyr::summarize)
conflicted::conflicts_prefer(ggpubr::get_legend)
conflicted::conflicts_prefer(stats::cor)
conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(stats::sd)
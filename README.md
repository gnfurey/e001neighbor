# e001neighbor
Reproducible analysis for e001 Field C metacommunity manuscript. Code pulls data from the LTER pasta database.  
1_makefile_e001Cneigh.R
This file pulls the raw data and calculates species richness for each experimental plot, neighborhood species richness for each plot, and neighborhood abundance for each species in each plot. This file contains custom functions to calculate the mean neighborhood value for species richness in each year for each plot and for each species in each year for each plot. It also includes several checking functions to ensure the custom function did not make any calculation errors. 
sr = focal plot species richness
mean.Sr.4 = mean neighborhood species richness in the four adjacent plots
x.neighbor = where x is a species. Each species' neighborhood biomass. Each species is named with the first four letters of genus and species. 
2_SupplementalTable1_e001Cneigh.R
This file runs a linear mixed effects model on species richness using package nlme. It makes Supplemental Table S1.
sr = focal plot species richness
ln.NAtm.NAdd = log of added nitrogen per plot + 1
means.Sr.4_mean = mean neighborhood species richness for each plot 1982-2004 
3_Figure1_e001Cneigh.R
This file runs a generalized least squares model on species richness averaged across 1982-2004. spdep function sp.correlogram is used to calculate Moran's I statistic on the residuals from the regression. The output from sp.correlogram is pulled to show standard errors and extract p-values. It generates Figure 1.
4_Figure2_map_e001Cneigh.R
This file takes the species richness of each plot and its experimental treatment and makes a 2d plot using ggplot2. It generates Figure 2.
5_Figure3_e001Cneigh.R
This file tests a linear mixed effect model on species richness in the six control plots. It generates Figure 3 and Supplemental Table S2. 
6_Figure4_e001Cneigh.R
This file runs a separate linear regression for each control plot through time. It then compares those slopes with each plot's average neighborhood species richness using major axis regression in package sma. It generates Figure 4 and Supplemental Figure S1.
7_Figure5_e001Cneigh.R
This file runs a separate linear mixed effects model for three species using their abundance as the response variable. The dependency of each species' abundance is tested on its neighborhood abundance, the experimental nutrient addition variable (Factor of 9 levels), and year as a continuous variable with a three-way interaction. It then generates an interaction plot for each species at three experimental treatments. It generates Figure 5 and Supplemental Figure 2. 
AppendixS2_gifs_e001Cneigh.R
This file takes the abundance data for each species and then uses package gganimate to generate a video for each species. Videos are rendered into .mov format using ffmpeg. Each figure was created with R version 4.0.4 using package ggplot2 (3.3.2) (Wickham 2016, ggplot2 authors 2020), animated using R package gganimate (1.0.7)
helperfunctions_e001Cneigh.R
This file contains a few helper functions for matching species' metadata, several long figure captions, a function to calculate the standard error, and the negation of the function %in%. 

  
[![DOI](https://zenodo.org/badge/274988837.svg)](https://zenodo.org/badge/latestdoi/274988837)

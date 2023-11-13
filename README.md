# GutMicrobiome


# Gut microbiota patterns associated with duration of diarrhea in children under five years of age in Ethiopia

# Hardware requirements 
R package requires only a standard computer with enough RAM to support the in-memory operations.
# Software requirement 
 R software and all libraries supported by macOS Version 14.1
We used free software R version 4.1.2 and libraries such as Vegan v2.5-7, Phyloseq v1.38.0, DESeq2 v1.34.0, and mixOmics.  All R libraries used for the analysis of the microbiome data are included in the respective R script files. 


These R-scripts are for the research project entitled “Gut microbiota patterns associated with duration of diarrhea in children under five years of age in Ethiopia” 
Clinical data(metadata) used in this study cannot be made freely available to protect the privacy of the participants, in accordance with the Danish Data Protection Act and European Regulation 2016/679 of the European Parliament and of the Council (GDPR). However, data can be made available upon reasonable request after an agreement is reached with the University of Copenhagen and other participating institutions via the corresponding authors Getnet Tesfaw (gettesfaw2@gmail.com).

# There are four scripts associated with this project

# Diarrhea_Alpha_Beta_diversity_microbiome_analysis.R
This script contains the relative abundance analysis, alpha diversity analysis, beta diversity analysis, distance-based redundancy (db-RDA) analysis, and logistic regression analysis. 
  
# Phyloseq_CSS_import_16S.R
This R file is used to import the Phyloseq object that is used for the “beta diversity analysis and distance-based redundancy (db-RDA) analysis”.
  
# Phyloseq_import_16S.R
   
This R file is used to import the Phyloseq object that is used for the “relative abundance analysis, alpha diversity analysis”. 
  
# Diarrhea_Differenatial_abundance_Microbiome_analysis.Rmd

This R  markdown file contains the differential microbiome analysis using the DESeq2 and sparse partial least square discriminatory analysis (sPLS-DA).

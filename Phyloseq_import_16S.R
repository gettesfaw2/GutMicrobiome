library(phyloseq)
library(ggplot2)      # graphics
library(readxl)       # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(stringr)
library(vegan)
library(ggpubr)
library(metagenomeSeq)
library(tidyr)
library(RColorBrewer)
library(reshape2)


#Set working directory to r-script directory location - will not be run when script is run through source command

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################################Construct Phyloseq objects##########################################



##Load raw data cryptopoc########## 




otu_mat <- read.delim('zOTU_table_sintax.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1, header = TRUE)
dim(otu_mat)

map_mat <- read.delim('metadata_cryptopoc_final.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1, header = TRUE)

dim(map_mat)



#Extract taxonomy from otu-table
tax_mat <- dplyr::select(otu_mat,taxonomy)

#Remove "k__", "p__" etc from taxonomy column

tax_mat[,1] <- str_remove_all(tax_mat[,1], str_c(c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), collapse="|"))


#Split taxonomy column by ";"

split <- colsplit(tax_mat[,1], ";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

tax_mat <- cbind(tax_mat,split)

#Remove old taxonomy column
tax_mat <- tax_mat[,-c(1)]

tax_mat <- as.matrix.data.frame(tax_mat)

otu_mat[] <- mutate_all(otu_mat, function(x) as.numeric(as.character(x)))
otu_mat <- as.matrix(otu_mat)

##Construct individual tables is phyloseq format

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(map_mat)


##Combine in phyloseq object

PSB <- phyloseq(OTU, TAX, samples)

#Remove left over objects

rm(tax_mat, map_mat, otu_mat, samples, split, OTU, TAX)



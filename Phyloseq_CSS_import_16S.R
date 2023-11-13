###CSS Nomralization###



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(metagenomeSeq)
library(dplyr)
library(stringr)
library(phyloseq)
library(reshape2)

#######################Cecum#######################

##Load raw data


otu_mat <- read.delim('zOTU_table_sintax.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1, header = TRUE)
dim(otu_mat)

map_mat <- read.delim('metadata_cryptopoc_final.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1, header = TRUE)

dim(map_mat)


#Extract taxonomy from otu-table
tax_mat <- dplyr::select(otu_mat,taxonomy)

#######Do CSS norm

rownames <- rownames(otu_mat)

otu_mat <- within(otu_mat,rm("taxonomy"))

otu_mat[] <- mutate_all(otu_mat, function(x) as.numeric(as.character(x)))

data.metagenomeSeq = newMRexperiment(otu_mat, #phenoData=phenotypeData, 
                                     featureData=NULL, libSize=NULL, normFactors=NULL) #using filtered data 
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
otu_mat = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)

row.names(otu_mat) <- rownames

otu_mat <- as.matrix(otu_mat)

#########Phyloseq import

##Format taxonomy table

#Remove "k__", "p__" etc from taxonomy column

tax_mat[,1] <- str_remove_all(tax_mat[,1], str_c(c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), collapse="|"))

#Split taxonomy column by ";"

split <- colsplit(tax_mat[,1], ";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

tax_mat <- cbind(tax_mat,split)

tax_mat <- tax_mat[,-c(1)]

#Set as matrix
tax_mat <- as.matrix.data.frame(tax_mat)

##Construct individual tables is phyloseq format

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(map_mat)

##Combine in phyloseq object

PSB.CSS <- phyloseq(OTU, TAX, samples)


#Remove leftover 

rm(tax_mat, map_mat, otu_mat, samples, split, OTU, TAX, data.cumnorm, data.metagenomeSeq, p, rownames)


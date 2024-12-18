library(MSstatsConvert)
library(MSstats)
library(artMS)
setwd("C:/Users/livei/Documents/Proteomics_Data_Analysis/ANAC017_PSB-689/totalMQ")

# Read in MaxQuant files
proteinGroups <- read.delim("proteinGroups.txt", sep="\t", header=TRUE)
infile <- read.delim("evidence.txt", sep="\t", header=TRUE)

# Read in annotation including condition and biological replicates per run.
## Users should make these annotation files. It is not the output from MaxQuant.
annot <- read.delim("Annot.txt",sep='\t', header=TRUE)
keys <- read.delim("Keys.txt",sep='\t', header=TRUE)
## MSstats and artMS need different types of annotation files, annot and keys respectively.

# QC using artMS package
artmsQualityControlEvidenceBasic(
  evidence_file = infile,
  keys_file = keys,
  prot_exp = "APMS") ## Need to describe which type of proteomics experiment it was: APMS for TurboID

artmsQualityControlEvidenceExtended(
  evidence_file = infile,
  keys_file = keys) ## Look in the Packages help to see which QC plots you want

# Change MQ output format to make it compatible with MSstats
raw <- MaxQtoMSstatsFormat(evidence=infile, 
                           annotation=annot, 
                           proteinGroups=proteinGroups)

# Data processing
QuantData <- dataProcess(raw,
                         normalization = 'equalizeMedians',
                         summaryMethod = 'TMP',
                         censoredInt = "NA",
                         MBimpute = TRUE,
                         maxQuantileforCensored=0.999)

# Create a contrast matrix to define the comparisons you want to make.
## Here I make 3 comparisons: NAC17_AA vs NAC17_mock, NAC17_mock vs GFP_mock and NAC17_AA vs GFP_AA.
NAC17comp <- matrix(c(1,-1,0,0),nrow=1)
NAC17_AAcomp <- matrix(c(1,0,-1,0),nrow=1)
NAC17_mockcomp <- matrix(c(0,1,0,-1),nrow=1)
comparison <- rbind(NAC17comp,NAC17_AAcomp,NAC17_mockcomp)


colnames(comparison) <- c("NAC17_AA","NAC17_mock","GFP_AA","GFP_mock") # Colnames need to be the same as the conditions in your annotation file
row.names(comparison) <- c("NAC17comp","NAC17_AAcomp","NAC17_mockcomp")
comparison

# Tests for differentialy abundant proteins:
testResultMultiComparisons <- groupComparison(contrast.matrix=comparison, data=QuantData)
write.table(testResultMultiComparisons$ComparisonResult, file= "TurboID_results.txt", sep='\t', quote=F)

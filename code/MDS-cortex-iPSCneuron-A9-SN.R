# MDS plots of Lipton's A9 marker genes comparing Vivek's cortex samples,
# Yuan's iPSC neurons, and Lipton's A9 cultures and human substantia nigra
# samples

# Inputs:
#   CQN GC and read length normalized and read depth filtered HTseq expression
#   values Metadata with RIN values for each sample

# Outputs:
#   MDS plot all genes
#   MDS plot of marker genes
#   MDS plot of CACNA1D

sessionInfo()

library(reshape2)
library(ggplot2)
library(biomaRt)

# Input file paths (choose 1)
# Load CQN normalized expression values for hESC A9 and human substantia nigra
# samples
load("../processed_data/HTSeqUnion_Gene_CQN_OutlierRemoved_cortex_iPSCneuron_humanSN_RDF5.rda")
load("../processed_data/HTSeqUnion_Gene_CQN_OutlierRemoved_A9_SN_cortex_iPSCneuron_RDF5_regRINtotalReads.rda")
exprDat <- as.data.frame(cqnDat)

# Output file paths and variables
outInfo <- " cortex iPSCneuron humanSN A9 readDF5 CQN regRINtotalReads"
outpathAllGenes <- paste(
  "../analysis/MDS all genes"
  , outInfo, ".pdf", sep="")
outpathMarks <- paste(
  "../analysis/MDS marker genes"
  , outInfo, ".pdf", sep="")
outpathCAC <- paste(
  "../analysis/MDS CACNA1D"
  , outInfo, ".pdf", sep="")

# Other variables
graphSubTitle <- paste("\nread depth filter: 5"
                     , "\nregressed out RIN and total Reads"
                     , "\nMDS-cortex-iPSCneuron-A9-SN.R", sep = "")
type <- as.factor(c(rep("cortex", 9), rep("2", 3)
                    , rep("human substantia nigra", 5), rep("7",3)
                    , rep("human substantia nigra", 4)
                    , rep("iPSC neuron", 8)))
# Data frame of Ensembl IDs, gene symbols, and status as marker or anti-marker
markersDF <- data.frame(
  ensembl = c(
    "ENSG00000180176",
    "ENSG00000157542",
    "ENSG00000165646",
    "ENSG00000157388",
    "ENSG00000162761",
    "ENSG00000125798",
    "ENSG00000153234",
    "ENSG00000165092",
    
    "ENSG00000104327",
    "ENSG00000172137",
    "ENSG00000165588",
    "ENSG00000165655",
    "ENSG00000148680"
  )
  , gene = c(
    #markers
    "TH", "KCNJ7","VMAT2","CACNA1D","LMX1A","FOXA2","NURR1","ALDH1",
    #anti markers
    "CALB1", "CALB2", "OTX2", "NOLZ1", "HTR7"
    #"NOLZ1"(also known as ZNF503)
  )
  , type = c(rep("marker",8),rep("anti-marker",5)))

# Functions
# Function to output data frame of MDS PCs values and PCs
calcMDS <- function (exprDF) {
  # dist calculates the distances between the rows of a data matrix
  # Transpose data frame so samples are rows and genes are columns
  mds = cmdscale(dist(t(exprDF)), eig = T)
  pc1 = mds$eig[1]^2 / sum(mds$eig^2)
  pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  mdsAndTreatmentLDF <- data.frame(mds$points
                                   , pc1 = pc1, pc2 = pc2)
  mdsAndTreatmentLDF
}
################################################################################

# MDS All Genes

# Calculate MDS
mdsDF <- calcMDS(exprDat)
# Add column with sample type info
mdsDF$type <- type

ggplot(mdsDF, aes(x = X1, y = X2)) +
  geom_point(aes(color = factor(type)), size = 4) +
  geom_text(aes(label = row.names(mdsDF)), vjust = -1) +
  scale_color_discrete(name = "Sample Type"
                       , labels = c("(2) High MEF2C", "(7) Low MEF2C"
                                    , "Human Cortex", "Human Substantia Nigra"
                                    , "iPSC neuron")) +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(
    "MDS plot: All genes expression - Vivek Cortex, Lipton iPSC A9 and human SN"
    , graphSubTitle, sep = "")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = outpathAllGenes, height = 9)

# MDS Marker Genes

# Merge with Lipton expression values
markExprDF <- merge(markersDF, exprDat, by.x = "ensembl", by.y = "row.names")

# Calculate MDS
mdsDF <- calcMDS(markExprDF[ ,4:ncol(markExprDF)])
# Add column with sample type info
mdsDF$type <- type

ggplot(mdsDF, aes(x = X1, y = X2)) +
  geom_point(aes(color = factor(type)), size = 4) +
  geom_text(aes(label = row.names(mdsDF)), vjust = -1) +
  scale_color_discrete(name = "Sample Type"
                       , labels = c("(2) High MEF2C", "(7) Low MEF2C"
                       , "Human Cortex", "Human Substantia Nigra"
                       , "iPSC neuron")) +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(
    "MDS plot: A9 markers expression - Vivek Cortex, Lipton iPSC A9 and human SN"
    , graphSubTitle, sep = "")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = outpathMarks, height = 9)

# MDS CACNA1D

# Subset to CACNA1D expression values
cACNA1Ddat <- markExprDF[markExprDF$gene == "CACNA1D", ]

# Calculate MDS
mdsDF <- calcMDS(cACNA1Ddat[ ,4:ncol(cACNA1Ddat)])
# Add column with sample type info
mdsDF$type <- type

ggplot(mdsDF, aes(x = X1, y = X2)) +
  geom_point(aes(color = factor(type)), size = 4) +
  geom_text(aes(label = row.names(mdsDF)), vjust = -1) +
  scale_color_discrete(name = "Sample Type"
                       , labels = c("(2) High MEF2C", "(7) Low MEF2C"
                                    , "Human Cortex", "Human Substantia Nigra"
                                    , "iPSC neuron")) +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(
    "MDS plot: CACNA1D expression - Vivek Cortex, Lipton iPSC A9 and human SN"
    , graphSubTitle, sep = "")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = outpathCAC, height = 9)

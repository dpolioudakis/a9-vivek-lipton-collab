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

# Parameters
# outpathInfo <- " cortex-iPSCneuron-humanSN-iPSCa9 readDF5 CQN-GC-geneLength-quantile regRIN"
outpathInfo <- " cortex-humanSN-iPSCa9 readDF5 CQN-GC-geneLength-quantile regRIN"
# graphSubTitle <- paste("\nVivek Cortex, Yuan iPSC Neurons, Lipton iPSC A9 and human SN"
graphSubTitle <- paste("\nVivek Cortex, Lipton iPSC A9 and human SN"                       
                       , "\nread depth filter: 5"
                       , "\nCQN GC, gene length, quantile"
                       , "\nregressed out RIN"
                       , "\nMDS-cortex-iPSCneuron-A9-SN.R", sep = "")
# For properly labeling samples after MDS calculation
type <- as.factor(c(rep("cortex", 9), rep("2", 3)
                    , rep("human substantia nigra", 5), rep("7",3)
                    , rep("human substantia nigra", 4)
                    , rep("iPSC neuron", 8)))
# type <- as.factor(c(rep("cortex", 9), rep("2", 3)
#                     , rep("human substantia nigra", 5), rep("7",3)
#                     , rep("human substantia nigra", 5)))


# Input file paths (choose 1)
# Load CQN normalized expression values for hESC A9 and human substantia nigra
# samples
load("../processed_data/HTSeqUnion_Gene_CQN_OutlierRemoved_cortex_iPSCneuron_humanSN_RDF5.rda")
load("../processed_data/HTSeqUnion_Gene_CQN_OutlierRemoved_A9_SN_cortex_iPSCneuron_RDF5_regRINtotalReads.rda")
load("../processed_data/HTSeqUnion_Gene_CQN_OutlierRemoved_A9_SN_cortex_iPSCneuron_RDF5_regRINalignedReads.rda")
load("../processed_data/HTSeqUnion_Gene_A9-SN-cortex_RDF5_CQN-geneLength-GC-quantile_OutlierRemoved_regRIN.rda")
load("../processed_data/HTSeqUnion_Gene_A9-SN-cortex-iPSCneuron_RDF5_CQN-geneLength-GC-quantile_OutlierRemoved_regRIN.rda")
exprDatDF <- as.data.frame(exprRegM)

# Output file paths and variables
outpathAllGenes <- paste(
  "../analysis/MDS all genes"
  , outpathInfo, ".pdf", sep="")
outDendroAllGenes <- paste(
  "../analysis/Dendro all genes"
  , outpathInfo, ".pdf", sep="")
outpathMarks <- paste(
  "../analysis/MDS marker genes"
  , outpathInfo, ".pdf", sep="")
outDendroMarks <- paste(
  "../analysis/Dendro marker genes"
  , outpathInfo, ".pdf", sep="")
outpathCAC <- paste(
  "../analysis/MDS CACNA1D"
  , outpathInfo, ".pdf", sep="")
outDendroCAC <- paste(
  "../analysis/Dendro CACNA1D genes"
  , outpathInfo, ".pdf", sep="")
outpathAllen <- paste(
  "../analysis/MDS Allen"
  , outpathInfo, ".pdf", sep="")
outDendroAllen <- paste(
  "../analysis/Dendro Allen genes"
  , outpathInfo, ".pdf", sep="")

# Other variables
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
mdsDF <- calcMDS(exprDatDF)
# Add column with sample type info
mdsDF$type <- type
centroids <- aggregate(cbind(X1, X2)~type, mdsDF, mean)

ggplot(mdsDF, aes(x = X1, y = X2, color = factor(type))) +
  geom_point(size = 3) +
  geom_point(data = centroids, size = 6, shape = 3) +
  # geom_text(aes(label = row.names(mdsDF)), vjust = -1) +
  scale_color_discrete(name = "Sample Type"
                       , labels = c("(2) High MEF2C", "(7) Low MEF2C"
                                    , "Human Cortex", "Human Substantia Nigra"
                                    , "iPSC neuron")) +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(
    "MDS plot: All genes expression", graphSubTitle, sep = "")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(aspect.ratio = 4/5)
ggsave(file = outpathAllGenes, height = 9)

formattedDatDF <- t(exprDatDF)
rownames(formattedDatDF) <- c(rep("Cortex", 9), rep("High_MEF2C", 3)
                      , rep("Substantia_Nigra", 5), rep("Low_MEF2C", 3)
                     , rep("Substantia_Nigra", 4), rep("iPSC_Neuron", 8))
colorCodes <- c(Cortex = "darkgreen", High_MEF2C = "red", Low_MEF2C = "gold"
                , Substantia_Nigra = "blue", iPSC_Neuron = "purple")
distObj <- dist(formattedDatDF)
# Function to set label color
labelCol <- function(x) {
  if (is.leaf(x)) {
    # fetch label
    label <- attr(x, "label")
    attr(x, "nodePar") <- list(lab.col = colorCodes[label], pch = NA)
  }
  return(x)
}
hC <- hclust(distObj, method = "average")
hC <- dendrapply(as.dendrogram(hC), labelCol)
pdf(file = outDendroAllGenes)
par(mar = c(8,4,9,2))
plot(hC, ylab = "Height", main = paste("Hierarchical Clustering: All Genes"
                                       , "\nEuclidean Distance, Average Linkage"
                                       , graphSubTitle, sep = ""))
dev.off()
dM <- as.matrix(distObj)
sapply(c("High_MEF2C", "Low_MEF2C", "iPSC_Neuron", "Cortex")
       , function(x) mean(dM[rownames(dM) == "Substantia_Nigra"
                             , colnames(dM) == x]))

# MDS Marker Genes

# Merge with Lipton expression values
markExprDF <- merge(markersDF, exprDatDF, by.x = "ensembl", by.y = "row.names")

# Calculate MDS
mdsDF <- calcMDS(markExprDF[ ,4:ncol(markExprDF)])
# Add column with sample type info
mdsDF$type <- type
centroids <- aggregate(cbind(X1, X2)~type, mdsDF, mean)

ggplot(mdsDF, aes(x = X1, y = X2, color = factor(type))) +
  geom_point(size = 3) +
  geom_point(data = centroids, size = 6, shape = 3) +
  # geom_text(aes(label = row.names(mdsDF)), vjust = -1) +
  scale_color_discrete(name = "Sample Type"
                       , labels = c("(2) High MEF2C", "(7) Low MEF2C"
                       , "Human Cortex", "Human Substantia Nigra"
                       , "iPSC neuron")) +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(
    "MDS plot: A9 markers expression", graphSubTitle, sep = "")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(aspect.ratio = 4/5)
ggsave(file = outpathMarks, height = 9)

formattedDatDF <- t(markExprDF[ ,4:ncol(markExprDF)])
rownames(formattedDatDF) <- c(rep("Cortex", 9), rep("High_MEF2C", 3)
                              , rep("Substantia_Nigra", 5), rep("Low_MEF2C", 3)
                              , rep("Substantia_Nigra", 4), rep("iPSC_Neuron", 8))
distObj <- dist(formattedDatDF)
hC <- hclust(distObj, method = "average")
hC <- dendrapply(as.dendrogram(hC), labelCol)
pdf(file = outDendroMarks)
par(mar = c(8,4,9,2))
plot(hC, ylab = "Height", main = paste("Hierarchical Clustering: Marker Genes"
                                       , "\nEuclidean Distance, Average Linkage"
                                       , graphSubTitle, sep = ""))
dev.off()
dM <- as.matrix(distObj)
sapply(c("High_MEF2C", "Low_MEF2C", "iPSC_Neuron", "Cortex")
       , function(x) mean(dM[rownames(dM) == "Substantia_Nigra"
                             , colnames(dM) == x]))

# MDS CACNA1D

# Subset to CACNA1D expression values
cACNA1Ddat <- markExprDF[markExprDF$gene == "CACNA1D", ]

# Calculate MDS
mdsDF <- calcMDS(cACNA1Ddat[ ,4:ncol(cACNA1Ddat)])
# Add column with sample type info
mdsDF$type <- type
centroids <- aggregate(cbind(X1, X2)~type, mdsDF, mean)

ggplot(mdsDF, aes(x = X1, y = X2, color = factor(type))) +
  geom_point(size = 3) +
  geom_point(data = centroids, size = 6, shape = 3) +
  # geom_text(aes(label = row.names(mdsDF)), vjust = -1) +
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
  theme(axis.text = element_text(color = "black")) +
  theme(aspect.ratio = 4/5)
ggsave(file = outpathCAC, height = 9)

formattedDatDF <- t(data.frame(cACNA1Ddat[ ,4:ncol(cACNA1Ddat)]))
rownames(formattedDatDF) <- c(rep("Cortex", 9), rep("High_MEF2C", 3)
                              , rep("Substantia_Nigra", 5), rep("Low_MEF2C", 3)
                              , rep("Substantia_Nigra", 4), rep("iPSC_Neuron", 8))
distObj <- dist(formattedDatDF)
hC <- hclust(distObj, method = "average")
hC <- dendrapply(as.dendrogram(hC), labelCol)
pdf(file = outDendroCAC)
par(mar = c(8,4,9,2))
plot(hC, ylab = "Height", main = paste("Hierarchical Clustering: CACNA1D"
                                       , "\nEuclidean Distance, Average Linkage"
                                       , graphSubTitle, sep = ""))
dev.off()
dM <- as.matrix(distObj)
sapply(c("High_MEF2C", "Low_MEF2C", "iPSC_Neuron", "Cortex")
       , function(x) mean(dM[rownames(dM) == "Substantia_Nigra"
                             , colnames(dM) == x]))

################################################################################

# MADS A9 marker genes from Allen

# Other variables
# Data frame of Ensembl IDs, gene symbols, and status as marker or anti-marker
markersDF <- data.frame(
  ensembl = c(
    "ENSG00000180176",
    "ENSG00000157542",
    "ENSG00000165646",
    "ENSG00000157388",
    "ENSG00000165092",
    
    "ENSG00000104327"
  )
  , gene = c(
    #markers
    "TH", "KCNJ7","VMAT2","CACNA1D","ALDH1",
    #anti markers
    "CALB1"
    #"NOLZ1"(also known as ZNF503)
  )
  , type = c(rep("marker",5),rep("anti-marker",1)))

# Merge with Lipton expression values
markExprDF <- merge(markersDF, exprDatDF, by.x = "ensembl", by.y = "row.names")

# Calculate MDS
mdsDF <- calcMDS(markExprDF[ ,4:ncol(markExprDF)])
# Add column with sample type info
mdsDF$type <- type

ggplot(mdsDF, aes(x = X1, y = X2)) +
  geom_point(aes(color = factor(type)), size = 4) +
  # geom_text(aes(label = row.names(mdsDF)), vjust = -1) +
  scale_color_discrete(name = "Sample Type"
                       , labels = c("(2) High MEF2C", "(7) Low MEF2C"
                                    , "Human Cortex", "Human Substantia Nigra"
                                    , "iPSC neuron")) +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(
    "MDS plot: Allen A9 markers expression", graphSubTitle, sep = "")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(aspect.ratio = 4/5)
ggsave(file = outpathAllen, height = 9)

formattedDatDF <- t(markExprDF[ ,4:ncol(markExprDF)])
rownames(formattedDatDF) <- c(rep("Cortex", 9), rep("High_MEF2C", 3)
                              , rep("Substantia_Nigra", 5), rep("Low_MEF2C", 3)
                              , rep("Substantia_Nigra", 4), rep("iPSC_Neuron", 8))
distObj <- dist(formattedDatDF)
hC <- hclust(distObj, method = "average")
hC <- dendrapply(as.dendrogram(hC), labelCol)
pdf(file = outDendroAllen)
par(mar = c(8,4,9,2))
plot(hC, ylab = "Height", main = paste("Hierarchical Clustering: Allen Marker Genes"
                                       , "\nEuclidean Distance, Average Linkage"
                                       , graphSubTitle, sep = ""))
dev.off()
dM <- as.matrix(distObj)
sapply(c("High_MEF2C", "Low_MEF2C", "iPSC_Neuron", "Cortex")
       , function(x) mean(dM[rownames(dM) == "Substantia_Nigra"
                             , colnames(dM) == x]))

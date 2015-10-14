# Make a synthetic eigengene of A9 marker genes and make boxplot of ME
# expression in Lipton's A9 cultures and human substantia nigra

print("#######################################################################")
print("Starting synthetic-A9-marker-eigengene.R script...")
sessionInfo()

library(WGCNA)
library(reshape2)
library(ggplot2)
library(biomaRt)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

# Vivek normalized RNAseq FPKMs
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_regSN.rda")
load("../processed_data/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_regRIN260280.rda")
exprDataA9sNdF <- as.data.frame(exprDataRegM)
load("HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_CQNtogether_regRIN.rda")
# Data frame of A9 samples and human SN samples, each column is a sample
exprDataA9sNdF <- as.data.frame(normExpr.reg)
load("../processed_data/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_regAgeSexPMiRIn260280.rda")
exprDataA9sNdF <- as.data.frame(exprDataRegM)
load("../processed_data/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_FDrRIn260280.rda")
exprDataA9sNdF <- as.data.frame(exprFDRfiltM)
load("../processed_data/HTSeqUnion_Gene_A9-SN-cortex-iPSCneuron_RDF5_CQN-geneLength-GC-quantile_OutlierRemoved_regRIN.rda")
exprDataA9sNdF <- as.data.frame(exprRegM)

readDepthFilt <- "5"

# Vivek normalized RNAseq FPKMs from hESC derived A9 neuronal cultures
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN_5.rda")
# Vivek normalized RNAseq FPKMs from hESC derived A9 neuronal cultures
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_5.rda")

print("#######################################################################")

# Boxplot of synthetic eigengene expressions in Lipton's hESC A9 cultures
# Damon's A9 markers

markerGenes <- c("ALDH1A1", "TH", "SLC18A2", "KCNJ6")
exprA9sNdF <- exprDataA9sNdF

AddEnsembl <- function (geneList) {
  moduleGenes <- data.frame(geneList)
  # bioMart manual:
  # http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  moduleEnsemblDF <- getBM(  attributes = c("ensembl_gene_id", "hgnc_symbol")
                             , filters = "hgnc_symbol"
                             , values = moduleGenes
                             , mart = ensembl)
  moduleEnsemblDF
}

# DF of Ensembl ID and Gene symbol
genesEnsemblDF <- AddEnsembl(markerGenes)

print("Marker genes in Lipton data:")
table(row.names(exprA9sNdF) %in% genesEnsemblDF[ ,1])
exprA9sNdF$marker <- "grey"
exprA9sNdF$marker[row.names(exprA9sNdF) %in% genesEnsemblDF[ ,1]] <- "marker"
exprA9sNdF <- exprA9sNdF[exprA9sNdF$marker == "marker", ]


markerMEa9sNdF <- moduleEigengenes(t(exprA9sNdF[ ,-(ncol(exprA9sNdF))])
                                  , exprA9sNdF$marker)$eigengenes
markerMEa9sNdF$biorep <- c( rep(2, 3)
                          , rep(7, 3)
                          , rep("human", 10))
markerMEa9sNdF <- melt(markerMEa9sNdF, id.vars = "biorep")
colnames(markerMEa9sNdF) <- c("biorep", "module", "MEexpression")

# Boxplot of marker modules

ggplot(data = markerMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Sample Type",
                      labels= c( "(2) High MEF2C"
                               , "(7) Low MEF2C"
                               , "Human Subsantia Nigra")) +
  ylab("ME Expression (arbitrary value)") +
  xlab("Sample Type") +
  labs(title = paste(
      "synthetic-A9-marker-eigengene.R"
    , "\nSynthetic ME A9 Damon's marker expression in"
    , "\nLipton A9 and human substantia nigra"
    , "\nMarker genes:", c(list(markerGenes))
    , "\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste(
  "../analysis/Synthetic Damon's A9 marker ME expression readDF", readDepthFilt
  , " regRINratio260280.pdf", sep=""), width = 14, height = 6)

# Code to check module eigengene expression without splitting by sample type
ggplot(data = markerMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot()
print("################################")

# Boxplot of synthetic eigengene expressions in Lipton's hESC A9 cultures
# Liptons's A9 markers

markerGenes <- c("ALDH1A1", "TH", "SLC18A2", "CACNA1D"
                 , "KCNJ6", "LMX1A", "FOXA2", "NR4A2", "ALDH1A1")
exprA9sNdF <- exprDataA9sNdF

# DF of Ensembl ID and Gene symbol
genesEnsemblDF <- AddEnsembl(markerGenes)

print("Marker genes in Lipton data:")
table(row.names(exprA9sNdF) %in% genesEnsemblDF[ ,1])
exprA9sNdF$marker <- "grey"
exprA9sNdF$marker[row.names(exprA9sNdF) %in% genesEnsemblDF[ ,1]] <- "marker"
exprA9sNdF <- exprA9sNdF[exprA9sNdF$marker == "marker", ]


markerMEa9sNdF <- moduleEigengenes(t(exprA9sNdF[ ,-(ncol(exprA9sNdF))])
                                   , exprA9sNdF$marker)$eigengenes
markerMEa9sNdF$biorep <- c( rep(2, 3)
                            , rep(7, 3)
                            , rep("human", 10))
markerMEa9sNdF <- melt(markerMEa9sNdF, id.vars = "biorep")
colnames(markerMEa9sNdF) <- c("biorep", "module", "MEexpression")

# Boxplot of marker modules

ggplot(data = markerMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Sample Type",
                      labels= c(  "(2) High MEF2C"
                                , "(7) Low MEF2C"
                                , "Human Subsantia Nigra")) +
  ylab("ME Expression (arbitrary value)") +
  xlab("Sample Type") +
  labs(title = paste(
    "synthetic-A9-marker-eigengene.R"
    , "\nSynthetic ME Lipton's A9 marker expression in"
    , "\nLipton A9 and human substantia nigra"
    , "\nMarker genes:", c(list(markerGenes))
    , "\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste(
  "../analysis/Synthetic Lipton's A9 marker ME expression readDF", readDepthFilt
  , " regRINratio260280.pdf", sep=""), width = 14, height = 6)
print("################################")

# Boxplot of synthetic eigengene expressions in Lipton's hESC A9 cultures
# Liptons's A9 anti-markers

markerGenes <- c("CALB1", "CALB2")

exprA9sNdF <- exprDataA9sNdF

# DF of Ensembl ID and Gene symbol
genesEnsemblDF <- AddEnsembl(markerGenes)

print("Marker genes in Lipton data:")
table(row.names(exprA9sNdF) %in% genesEnsemblDF[ ,1])
exprA9sNdF$marker <- "grey"
exprA9sNdF$marker[row.names(exprA9sNdF) %in% genesEnsemblDF[ ,1]] <- "marker"
exprA9sNdF <- exprA9sNdF[exprA9sNdF$marker == "marker", ]


markerMEa9sNdF <- moduleEigengenes(t(exprA9sNdF[ ,-(ncol(exprA9sNdF))])
                                   , exprA9sNdF$marker)$eigengenes
markerMEa9sNdF$biorep <- c( rep(2, 3)
                            , rep(7, 3)
                            , rep("human", 10))
markerMEa9sNdF <- melt(markerMEa9sNdF, id.vars = "biorep")
colnames(markerMEa9sNdF) <- c("biorep", "module", "MEexpression")

# Boxplot of marker modules

ggplot(data = markerMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Sample Type",
                      labels= c(  "(2) High MEF2C"
                                  , "(7) Low MEF2C"
                                  , "Human Subsantia Nigra")) +
  ylab("ME Expression (arbitrary value)") +
  xlab("Sample Type") +
  labs(title = paste(
    "synthetic-A9-marker-eigengene.R"
    , "\nSynthetic ME Lipton's A9 anti-marker expression in"
    , "\nLipton A9 and human substantia nigra"
    , "\nMarker genes:", c(list(markerGenes))
    , "\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste(
  "../analysis/Synthetic Lipton's A9 anti-marker ME expression readDF", readDepthFilt
  , " regRINratio260280.pdf.pdf", sep=""), width = 14, height = 6)

print("#######################################################################")

# Developmental code checking median expression of marker genes

exprA9sNdF <- as.data.frame(datExpr.HTSC.A9SN)

head(exprA9sNdF)
exprMarkerDF <- data.frame(exprA9sNdF)
# Melt as matrix to keep row names
head(melt(exprMarkerDF), 20)
exprMarkerDF <- melt(as.matrix(exprMarkerDF))

exprMarkerDF$biorep <- c(  rep(2, 48240)
                             , rep(7, 48240)
                             , rep("human", 160800))

colnames(exprMarkerDF) <- c("ensembl_id", "sample", "expression", "biorep")

print("Marker genes in Lipton data:")
table(exprMarkerDF$ensembl_id %in% genesEnsemblDF[ ,1])
exprMarkerDF$marker <- "none"
exprMarkerDF$marker[exprMarkerDF$ensembl_id %in% genesEnsemblDF[ ,1]] <- "marker"
exprMarkerLDF <- split(exprMarkerDF, exprMarkerDF$marker)
exprMarkerLDF <- split(exprMarkerLDF$marker, exprMarkerLDF$marker$biorep)
sapply(exprMarkerLDF, function(x) summary(x$expression))

ggplot(data = exprMarkerDF, aes(x = module, y = MEexpression)) +
  geom_boxplot()
ggplot(data = exprMarkerDF, aes(x = marker, y = expression)) +
  geom_boxplot(aes(fill=as.factor(biorep)))
print("#######################################################################")

# Boxplot of synthetic eigengene expressions in Lipton's hESC A9 cultures and
# Human Substantia Nigra samples, Vivek's human cortex samples, and Yuan's
# iPSC Neurons

# A9 markers from Allen Brain Atlas
markerGenes <- c("ALDH1A1", "TH", "SLC18A2", "KCNJ6", "CACNA1D", "CALB1")
# Liptons's A9 markers
markerGenes <- c("ALDH1A1", "TH", "SLC18A2", "CACNA1D"
                 , "KCNJ6", "LMX1A", "FOXA2", "NR4A2")

exprA9sNdF <- exprDataA9sNdF

# DF of Ensembl ID and Gene symbol
genesEnsemblDF <- AddEnsembl(markerGenes)

print("Marker genes in Lipton data:")
table(row.names(exprA9sNdF) %in% genesEnsemblDF[ ,1])
exprA9sNdF$marker <- "grey"
exprA9sNdF$marker[row.names(exprA9sNdF) %in% genesEnsemblDF[ ,1]] <- "marker"
exprA9sNdF <- exprA9sNdF[exprA9sNdF$marker == "marker", ]


markerMEa9sNdF <- moduleEigengenes(t(exprA9sNdF[ ,-(ncol(exprA9sNdF))])
                                   , exprA9sNdF$marker)$eigengenes
markerMEa9sNdF$biorep <- as.factor(c(rep("cortex", 9), rep("2", 3)
                                     , rep("human substantia nigra", 5), rep("7",3)
                                     , rep("human substantia nigra", 4)
                                     , rep("iPSC neuron", 8)))
markerMEa9sNdF <- melt(markerMEa9sNdF, id.vars = "biorep")
colnames(markerMEa9sNdF) <- c("biorep", "module", "MEexpression")

# Boxplot of marker modules

ggplot(data = markerMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Sample Type",
                      labels= c(  "(2) High MEF2C"
                                  , "(7) Low MEF2C"
                                  , "Human Cortex"
                                  , "Human Subsantia Nigra"
                                  , "iPSC Neuron")) +
  ylab("ME Expression (arbitrary value)") +
  xlab("Sample Type") +
  labs(title = paste(
    "synthetic-A9-marker-eigengene.R"
    , "\nSynthetic ME Allen A9 marker expression in"
    , "\nLipton A9 and human substantia nigra, Vivek Cortex, Yuan iPSC neuron"
    , "\nMarker genes:", c(list(markerGenes))
    , "\nCQN normalized: GC, gene length, quantile"
    , "\nRIN regressed out"
    , "\nread depth filter: ", readDepthFilt
    , sep = "")) +
  theme_grey(base_size = 18) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text = element_text(color = "black")) +
  theme(plot.title = element_text(size = rel(0.6))) +
  theme(aspect.ratio = 4/4)
ggsave(file = paste(
  "../analysis/Synthetic Allen A9 marker ME expression readDF"
  , "cortex-humanSN-iPSCa9 readDF5 CQN-GC-geneLength-quantile regRIN.pdf"
  , sep=""))

# Boxplot for antimarkers
antiMarkGenes <- c("CALB1", "CALB2")
exprA9sNdF <- exprDataA9sNdF

# DF of Ensembl ID and Gene symbol
genesEnsemblDF <- AddEnsembl(antiMarkGenes)

print("Marker genes in Lipton data:")
table(row.names(exprA9sNdF) %in% genesEnsemblDF[ ,1])
exprA9sNdF$marker <- "grey"
exprA9sNdF$marker[row.names(exprA9sNdF) %in% genesEnsemblDF[ ,1]] <- "marker"
exprA9sNdF <- exprA9sNdF[exprA9sNdF$marker == "marker", ]


markerMEa9sNdF <- moduleEigengenes(t(exprA9sNdF[ ,-(ncol(exprA9sNdF))])
                                   , exprA9sNdF$marker)$eigengenes
markerMEa9sNdF$biorep <- as.factor(c(rep("cortex", 9), rep("2", 3)
                                     , rep("human substantia nigra", 5), rep("7",3)
                                     , rep("human substantia nigra", 4)
                                     , rep("iPSC neuron", 8)))
markerMEa9sNdF <- melt(markerMEa9sNdF, id.vars = "biorep")
colnames(markerMEa9sNdF) <- c("biorep", "module", "MEexpression")

# Boxplot of marker modules

ggplot(data = markerMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Sample Type",
                      labels= c(  "(2) High MEF2C"
                                  , "(7) Low MEF2C"
                                  , "Human Cortex"
                                  , "Human Subsantia Nigra"
                                  , "iPSC Neuron")) +
  ylab("ME Expression (arbitrary value)") +
  xlab("Sample Type") +
  labs(title = paste(
    "synthetic-A9-marker-eigengene.R"
    , "\nSynthetic ME Allen A9 marker expression in"
    , "\nLipton A9 and human substantia nigra, Vivek Cortex, Yuan iPSC neuron"
    , "\nMarker genes:", c(list(markerGenes))
    , "\nCQN normalized: GC, gene length, quantile"
    , "\nRIN regressed out"
    , "\nread depth filter: ", readDepthFilt
    , sep = "")) +
  theme_grey(base_size = 18) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text = element_text(color = "black")) +
  theme(plot.title = element_text(size = rel(0.6))) +
  theme(aspect.ratio = 4/4)
ggsave(file = paste(
  "../analysis/Synthetic Allen A9 anti-marker ME expression readDF"
  , "cortex-humanSN-iPSCa9 readDF5 CQN-GC-geneLength-quantile regRIN.pdf"
  , sep=""))

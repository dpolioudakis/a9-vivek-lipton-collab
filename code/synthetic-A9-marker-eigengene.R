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

readDepthFilt <- "5"

# Vivek normalized RNAseq FPKMs from hESC derived A9 neuronal cultures
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN_5.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN_0.rda")
# Vivek normalized RNAseq FPKMs from hESC derived A9 neuronal cultures
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_0.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_5.rda")

print("#######################################################################")

# Boxplot of synthetic eigengene expressions in Lipton's hESC A9 cultures
# Damon's A9 markers

markerGenes <- c("ALDH1A1", "TH", "SLC18A2", "KCNJ6")

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

# Data frame of A9 samples and human SN samples, each column is a sample
exprA9sNdF <- as.data.frame(datExpr.HTSC.A9SN)

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

markerMEinA9$module <- factor(markerMEinA9$module
                              , levels = as.character(unique(MEtoUse)))

# Boxplot of marker modules

ggplot(data = markerMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Sample Type",
                      labels= c( "(2) HighMEF2C"
                               , "(7) LowMEF2C"
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
  , ".pdf", sep=""), width = 14, height = 6)

# Code to check module eigengene expression without splitting by sample type
ggplot(data = markerMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot()
print("################################")

# Boxplot of synthetic eigengene expressions in Lipton's hESC A9 cultures
# Liptons's A9 markers

markerGenes <- c("ALDH1A1", "TH", "SLC18A2", "CACNA1D"
                 , "KCNJ6", "LMX1A", "FOXA2", "NR4A2", "ALDH1A1")

# DF of Ensembl ID and Gene symbol
genesEnsemblDF <- AddEnsembl(markerGenes)

# Data frame of A9 samples and human SN samples, each column is a sample
exprA9sNdF <- as.data.frame(datExpr.HTSC.A9SN)

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

markerMEinA9$module <- factor(markerMEinA9$module
                              , levels = as.character(unique(MEtoUse)))

# Boxplot of marker modules

ggplot(data = markerMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Sample Type",
                      labels= c(  "(2) HighMEF2C"
                                , "(7) LowMEF2C"
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
  , ".pdf", sep=""), width = 14, height = 6)

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


# Allen module eigengene expression in Lipton A9 and human substantia nigra

print("#######################################################################")
print("Starting allen-ME-expr-A9-SN.R script...")
sessionInfo()

library(WGCNA)
library(reshape2)
library(ggplot2)
library(biomaRt)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/allen_BW_modules.rda")
load("../processed_data/array_data_subset_avg_probes.rda")

load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN_5.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN_0.rda")
# Vivek normalized RNAseq FPKMs from hESC derived A9 neuronal cultures
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_0.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_5.rda")

# Normalized together RNA hESC A9 and human substantia nigra
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_regSN.rda")
datExpr.HTSC.A9 <- datExpr.HTSC.A9SN[ ,1:6]
datExpr.HTSC.SN <- datExpr.HTSC.A9SN[ ,7:16]
load("HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_CQNtogether_reg.rda")
datExpr.HTSC.A9 <- normExpr.reg[ ,1:6]
datExpr.HTSC.SN <- normExpr.reg[ ,7:16]

# variable for read depth filter to record in output graph titles
readDepthFilt <- "5"
minModSize <- "30"

# bwModulesLL is list of modules from different blockwiseModules parameters used
# 12 corresponds to softPower 9, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 12
modsToUse <- c("saddlebrown", "salmon", "red", "pink", "black", "green"
               , "grey60")
# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 11
modsToUse <- c("plum1", "grey60", "brown", "red", "cyan", "yellowgreen"
               , "sienna3", "royalblue")
# 12 corresponds to softPower 9, minModSize 100, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 18
modsToUse <- c("salmon", "brown", "blue", "purple", "greenyellow")
print("#######################################################################")

# Boxplot of Allen module eigengene expressions in Lipton's hESC A9 cultures 

genesColorsDF <- data.frame(colnames(exprData)
                            , bwModulesLL[[modNetworkToUse]]$colors)
colnames(genesColorsDF) <- c("gene", "module")

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

# List of data frames of Ensembl ID and Gene symbol
genesEnsemblDF <- AddEnsembl(colnames(exprData))
ensemblColorsDF <- merge(genesEnsemblDF, genesColorsDF
                         , by.x = "hgnc_symbol", by.y = "gene")

# Data frame of A9 samples and human SN samples, each column is a sample
exprA9sNdF <- merge(datExpr.HTSC.A9, datExpr.HTSC.SN
                    , by.x = "row.names", by.y = "row.names")
exprA9sNdF <- data.frame(exprA9sNdF[ ,-1], row.names = exprA9sNdF[ ,1])

modsA9sNdF <- merge(ensemblColorsDF, exprA9sNdF
                   , by.x = "ensembl_gene_id", by.y = "row.names" )

allenMEa9sNdF <- moduleEigengenes(t(modsA9sNdF[ ,4:18])
                                  , modsA9sNdF$module)$eigengenes
allenMEa9sNdF$biorep <- c(rep(2, 3), rep(7, 3), rep("human", 9))
allenMEa9sNdF <- melt(allenMEa9sNdF, id.vars = "biorep")
colnames(allenMEa9sNdF) <- c("biorep", "module", "MEexpression")

MEtoUse <- sapply(modsToUse, function(mod) paste("ME", mod, sep=""))
markerMEinA9 <- allenMEa9sNdF[allenMEa9sNdF$module %in% MEtoUse, ]

markerMEinA9$module <- factor(markerMEinA9$module
                              , levels = as.character(unique(MEtoUse)))

# Boxplot of marker modules
ggplot(data = markerMEinA9, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Biological\nReplicate",
                      labels= c(    "(2) High MEF2C"
                                  , "(7) Low MEF2C"
                                  , "Human Subsantia Nigra")) +
  ylab("ME Expression (arbitrary value)") +
  xlab("Module") +
  labs(title = paste(
      "allen-ME-expr-A9-SN.R"
    , "\nAllen marker ME expression in Lipton A9 and SN"
    , "\nA9 and SN samples CQN normalized together"
    , "\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste(
  "../analysis/Allen marker ME expr in A9 SN readDF", readDepthFilt
  , " ModSize", minModSize, " CQN together.pdf", sep=""))

# Boxplot all modules
ggplot(data = allenMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Biological\nReplicate",
                      labels= c(      "(2) High MEF2C"
                                    , "(7) Low MEF2C"
                                    , "Human Subsantia Nigra")) +
  ylab("ME Expression (arbitrary value)") +
  xlab("Module") +
  labs(title = paste(
      "allen-ME-expr-A9-SN.R"
    , "\nAllen ME expression in Lipton A9 and SN"
    , "\nA9 and SN samples CQN normalized together"
    , "\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste(
  "../analysis/Allen ME expr in A9 SN readDF", readDepthFilt
  , " ModSize", minModSize, " CQN together.pdf", sep=""))

# Boxplot all modules - not separated by treatment group
ggplot(data = allenMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot() +
  ylab("ME Expression (arbitrary value)") +
  xlab("Module") +
  labs(title = paste(
     "allen-ME-expr-A9-SN.R"
    , "\nAllen ME expression in Lipton A9 and SN"
    , "\nNot separated by treatment group"
    , "\nA9 and SN samples CQN normalized together"
    , "\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste(
  "../analysis/Allen marker ME expr in A9 SN combined Tx readDF", readDepthFilt
  , " ModSize", minModSize, " CQN together.pdf", sep=""))

# Developmental code to assign genes to random modules
modsA9sNdF$module <- sample(1:15, nrow(modsA9sNdF), replace=T)

allenMEa9sNdF <- moduleEigengenes(t(modsA9sNdF[ ,4:18])
                                  , modsA9sNdF$module)$eigengenes
allenMEa9sNdF$biorep <- c(rep(2, 3), rep(7, 3), rep("human", 9))
allenMEa9sNdF <- melt(allenMEa9sNdF, id.vars = "biorep")
colnames(allenMEa9sNdF) <- c("biorep", "module", "MEexpression")

MEtoUse <- sapply(modsToUse, function(mod) paste("ME", mod, sep=""))
markerMEinA9 <- allenMEa9sNdF[allenMEa9sNdF$module %in% MEtoUse, ]

markerMEinA9$module <- factor(markerMEinA9$module
                              , levels = as.character(unique(MEtoUse)))

# Boxplot all modules - genes randomly assigned to each module
ggplot(data = allenMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Biological\nReplicate",
                      labels= c(      "(2) High MEF2C"
                                      , "(7) Low MEF2C"
                                      , "Human Subsantia Nigra")) +
  ylab("ME Expression (arbitrary value)") +
  xlab("Module") +
  labs(title = paste(
    "allen-ME-expr-A9-SN.R"
    , "\nGenes randomly assigned to modules"
    , "\nA9 and SN samples CQN normalized together"
    , "\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste(
  "../analysis/Allen random ME expr in A9 SN readDF", readDepthFilt
  , " ModSize", minModSize, " CQN together.pdf", sep=""))
print("#######################################################################")


high <- moduleEigengenes(t(modsA9sNdF[ ,4:6])
                                   , modsA9sNdF$module)$eigengenes
allenMEa9sNdF$biorep <- c(rep(2, 3), rep(7, 3), rep("human", 9))
allenMEa9sNdF <- melt(allenMEa9sNdF, id.vars = "biorep")
colnames(allenMEa9sNdF) <- c("biorep", "module", "MEexpression")

MEtoUse <- sapply(modsToUse, function(mod) paste("ME", mod, sep=""))
markerMEinA9 <- allenMEa9sNdF[allenMEa9sNdF$module %in% MEtoUse, ]

markerMEinA9$module <- factor(markerMEinA9$module
                              , levels = as.character(unique(MEtoUse)))


low <- moduleEigengenes(t(modsA9sNdF[ ,7:9])
                                   , modsA9sNdF$module)$eigengenes
allenMEa9sNdF$biorep <- c(rep(2, 3), rep(7, 3), rep("human", 9))
allenMEa9sNdF <- melt(allenMEa9sNdF, id.vars = "biorep")
colnames(allenMEa9sNdF) <- c("biorep", "module", "MEexpression")

MEtoUse <- sapply(modsToUse, function(mod) paste("ME", mod, sep=""))
markerMEinA9 <- allenMEa9sNdF[allenMEa9sNdF$module %in% MEtoUse, ]

markerMEinA9$module <- factor(markerMEinA9$module
                              , levels = as.character(unique(MEtoUse)))


all <- moduleEigengenes(t(modsA9sNdF[ ,4:9])
                                   , modsA9sNdF$module)$eigengenes
allenMEa9sNdF$biorep <- c(rep(2, 3), rep(7, 3), rep("human", 9))
allenMEa9sNdF <- melt(allenMEa9sNdF, id.vars = "biorep")
colnames(allenMEa9sNdF) <- c("biorep", "module", "MEexpression")

MEtoUse <- sapply(modsToUse, function(mod) paste("ME", mod, sep=""))
markerMEinA9 <- allenMEa9sNdF[allenMEa9sNdF$module %in% MEtoUse, ]

markerMEinA9$module <- factor(markerMEinA9$module
                              , levels = as.character(unique(MEtoUse)))



# Data frame of A9 samples and human SN samples, each column is a sample
exprA9sNdF <- merge(datExpr.HTSC.A9, datExpr.HTSC.SN
                    , by.x = "row.names", by.y = "row.names")
exprA9sNdF <- data.frame(exprA9sNdF[ ,-1], row.names = exprA9sNdF[ ,1])

modsA9sNdF <- merge(ensemblColorsDF, exprA9sNdF
                    , by.x = "ensembl_gene_id", by.y = "row.names" )

mean(unlist(modsA9sNdF[, 4:9]))
mean(unlist(modsA9sNdF[, 4:6]))
mean(unlist(modsA9sNdF[, 7:9]))

summary(unlist(modsA9sNdF[, 4:9]))
summary(unlist(modsA9sNdF[, 4:6]))
summary(unlist(modsA9sNdF[, 7:9]))



modsA9sNdF$module <- sample(1:15, nrow(modsA9sNdF), replace=T)

allenMEa9sNdF <- moduleEigengenes(t(modsA9sNdF[ ,4:9])
                                   , modsA9sNdF$module)$eigengenes

high <- moduleEigengenes(t(modsA9sNdF[ ,4:6])
                                  , modsA9sNdF$module)$eigengenes
low <- moduleEigengenes(t(modsA9sNdF[ ,7:9])
                         , modsA9sNdF$module)$eigengenes
allenMEa9sNdF <- rbind(high, low)
allenMEa9sNdF$biorep <- c(rep(2, 3), rep(7, 3))

allenMEa9sNdF <- melt(allenMEa9sNdF, id.vars = "biorep")
colnames(allenMEa9sNdF) <- c("biorep", "module", "MEexpression")

MEtoUse <- sapply(modsToUse, function(mod) paste("ME", mod, sep=""))
markerMEinA9 <- allenMEa9sNdF[allenMEa9sNdF$module %in% MEtoUse, ]

markerMEinA9$module <- factor(markerMEinA9$module
                              , levels = as.character(unique(MEtoUse)))

# Boxplot of modules
ggplot(data = allenMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot(aes(fill=as.factor(biorep))) +
  scale_fill_discrete(name= "Biological\nReplicate",
                      labels= c(  "2_HighMEF2C"
                                  , "7_LowMEF2C"
                                  )) +
  ylab("ME Expression (arbitrary value)") +
  xlab("Module") +
  labs(title = paste(
    "allen-A9-marker-expr-in-hESC-A9.R\nAllen marker ME expression in"
    , "Lipton A9\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(color = "black"))



allenMEa9sNdF <- moduleEigengenes(t(modsA9sNdF[ ,10:18])
                                   , modsA9sNdF$module)$eigengenes
allenMEa9sNdF$biorep <- c(rep("human", 9))
allenMEa9sNdF <- melt(allenMEa9sNdF, id.vars = "biorep")
colnames(allenMEa9sNdF) <- c("biorep", "module", "MEexpression")

MEtoUse <- sapply(modsToUse, function(mod) paste("ME", mod, sep=""))
markerMEinA9 <- allenMEa9sNdF[allenMEa9sNdF$module %in% MEtoUse, ]

markerMEinA9$module <- factor(markerMEinA9$module
                              , levels = as.character(unique(MEtoUse)))

# Boxplot of marker modules
ggplot(data = allenMEa9sNdF, aes(x = module, y = MEexpression)) +
  geom_boxplot() +
  scale_fill_discrete(name= "Biological\nReplicate",
                      labels= c(  "Human Subsantia Nigr"
                      )) +
  ylab("ME Expression (arbitrary value)") +
  xlab("Module") +
  labs(title = paste(
    "allen-A9-marker-expr-in-hESC-A9.R\nAllen marker ME expression in"
    , "Lipton A9\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(color = "black"))




modsA9dF <- merge(ensemblColorsDF, datExpr.HTSC.A9
                    , by.x = "ensembl_gene_id", by.y = "row.names" )

mean(unlist(modsA9dF[modsA9dF$module == "plum1", 4:9]))
mean(unlist(modsA9dF[modsA9dF$module == "plum1", 4:6]))
mean(unlist(modsA9dF[modsA9dF$module == "plum1", 7:9]))

allenMEa9dF <- moduleEigengenes(t(modsA9dF[4:9])
                                   , modsA9dF$module)$eigengenes
allenMEa9sNdF$biorep <- c(rep(2, 3), rep(7, 3), rep("human", 9))
allenMEa9sNdF <- melt(allenMEa9sNdF, id.vars = "biorep")
colnames(allenMEa9sNdF) <- c("biorep", "module", "MEexpression")

MEtoUse <- sapply(modsToUse, function(mod) paste("ME", mod, sep=""))
markerMEinA9 <- allenMEa9sNdF[allenMEa9sNdF$module %in% MEtoUse, ]

markerMEinA9$module <- factor(markerMEinA9$module
                              , levels = as.character(unique(MEtoUse)))
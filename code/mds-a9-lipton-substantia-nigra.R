# MDS plots of Lipton hESC A9 neuronal cultures and human substantia nigra (SN)
# samples

# Inputs:
#   CQN normalized and read depth filtered expression (normalized FPKM) tables
#   of gene expression for hESC A9 and human SN samples
#     Cufflinks&HTseqCounts_LiptonSN.R - Vivek's script

# Outputs:
#   MDS plot of all genes for hESC A9 and human SN
#   MDS plot of all genes for hESC A9 and human SN subset by Allen modules


sessionInfo()

library(reshape2)
library(ggplot2)
library(biomaRt)

load("../processed_data/allen_BW_modules.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_5.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_0.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN_5.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN_0.rda")

# variable for read depth filter to record in output graph titles
readDepthFilt <- "5"

# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 11
modsToUse <- c("plum1", "grey60", "brown", "red", "cyan", "yellowgreen"
               , "sienna3", "royalblue")

# Data frame of A9 samples and human SN samples, each column is a sample
exprA9sNdF <- merge(datExpr.HTSC.A9, datExpr.HTSC.SN
                    , by.x = "row.names", by.y = "row.names")
exprA9sNdF <- data.frame(exprA9sNdF[ ,-1], row.names = exprA9sNdF[ ,1])
print("#######################################################################")

# MDS of all genes for hESC A9 and human SN

# Returns data frame with distances as columns and samples as rows
calcMDS <- function (exprDF) {
  # dist calculates the distances between the rows of a data matrix
  # Transpose data frame so samples are rows and genes are columns
  mds = cmdscale(dist(t(exprDF)), eig = T)
  pc1 = mds$eig[1]^2 / sum(mds$eig^2)
  pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  mdsAndTreatmentLDF <- list(distVals = data.frame(mds$points)
                             , pc1 = pc1, pc2 = pc2)
  mdsAndTreatmentLDF
}

mdsLDF <- calcMDS(exprA9sNdF)
mdsLDF$distVals$type <- as.factor(c(rep("2",3), rep("7",3), rep("human",9)))

ggplot(mdsLDF$distVals, aes(x = X1, y = X2)) +
  geom_point(aes(color = factor(type)), size = 4) +
  geom_text(aes(label = row.names(mdsLDF$distVals)), vjust = -1) +
  scale_color_discrete(name = "Sample Type"
                       , labels = c("2_HighMEF2C", "7_LowMEF2C"
                                   , "Human Substantia Nigra")) +
  xlab(paste("PC1 (", signif(100*mdsLDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsLDF$pc2, 3),"%)", sep = "")) +
  labs(title = paste(
    "mds-a9-human-substantia-nigra.R"
    , "\nMDS plot: hESC A9 and human substantia nigra"
    , "\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 20) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste(
  "../analysis/MDS - hESC A9 and human substantia nigra readDF"
  , readDepthFilt, ".pdf", sep=""))
print("#######################################################################")

# Data frame of gene symbols and corresponding Allen module color
genesColorsDF <- data.frame(colnames(exprData)
                            , bwModulesLL[[modNetworkToUse]]$colors)
colnames(genesColorsDF) <- c("gene", "module")

AddEnsembl <- function (geneList) {
  geneListDF <- data.frame(geneList)
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
                             , values = geneListDF
                             , mart = ensembl)
  moduleEnsemblDF
}

# Add Ensembl gene ID to gene symbol and module color data frame
genesEnsemblDF <- AddEnsembl(colnames(exprData))
ensemblColorsDF <- merge(genesEnsemblDF, genesColorsDF
                         , by.x = "hgnc_symbol", by.y = "gene")

# Merge with Lipton expression values
modsExprA9sNdF <- merge(ensemblColorsDF, exprA9sNdF
                       , by.x = "ensembl_gene_id", by.y = "row.names")

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

# Split by module
modsExprA9sNlDF <- split(modsExprA9sNdF, modsExprA9sNdF$module)
# Find any modules that have no genes left after merging with Lipton
sort(sapply(modsExprA9sNlDF, function(df) nrow(df) != "0"))
# Remove any modules that have no genes left after merging with Lipton
modsExprA9sNlDF <- modsExprA9sNlDF[sapply(modsExprA9sNlDF
                                          , function(df) nrow(df) != "0")]

# Call function to output data frame of MDS PCs values and PCs on each module
# group of genes
mdsLDF <- lapply(modsExprA9sNlDF, function(df) calcMDS(df[ ,c(-1,-2,-3)]))
mdsLDF <- melt(mdsLDF, id.vars = c("X1", "X2", "pc1", "pc2"))
# Add column with sample type information to MDS data frame
mdsLDF$type <- rep(c(rep("2",3), rep("7",3), rep("human",9)))
# Add column to mdsLDF to facet in ggplot2 with desired facet titles
titles <- apply(mdsLDF, 1, function(x) paste(  
            x[5]
            , "\nPC1: ", signif(100*as.numeric(x[3]), 3), "%"
            , "\nPC2: ", signif(100*as.numeric(x[4]), 3), "%", sep = ""))
mdsLDF$title <- titles

# MDS plots of Lipton data subset by each Allen module
ggplot(mdsLDF, aes(x = X1, y = X2)) +
  geom_point(aes(color = factor(type)), size = 1.5) +
  facet_wrap(~title, ncol = 8, scales = "free") +
  scale_color_discrete(name = "Sample Type"
                       , labels = c("2_ High MEF2C", "7_ Low MEF2C"
                                    , "Human Substantia Nigra")) +
  xlab("PC1") +
  ylab("PC2") +
  labs(title = paste(
    "mds-a9-human-substantia-nigra.R"
    , "\nMDS plot: hESC A9 and human substantia nigra subset by Allen modules"
    , "\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 14) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(file = paste(
  "../analysis/MDS - Allen modules hESC A9 and human substantia nigra readDF"
  , readDepthFilt, ".pdf", sep=""))

# Subset MDS data down to marker modules
mdsMarkerLDF <- mdsLDF[mdsLDF$L1 %in% modsToUse, ]
# Set factor levels to order marker modules like the order in modsToUse variable
mdsMarkerLDF$L1 <- factor(mdsMarkerLDF$L1, levels = modsToUse)
mdsMarkerLDF <- mdsMarkerLDF[order(mdsMarkerLDF$L1), ]
mdsMarkerLDF$title <- factor(mdsMarkerLDF$title
                             , levels = as.character(unique(mdsMarkerLDF$title)))

# MDS plots of Lipton data subset by each Allen marker module
ggplot(mdsMarkerLDF, aes(x = X1, y = X2)) +
  geom_point(aes(color = factor(type)), size = 2) +
  facet_wrap(~title, ncol = 4, scales = "free") +
  scale_color_discrete(name = "Sample Type"
                       , labels = c("2_ High MEF2C", "7_ Low MEF2C"
                                    , "Human Substantia Nigra")) +
  xlab("PC1") +
  ylab("PC2") +
  labs(title = paste(
    "mds-a9-human-substantia-nigra.R"
    , "\nMDS plot: hESC A9 and human substantia"
    , "\nnigra subset by Allen marker modules"
    , "\nread depth filter: ", readDepthFilt)
    , sep = "") +
  theme_grey(base_size = 15) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(file = paste(
  "../analysis/MDS - Allen marker modules hESC A9 and human substantia nigra readDF"
  , readDepthFilt, ".pdf", sep=""))
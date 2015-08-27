# MDS plots of Lipton hESC A9 neuronal cultures and human substantia nigra (SN)
# samples

sessionInfo()

library(WGCNA)
library(reshape2)
library(ggplot2)
library(biomaRt)


load("../processed_data/allen_BW_modules.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_5.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN.rda")

# variable for read depth filter to record in output graph titles
readDepthFilt <- "5"

# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 11
modsToUse <- c("plum1", "grey60", "brown", "red", "cyan", "yellowgreen"
               , "sienna3", "royalblue")
print("#######################################################################")

# MDS of all genes for hESC A9 and human SN

# Data frame of A9 samples and human SN samples, each column is a sample
exprA9sNdF <- merge(datExpr.HTSC.A9, datExpr.HTSC.SN
                    , by.x = "row.names", by.y = "row.names")
exprA9sNdF <- data.frame(exprA9sNdF[ ,-1], row.names = exprA9sNdF[ ,1])

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



plotMDS <- function (mds) {
  plot(x = mds$distVals[,1], y = mds$distVals[,2]
       , col = as.numeric(as.factor(mds$distVals[,3]))
       , pch = 16, asp=1
       , main = paste("allen-A9-marker-expr-in-hESC-A9.R\nMDS Plot By Module"
                      , " Marker Gene Expression\nmodule:", module
                      , " read depth filter:", readDepthFilt, sep = "")
       , xlab = paste("PC1 (", signif(100*mds$pc1,3), "%)", sep="")
       , ylab = paste("PC2 (", signif(100*mds$pc2,3),"%)",sep=""))
  legend("bottomright", levels(mdsAndTreatment[,3])
         , col=1:length(levels(mdsAndTreatment[,3])), pch=16, cex=0.8)
  # dev.off()
}


plotMDS(mdsLDF)


# variable for read depth filter to record in output graph titles
readDepthFilt <- "5"

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



SelectModule <- function (modNetworkToUse, modToUse) {
  blockwiseMEs <- moduleEigengenes(
    exprData, bwModulesLL[[modNetworkToUse]]$colors)$eigengenes
  
  # Write table of gene names in module - is correlation even necessary?
  # Pretty sure I copied this from WGCNA tutorial?
  geneModuleMembership <- as.data.frame(cor(exprData, blockwiseMEs, use = "p"))
  moduleGenes <- bwModulesLL[[modNetworkToUse]]$colors==modToUse
  moduleGenes <- data.frame(row.names(geneModuleMembership[moduleGenes, ]))
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

SubsetMarkerModInA9 <- function (moduleEnsemblDF) {
  # Subset Vivek produced expression data by marker module genes present
  markersExprDF <- datExpr.HTSC.A9[
    (rownames(datExpr.HTSC.A9) %in% moduleEnsemblDF[,1]), ]
  # Combine data frames to add gene symbol and module information
  markersExprDF <- merge(markersExprDF, moduleEnsemblDF
                         , by.x="row.names", by.y=1)
  # Changing column name like this is not working
  colnames(markersExprDF[1]) <- "ensembl"
  # Reshape for ggplot2 using Reshape2 package
  markersExprDF <- melt(markersExprDF
                        , value.name="expression", variable.name="sample")
  # Add column with treatment group label
  numGenesInTreatmentGroup <- (nrow(markersExprDF) / 2)
  markersExprDF <- cbind(markersExprDF
                         , bio.rep = as.factor(
                           c(rep(2,numGenesInTreatmentGroup)
                             , rep(7,numGenesInTreatmentGroup))))
  markersExprDF
}


# MDS plots

PlotMDS <- function (exprDF, module, i) {
  # mds = cmdscale(dist(markerModulesA9LDF), eig = T)
  # test <- matrix(runif(2502, 1, 10),417,6)
  # markerModulesA9LDF <- test
  mds = cmdscale(dist(t(exprDF)), eig = T)
  pc1 = mds$eig[1]^2 / sum(mds$eig^2)
  pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  # lapply((mds$eig^2 / sum(mds$eig^2)*100), function(x) signif(x,3))
  mdsAndTreatment <- data.frame(mds$points
                                , as.factor(c(rep("2",3), rep("7",3))))
  #   pdf(paste("../analysis/Allen hESC A9 PCA readDF", readDepthFilt
  #             , " ModSize30 mod-", module, "-", i, ".pdf", sep = "")
  #       , height=8, width=8)
  plot(x = mdsAndTreatment[,1], y = mdsAndTreatment[,2]
       , col = as.numeric(as.factor(mdsAndTreatment[,3]))
       , pch = 16, asp=1
       , main = paste("allen-A9-marker-expr-in-hESC-A9.R\nMDS Plot By Module"
                      , " Marker Gene Expression\nmodule:", module
                      , " read depth filter:", readDepthFilt, sep = "")
       , xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep="")
       , ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""))
  legend("bottomright", levels(mdsAndTreatment[,3])
         , col=1:length(levels(mdsAndTreatment[,3])), pch=16, cex=0.8)
  # dev.off()
}

# PCA for A9 modules identified from Allen
i <- 0
pdf(paste("../analysis/Allen hESC A9 PCA readDF", readDepthFilt
          , " minMod", minModSize, ".pdf", sep = "")
    , height=8, width=8)
for (modToUse in modsToUse) {
  # Empty variable that was used in random module PCA above
  i <- ""
  print(modToUse)
  markerModulesLDF <- SelectModule(modNetworkToUse, modToUse)
  markerModulesA9DF <- SubsetMarkerModInA9(markerModulesLDF)
  markerModulesA9DF <- dcast(markerModulesA9DF, Row.names~sample
                             , value.var = "expression")
  # Move column 1 to row names
  row.names(markerModulesA9DF) <- markerModulesA9DF[ ,1]
  markerModulesA9DF <- markerModulesA9DF[ ,-1]
  PlotMDS(markerModulesA9DF, modToUse, i)
}
dev.off()
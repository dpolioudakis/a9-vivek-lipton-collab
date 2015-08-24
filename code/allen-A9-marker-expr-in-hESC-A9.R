# Expression of A9 module 28 genes in hESC derived A9 neuronal cultures

print("#######################################################################")
print("Starting allen-A9-marker-expr-in-hESC-A9.R script...")
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
# Vivek normalized RNAseq FPKMs from hESC derived A9 neuronal cultures
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_0.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_5.rda")

# variable for read depth filter to record in output graph titles
readDepthFilt <- "5"
minModSize <- "100"

# bwModulesLL is list of modules from different blockwiseModules parameters used
# 12 corresponds to softPower 9, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 12
modsToUse <- c("saddlebrown", "salmon", "red", "pink", "black", "green"
               , "grey60")

modNetworkToUse <- 11
modsToUse <- c("plum1", "grey60", "brown", "red", "cyan", "yellowgreen"
               , "sienna3", "royalblue")

modNetworkToUse <- 18
modsToUse <- c("salmon", "brown", "blue", "purple", "greenyellow")
print("#######################################################################")

# Two general functions to load:

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
print("#######################################################################")

# # Plot expression marker genes high MEF2C and low MEF2C samples
# # Separate plot for each module
# for (modToUse in modsToUse) {
#   print(modToUse)
#   markerModulesLDF <- SelectModule(modNetworkToUse, modToUse)
#   markerModulesA9LDF <- SubsetMarkerModInA9(markerModulesLDF)
#   print(head(markerModulesA9LDF))
#   # Plot expression marker genes high MEF2C and low MEF2C samples
#   quartz()
#   ggplot(markerModulesA9LDF, aes(x=hgnc_symbol, y=expression)) + 
#     geom_boxplot(aes(fill=bio.rep)) +
#     scale_fill_discrete(name= "Biological\nReplicate",
#                         labels= c("2_(high MEF2C)", "7_(low MEF2C)")) +
#     labs(title = paste(
#      "Allen derived A9 marker expression in Lipton A9\nmodule:", modToUse)) +
#     ylab("Expression (normalized FPKM)") +
#     xlab("Gene") +
#     theme_grey(base_size = 18)
#   ggsave(file= paste(
#     "../analysis/Allen A9 marker expression - minModSize 30 - module -", modToUse, Sys.Date(), ".pdf"))
# }
# print("#######################################################################")

# PCA plots

RandomModule <- function (modNetworkToUse, modSizes) {
  blockwiseMEs <- moduleEigengenes(
    exprData, bwModulesLL[[modNetworkToUse]]$colors)$eigengenes
  geneModuleMembership <- as.data.frame(cor(exprData, blockwiseMEs, use = "p"))
  moduleGenes <- c(sample(
      1:length(bwModulesLL[[modNetworkToUse]]$colors)
      , size = modSizes, replace=F))
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

# PCA for some randomly drawn modules as well as all genes in the network
modSizes <- c(30, 30, 30, 100, 100, 100
  , length(bwModulesLL[[modNetworkToUse]]$colors))
i <- 0
pdf(paste("../analysis/Allen hESC A9 PCA readDF", readDepthFilt
          , " ", minModSize, ".pdf", sep = "")
    , height=8, width=8)
for (modSize in modSizes) {
  print(modSize)
  # Counter used to save graphs with different names so as not to overwrite
  i <- i + 1
  print(i)
  randomModuleDF <- RandomModule(modNetworkToUse, modSize)
  randomModuleDF <- SubsetMarkerModInA9(randomModuleDF)
  randomModuleDF <- dcast(randomModuleDF, Row.names~sample
                          , value.var = "expression")
  # Move column 1 to row names
  row.names(randomModuleDF) <- randomModuleDF[ ,1]
  randomModuleDF <- randomModuleDF[ ,-1]
  PlotMDS(randomModuleDF, modSize, i)
}

# PCA for A9 modules identified from Allen
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
print("#######################################################################")

# Mean of expression fold changes for each module marker gene in high MEF2C
# versus low MEF2C
# Make list of data frames of expression (normalized FPKM) for each gene
# Each list element is a module
markerModulesA9LDF <- NULL
for (modToUse in modsToUse) {
  print(modToUse)
  markerModulesLDF <- SelectModule(modNetworkToUse, modToUse)
  # markerModulesLDF <- SelectModule(modNetworkToUse, "plum1")
  # Subset genes in module to only those found in Lipton hESC A9 data
  markerModulesA9LDF[[modToUse]] <-SubsetMarkerModInA9(markerModulesLDF)
}

# Calculate ratio of expression in high MEF2C samples versus low MEF2C samples
# for each gene in each module
# List (modules) of lists (each gene in that module)
ratioExprLL <- lapply(markerModulesA9LDF, 
                      function(x) {
                        markExpr <- dcast(x, Row.names~sample, value.var="expression")
                        row.names(markExpr) <- markExpr[ ,1]
                        markExpr <- markExpr[ ,-1]
                        apply(markExpr, 1, function(x) (sum(x[1:3]) / sum(x[4:6])))
                      }
)
# Mean ratio of expression in each module
sapply(ratioExprLL, mean)
# Paired T-test of high MEF2C versus low MEF2C expression (not ratios) for each
# module
sapply(markerModulesA9LDF, function(module) {
  highMEF2C <- module[module$bio.rep=="2", ]$expression
  lowMEF2C <-  module[module$bio.rep=="7", ]$expression
  round(t.test(highMEF2C, lowMEF2C, paired = TRUE)$p.value, 5)
})

# Boxplot of expression ratios for each module
ratioExprDF <- melt(ratioExprLL)
colnames(ratioExprDF) <- c("ratio.expr", "module")
ggplot(data = ratioExprDF, aes(x=module, y=ratio.expr)) + 
  geom_boxplot() +
  # geom_boxplot(aes(fill=module)) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(title = paste(
    "allen-A9-marker-expr-in-hESC-A9.R\nAllen derived A9 marker expression in"
    , "Lipton A9\nread depth filter: ", readDepthFilt)
    , sep = "") +
  ylab("Mean Expression (normalized FPKM)") +
  xlab("Treatment") +
  theme_grey(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  ggsave(file = paste(
    "../analysis/Allen hESC A9 ratio expr readDF", readDepthFilt
    , " ModSize30 mod-", modToUse, "-", Sys.Date(), ".pdf", sep=""))
print("#######################################################################")

# Mean expression
# Loop through models of interest in network
for (modToUse in modsToUse) {
  print(modToUse)
  markerModulesLDF <- SelectModule(modNetworkToUse, modToUse)
  markerModulesLDF <- SelectModule(modNetworkToUse, "plum1")
  # Subset genes in module to only those found in Lipton hESC A9 data
  markerModulesA9LDF <- SubsetMarkerModInA9(markerModulesLDF)
  # Make data frame of expression for each Tx group as col 1 and col 2
  txMarkersMeansDF <- data.frame(
      markerModulesA9LDF[markerModulesA9LDF$bio.rep == 2, ]$expression
    , markerModulesA9LDF[markerModulesA9LDF$bio.rep == 7, ]$expression)
  colnames(txMarkersMeansDF) <- c("Tx_2_highMEF2C", "Tx_7_lowMEF2C")
  # Paired T-test comparing high MEF2C Tx to low MEF2C Tx
  pval <- t.test(txMarkersMeansDF$Tx_2_highMEF2C, txMarkersMeansDF$Tx_7_lowMEF2C
                 , paired=TRUE)$p.value
  # Format expression data for ggplot2
  txMarkersMeansDF <- melt(txMarkersMeansDF
                           , measure.vars=c("Tx_2_highMEF2C", "Tx_7_lowMEF2C"))
  colnames(txMarkersMeansDF) <- c("treatment", "expression")
  print(txMarkersMeansDF)
  
  # Make dataframe of pval to keep ggplot2 happy
  formattedPvalggplotDF <- data.frame(
      label = paste("Paired T-test p-value:", as.character(signif(pval, 3)))
  )
  print(ggplot(data = txMarkersMeansDF, aes(x = treatment, y = expression)) + 
    geom_boxplot(aes(fill = treatment)) +
    geom_text(data = formattedPvalggplotDF, aes(1.5, 8, label = label), type = "NA*") +
    scale_fill_discrete(name = "Biological\nReplicate",
                        labels = c("Tx_2_(high MEF2C)", "Tx_7_(low MEF2C)")) +
    labs(title = paste(
      "allen-A9-marker-expr-in-hESC-A9.R\nAllen derived A9 marker expression in"
      , "Lipton A9\nmodule: ", modToUse, sep = "")) +
    ylab("Mean Expression (normalized FPKM)") +
    xlab("Treatment") +
    theme_bw(base_size = 18) +
    ggsave(file = paste(
      "../analysis/Allen A9 marker expr in hESC A9 minModSize30 mod-"
      , modToUse, "-", Sys.Date(), ".pdf", sep = ""))
  )
}
print("#######################################################################")
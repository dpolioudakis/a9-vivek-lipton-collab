# Expression of Allen modules in Lipton Human SN versus all genes in Lipton
# Human SN

library(biomaRt)
library(reshape2)
library(ggplot2)

load("../Vivek_WGCNA_Lipton_A9_SN/Vivek_SN_mod_colors.rda")
load("../processed_data/allen_BW_modules.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN.rda")

# Variable for minModSize parameter used in blockwiseModules WGCNA function
# to record in output graph titles
minModSize <- "30"
# bwModulesLL is list of modules from different blockwiseModules parameters used
# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
# 18corresponds to softPower 9, minModSize 100, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
# modNetworkToUse <- 11
modNetworkToUse <- 11

# Split by model into lists or dataframes of genes in each model
allenGenesModelsLL <- split(colnames(exprData)
                          , bwModulesLL[[modNetworkToUse]]$colors)

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
geneListsLDF <- lapply(allenGenesModelsLL, AddEnsembl)

# Add expression data to Ensembl ID and Gene symbol module lists
humanSNexprLL <- lapply(geneListsLDF
                        , function(geneList) merge(
                          x = geneList
                          , y = datExpr.HTSC.SN
                          , by.x = "ensembl_gene_id"
                          , by.y = "row.names"))
print("#######################################################################")

# Calculate ratio of mean expression of genes in module versus mean expression
# of all genes profiled
ratiosExprLL <- sapply(humanSNexprLL, function(moduleExpr)
  mean(as.matrix(moduleExpr[ , c(-1,-2)])) / mean(datExpr.HTSC.SN))

# Order ratios of mean expression by least to greatest
ratiosExprDF <- data.frame(sort(ratiosExprLL))
ratiosExprDF <- cbind(ratiosExprDF, row.names(ratiosExprDF))
colnames(ratiosExprDF) <- c("ratios", "modules")

# Set factor levels to maintain order of ratios in ggplot2
ratiosExprDF <- within(ratiosExprDF,
                       modules <- factor(modules
                             , levels = modules))

ggplot(ratiosExprDF, aes(x = modules, y = ratios)) +
     geom_bar(stat="identity") +
     geom_hline(aes(yintercept = 1)) +
     labs(title = paste(
          "allen-module-expr-in-human-SN.R"
          ,"\nA9 marker module ratio of expression in Lipton human substantia"
          ,"\nnigra versus all genes profiled in human substantia nigra"
          ,"\n minModSize", minModSize
          , sep = "")) +
     theme_grey(base_size = 21) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     ylab("Ratio Expression of Module versus All Genes") +
     xlab("Modules") +
     theme(axis.text = element_text(color = "black"))
     ggsave(file = paste(
         "../analysis/Allen module expression in human substantia nigra"
       , "minModSize", minModSize, ".pdf", sep=""))
print("#######################################################################")

# Calculate ratio of expression of each gene in module versus mean expression
# of all genes profiled
ratiosExprLL <- sapply(humanSNexprLL, function(moduleExpr) {
  ratiosExprLL <- apply(as.matrix(moduleExpr[c(-1,-2)]), 1, function(geneExpr)
    (mean(geneExpr) / mean(datExpr.HTSC.SN))
  )
})

ratiosExprDF <- melt(ratiosExprLL)
colnames(ratiosExprDF) <- c("ratios", "modules")
# Order ratios of mean expression by least to greatest

ggplot(ratiosExprDF, aes(x = modules, y = ratios)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



ratioExprDF$module <- factor(ratioExprDF$module
                             , levels = as.character(unique(ratioExprDF$module)))
ggplot(data = ratioExprDF, aes(x=module, y=ratio.expr)) + 
  geom_boxplot() +
  scale_x_discrete(labels = c(  "30 Gene Marker"
                                , "100 Gene Marker"
                                , "30 Gene Anti-marker"
                                , "100 Gene Anti-marker")) +
  # geom_boxplot(aes(fill=module)) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(title = paste(
    "multiTOM-module-expr-A9.R
    Mean ratio of expression of multiTOM module genes in hESC A9
    high MEF2C samples versus low MEF2C samples"
    , "\nMarker Module Seeds: ALDH1A1, TH, SLC18A2, KCNJ6"
    , "\nAnti-marker Module Seeds: CALB1"
    , sep = "")) +
  ylab("Mean Expression Ratio") +
  xlab("multiTOM module") +
  theme_grey(base_size = 20) +
  theme(axis.text = element_text(color = "black")) +
  ggsave(file ="../analysis/multiTOM module ratio expression in hESC A9.pdf")
print("#######################################################################")

# Developmental script examing Allen module eigengene expression in Lipton's 
# human Substantia Nigra samples

load("../processed_data/allen_BW_modules.rda")
load("../processed_data/array_data_subset_avg_probes.rda")

# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 11
modsToUse <- c("plum1", "grey60", "brown", "red", "cyan", "yellowgreen"
               , "sienna3", "royalblue")

genesColorsDF <- data.frame(colnames(exprData), bwModulesLL[[modNetworkToUse]]$colors)
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

exprSNDF <- data.frame(datExpr.HTSC.SN)
modulesA9 <- merge(ensemblColorsDF, exprSNDF, by.x = "ensembl_gene_id", by.y = "row.names" )

allenMEinSN <- moduleEigengenes(t(modulesA9[ ,4:12]), modulesA9$module)$eigengenes
allenMEinSN$biorep <- c(rep(2, 3), rep(7, 3))
allenMEinSN <- melt(allenMEinSN)
colnames(allenMEinSN) <- c("module", "MEexpression")

MEtoUse <- sapply(modsToUse, function(mod) paste("ME", mod, sep=""))
markerMEinSN <- allenMEinSN[allenMEinSN$module %in% MEtoUse, ]

markerMEinSN$module <- factor(markerMEinSN$module
                              , levels = as.character(unique(MEtoUse)))


ggplot(data = markerMEinSN, aes(x = module, y = MEexpression)) +
  geom_boxplot() +
  scale_fill_discrete(name= "Biological\nReplicate",
                      labels= c("2_HighMEF2C", "7_LowMEF2C")) +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(color = "black"))


ggplot(data = allenMEinSN, aes(x = module, y = MEexpression)) +
  geom_boxplot() +
  scale_fill_discrete(name= "Biological\nReplicate",
                      labels= c("2_HighMEF2C", "7_LowMEF2C")) +
  theme_grey(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(color = "black"))
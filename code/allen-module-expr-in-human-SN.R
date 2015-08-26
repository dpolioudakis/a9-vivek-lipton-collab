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
modNetworkToUse <- 18

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

humanSNexprLL <- lapply(geneListsLDF
                        , function(geneList) merge(
                            x = geneList
                          , y = datExpr.HTSC.SN
                          , by.x = "ensembl_gene_id"
                          , by.y = "row.names"))
ratiosExprLL <- sapply(humanSNexprLL, function(moduleExpr)
  mean(as.matrix(moduleExpr[ , c(-1,-2)])) / mean(datExpr.HTSC.SN))

ratiosExprDF <- data.frame(sort(ratiosExprLL))
ratiosExprDF <- cbind(ratiosExprDF, row.names(ratiosExprDF))
colnames(ratiosExprDF) <- c("ratios", "modules")

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
          ,"\n minModSize30)"
          , sep = "")) +
     theme_grey(base_size = 21) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     ylab("Ratio Expression of Module versus All Genes") +
     xlab("Modules") +
     theme(axis.text = element_text(color = "black"))
     ggsave(file = paste(
     "../analysis/Allen module expression in human substantia nigra.pdf", sep=""))
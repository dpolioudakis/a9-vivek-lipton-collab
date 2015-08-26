# multiTOM module expression in Lipton hESC A9 cells

library(reshape2)
library(ggplot2)
library(biomaRt)

# Vivek normalized RNAseq FPKMs from hESC derived A9 neuronal cultures
# load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells.rda")
# load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_0.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_5.rda")
dataMTmark100DF <- read.csv(
  "../processed_data/multiTOM_allen_neighbors100_ALDH1A1_TH_SLC18A2_KCNJ6.csv")
dataMTmark30DF <- read.csv(
  "../processed_data/multiTOM_allen_neighbors30_ALDH1A1_TH_SLC18A2_KCNJ6.csv")

dataMTLDF <- list(  dataMTmark30DF = dataMTmark30DF[ ,1]
                  , dataMTmark100DF = dataMTmark100DF[, 1])

print("#######################################################################")

# Look up Ensembl IDs
AddEnsembl <- function (moduleGenes) {
  # bioMart manual:
  # http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  moduleEnsembDF <- getBM(  attributes = c("ensembl_gene_id", "hgnc_symbol")
                             , filters = "hgnc_symbol"
                             , values = moduleGenes
                             , mart = ensembl)
  moduleEnsembDF
}

ensemblMTLDF <- lapply(dataMTLDF, AddEnsembl)

# Subset hESC A9 data by multiTOM modules
  # Format data frames to have column names:
    #Row.names hgnc_symbol #sample #expression #bio.rep
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

MTmodA9exprLDF <- lapply(ensemblMTLDF, SubsetMarkerModInA9)


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
ratioExprLL <- lapply(MTmodA9exprLDF, 
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
# Preserve module order in boxplot
ratioExprDF$module <- factor(ratioExprDF$module
                             , levels = as.character(unique(ratioExprDF$module)))
ggplot(data = ratioExprDF, aes(x=module, y=ratio.expr)) + 
  geom_boxplot() +
  # geom_boxplot(aes(fill=module)) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(title = paste(
    "multiTOM-module-expr-A9.R
    Mean ratio of expression of muliTOM module genes in hESC A9
    high MEF2C samples versus low MEF2C samples")
    , sep = "") +
  ylab("Mean Expression Ratio") +
  xlab("multiTOM module") +
  theme_grey(base_size = 20) +
  theme(axis.text = element_text(color = "black")) +
  ggsave(file ="../analysis/multiTOM module ratio expression in hESC A9.pdf")  #  "-", Sys.Date(),
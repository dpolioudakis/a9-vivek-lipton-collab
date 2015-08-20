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

# bwModulesLL is list of modules from different blockwiseModules parameters used
# 12 corresponds to softPower 9, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 12
modToUse <- c("saddlebrown", "salmon", "red", "pink", "black", "green", "grey60")
blockwiseMEs <- moduleEigengenes(
  exprData, bwModulesLL[[modNetworkToUse]]$colors)$eigengenes

# Write table of gene names in module 28 (ME correlated with A9 markers)
geneModuleMembership <- as.data.frame(cor(exprData, blockwiseMEs, use = "p"))
moduleGenes <- bwModulesLL[[modNetworkToUse]]$colors=="saddlebrown"
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

# Plot expression marker genes high MEF2C and low MEF2C samples
ggplot(markersExprDF, aes(x=hgnc_symbol, y=expression)) + 
  geom_boxplot(aes(fill=bio.rep)) +
  scale_fill_discrete(name= "Biological\nReplicate",
                      labels= c("2_(high MEF2C)", "7_(low MEF2C)")) +
  labs(title= "Allen derived A9 marker expression") +
  ylab("Expression (normalized FPKM)") +
  xlab("Gene") +
  theme_grey(base_size = 18)
ggsave(file= "A9 marker gene expression.pdf")

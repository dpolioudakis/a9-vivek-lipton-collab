# To do:
# Cannot output dendrogram of combined blockwise modules, can only access
# dendrograms of the separate blocks?



# Cluster samples and construct modules blockwise

print("#######################################################################")
print("Starting allen-construct-network-modules-blockwise.R script...")
sessionInfo()

library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/array_data_subset_avg_probes.rda")

# Transpose expression data and assign column names from column 1
exprData <- t(arrayDataSubsetAvgProbesDF)
# # Selecting 5000 genes with highest expression values (avg across samples)
# exprDataTop5000 <- exprData[,rank(-colMeans(exprData))<=5000]
# exprData <- exprDataTop5000
# # Sampling 5000 genes randomly
# exprDataRandom5000 <- exprData[ , sample(ncol(exprData), 5000)]
# exprData <- exprDataRandom5000

softPower = 5

# Module identification using blockwise:

# Function to construct modules
# Args: minimum module size, deep split
MakeBWModules <- function (minModSize, deepSplit, mergeCutHeight) {
  print("Treecut arguments:")
  print(c(minModSize, deepSplit, mergeCutHeight))
  bwNet = blockwiseModules(exprData, maxBlockSize = 10000, power = softPower
                           , TOMType = "signed", minModuleSize = minModSize
                           , reassignThreshold = 0, pamStage = FALSE,
                           , numericLabels = TRUE, deepSplit = deepSplit
                           , mergeCutHeight = mergeCutHeight
                           ,  verbose = 3)
  print("Table of genes per module:")
  print(table(bwNet$colors))
  bwNet
}


minModSize <- c(30)
deepSplit <- c(2)
# Cut height of 0.25, corresponds to a correlation of 0.75, to merge ME:
MEmergeCutHeight <- c(0.2)
bwModules <- MakeBWModules(minModSize, deepSplit, MEmergeCutHeight)

save(bwModules, exprData, file="../processed_data/allen_BW_modules.rda")



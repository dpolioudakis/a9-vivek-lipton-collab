# To do:
# Cannot output dendrogram of combined blockwise modules, can only access
# dendrograms of the separate blocks?


# Workflow
#   1a-allen-subset-to-basal-ganglia.R
#   1b-allen-combine-probes.R
#   2a-allen-soft-thresholding-power.R
#   2b-allen-adjacency-TOM.R
#   3-allen-construct-network-modules.R
#   4-allen-compare-modules-metadata.R


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
  bwNet = blockwiseModules(exprData, maxBlockSize = 12000, power = softPower
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
MEmergeCutHeight <- c(0.2, 0.25, 0.3)
# bwModules <- MakeBWModules(minModSize, deepSplit, MEmergeCutHeight)
# blockwiseMEs <- moduleEigengenes(exprData, bwModules$colors)$eigengenes

bwModules <- mapply(MakeBWModules
                    , MEmergeCutHeight
                    , MoreArgs = list(
                      minModSize = minModSize, deepSplit = deepSplit)
                    , SIMPLIFY = FALSE)

blockwiseMEs <- lapply(bwModules, function (x)
                       moduleEigengenes(exprData, x$colors)$eigengenes)


# signif(cor(blockwiseMEs, blockwiseMEs), 3)

pdf("../analysis/3a_BW_ME_correlation.pdf", height=8, width=8)
par(cex = 1.0)
# sizeGrWindow(6,7)
par(mfrow = c(3,1))
for (i in 1:length(blockwiseMEs)) {
  plotEigengeneNetworks(orderMEs(blockwiseMEs[[i]])
    , paste("Eigengene correlation - merge cut height:", MEmergeCutHeight[i])
    , signed=TRUE, colorLabels=TRUE,
    , marHeatmap = c(1,4,1,2)
    , marDendro = c(0,4,2,0))
}
dev.off()

save(bwModules, exprData, file="../processed_data/allen_BW_modules.rda")



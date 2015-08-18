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

softPower = c(5,7,9)
minModSize <- c(30, 50, 100)
deepSplit <- c(2)
# Cut height of 0.25, corresponds to a correlation of 0.75, to merge ME:
MEmergeCutHeight <- c(0.2, 0.25, 0.3)
maxBlockSize <- c(12000)
ArgsbwModulesLLDF <- expand.grid(softPower = softPower
                             , minModSize = minModSize
                             , deepSplit = deepSplit
                             , MEmergeCutHeight = MEmergeCutHeight
                             , maxBlockSize = maxBlockSize)


load("../processed_data/array_data_subset_avg_probes.rda")

# Transpose expression data and assign column names from column 1
exprData <- t(arrayDataSubsetAvgProbesDF)
# # Selecting 5000 genes with highest expression values (avg across samples)
# exprDataTop5000 <- exprData[,rank(-colMeans(exprData))<=5000]
# exprData <- exprDataTop5000
# # Sampling 5000 genes randomly
# exprDataRandom5000 <- exprData[ , sample(ncol(exprData), 5000)]
# exprData <- exprDataRandom5000



# Module identification using blockwise:

# Function to construct modules
# blockwiseModules identifies modules using Dynamic Hybrid tree cut
# Args: minimum module size, deep split
MakebwModulesLL <- function (softPower
                           , minModSize
                           , deepSplit
                           , MEmergeCutHeight
                           , maxBlockSize) {
  print("WGCNA blockwiseModules function arguments:")
  print(c(softPower, minModSize, deepSplit, MEmergeCutHeight, maxBlockSize))
  bwNet <- blockwiseModules(  exprData
                           , power = softPower
                           , minModuleSize = minModSize
                           , deepSplit = deepSplit
                           , mergeCutHeight = MEmergeCutHeight
                           , maxBlockSize = maxBlockSize
                           , TOMType = "signed"
                           , reassignThreshold = 0, pamStage = FALSE
                           , numericLabels = TRUE
                           , verbose = 3)
  print("Table of genes per module:")
  print(table(bwNet$colors))
  bwNet
}

# ArgsbwModulesLLDF is all combinations of arguments for MakebwModulesLL
# Each row of ArgsbwModulesLLDF is a set of arguments
# bwModuleLL is a list, each item is an output of MakebwModulesLL that is
# calling the WGCNA blockwiseModules function.  blockwiseModules output is a
# list
bwModulesLL <- apply(ArgsbwModulesLLDF, 1
                     , function(x) MakebwModulesLL(x[1],x[2],x[3],x[4],x[5]))

# DELETE THIS BLOCK IF SCRIPT UPDATE WORKS
# # bwModulesLL <- MakebwModulesLL(minModSize, deepSplit, MEmergeCutHeight)
# # blockwiseMEs <- moduleEigengenes(exprData, bwModulesLL$colors)$eigengenes
# bwModulesLL <- mapply(MakebwModulesLL
#                     , MEmergeCutHeight
#                     , MoreArgs = list(
#                       minModSize = minModSize, deepSplit = deepSplit)
#                     , SIMPLIFY = FALSE)



# Graph module eigengene correlations for each set of arguments used in
# blockwiseModules WGCNA function
blockwiseMEs <- lapply(bwModulesLL, function (x)
  moduleEigengenes(exprData, x$colors)$eigengenes)
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

save(bwModulesLL, exprData, file="../processed_data/allen_BW_modules.rda")



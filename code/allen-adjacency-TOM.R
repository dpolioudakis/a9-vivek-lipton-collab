# Cluster samples and construct modules

# Workflow
#   1a-allen-subset-to-basal-ganglia.R
#   1b-allen-combine-probes.R
#   2a-allen-soft-thresholding-power.R
#   2b-allen-adjacency-TOM.R
#   3-allen-construct-network-modules.R
#   4-allen-compare-modules-metadata.R

print("#######################################################################")
print("Starting allen-adjacency-TOM.R script...")
sessionInfo()

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/array_data_subset_avg_probes.rda")

# Transpose expression data and assign column names from column 1
exprData <- t(arrayDataSubsetAvgProbesDF)
# colnames(exprData) <- exprData[1, ]
# exprData <- exprData[-1, ]
# Selecting 5000 genes with highest expression values (avg across samples)
# exprDataTop5000 <- exprData[,rank(-colMeans(exprData))<=5000]
# exprData= exprDataTop5000
# Selecting 5000 random genes
# exprDataRandom5000 <- exprData[ , sample(ncol(exprData), 5000)]
# exprData <- exprDataRandom5000

softPower = 5
# Biweight midcorrelation is considered to be a good alternative to Pearson
# correlation since it is more robust to outliers.
adjacency = adjacency(exprData, power= softPower, corFnc= "bicor")

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

save(adjacency, TOM, dissTOM
     , file="../processed_data/allen_adjacency_TOM.rda")
# Cluster samples and construct modules blockwise

print("#######################################################################")
print("Starting allen-ME-region-corr.R script...")
sessionInfo()

library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/allen_BW_modules.rda")
load("../processed_data/array_data_subset_avg_probes.rda")

blockwiseMEs <- moduleEigengenes(exprData, bwModules$colors)$eigengenes

signif(cor(blockwiseMEs, blockwiseMEs), 3)

pdf("../analysis/3b_BW_ME_correlation.pdf", height=8, width=8)
par(cex = 1.0)
sizeGrWindow(6,7)
plotEigengeneNetworks(orderMEs(blockwiseMEs), "Eigengene correlation"
                               , signed=TRUE, colorLabels=TRUE,
                               , marHeatmap = c(1,4,1,2)
                               , marDendro = c(0,4,2,0))
dev.off()

# Each list of metaDataSubsetLDF is a DF of metadata corresponding to each brain
# Brains and samples are listed in the order of the observations in blockwiseMEs
# Make vector of structure acronyms in order of blockwiseMEs
brainRegionV <- NULL
brainRegionV <- lapply(metaDataSubsetLDF
                       , function(x) c(brainRegionV, x$structure_acronym))
brainRegionV <- unlist(brainRegionV)
brainRegionV <- as.factor(brainRegionV)
# Boxplots of the ME expression by brain region for each ME
sizeGrWindow(12,12)
pdf("../analysis/3b_BW_ME_region_corr.pdf", height=12, width=12)
par(mfrow = c(4,5))
par(las=2)
# Make list of DFs
#   Col 1: ME expression
#   Col 2: brain region
#   Rows: Samples
MEbrainRegionLDF <- lapply(blockwiseMEs
                           , function(ME) data.frame(ME, brainRegionV))
MEnames <- names(MEbrainRegionLDF)
# Loop through list of ME expression and brain region and list of ME name
# and plot
for(i in 1:length(MEbrainRegionLDF)) {
  boxplot(ME~brainRegionV, data=MEbrainRegionLDF[[i]], main=MEnames[[i]])
  }
dev.off()
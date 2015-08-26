# Cluster samples and construct modules blockwise

# Workflow
#   1a-allen-subset-to-basal-ganglia.R
#   1b-allen-combine-probes.R
#   2a-allen-soft-thresholding-power.R
#   2b-allen-adjacency-TOM.R
#   3-allen-construct-network-modules.R
#   4-allen-compare-modules-metadata.R


print("#######################################################################")
print("Starting allen-ME-expr-by-brain-region.R script...")
sessionInfo()

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/allen_BW_modules.rda")
load("../processed_data/array_data_subset_avg_probes.rda")

# bwModules is list of modules from 3 different ME merge cut heights
blockwiseMEs <- moduleEigengenes(exprData, bwModules[[3]]$colors)$eigengenes

# Each list of metaDataSubsetLDF is a DF of metadata corresponding to each brain
# Brains and samples are listed in the order of the observations in blockwiseMEs
# Make vector of structure acronyms in order of blockwiseMEs
brainRegionV <- NULL
brainRegionV <- lapply(metaDataSubsetLDF
                       , function(x) c(brainRegionV, as.character(x$structure_acronym)))
brainRegionV <- unlist(brainRegionV)
brainRegionV <- as.factor(brainRegionV)
# Boxplots of the ME expression by brain region for each ME
# sizeGrWindow(12,12)
# par(mfrow = c(12,3))
# par(las=2)
# Make list of DFs
#   Col 1: ME expression
#   Col 2: brain region
#   Rows: Samples
MEbrainRegionLDF <- lapply(blockwiseMEs
                           , function(ME) data.frame(ME, brainRegionV))
MEnames <- names(MEbrainRegionLDF)
# Loop through list of ME expression and brain region and list of ME name
# and plot separate graph for each ME
for(i in 1:length(MEbrainRegionLDF)) {
  pdf(paste("../analysis/3b_BW_ME_region_corr.pdf",i))
  boxplot(ME~brainRegionV, data=MEbrainRegionLDF[[i]], main=MEnames[[i]]
          , ylab = "ME Expression (arbitrary value)"
          , xlab = "Brain region")
  dev.off()
}







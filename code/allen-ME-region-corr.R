# Cluster samples and construct modules blockwise

# Workflow
#   1a-allen-subset-to-basal-ganglia.R
#   1b-allen-combine-probes.R
#   2a-allen-soft-thresholding-power.R
#   2b-allen-adjacency-TOM.R
#   3-allen-construct-network-modules.R
#   4-allen-compare-modules-metadata.R


print("#######################################################################")
print("Starting allen-ME-region-corr.R script...")
sessionInfo()

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/allen_BW_modules.rda")
load("../processed_data/array_data_subset_avg_probes.rda")

# bwModules is list of modules from 3 different ME merge cut heights
blockwiseMEs <- moduleEigengenes(exprData, bwModules[[3]]$colors)$eigengenes


pdf("../analysis/3b_BW_ME_correlation.pdf", height=8, width=8)
par(cex = 1.0)
sizeGrWindow(6,7)
plotEigengeneNetworks(orderMEs(blockwiseMEs), "Eigengene correlation"
                               , signed=TRUE, colorLabels=TRUE,
                               , marHeatmap = c(1,4,1,2)
                               , marDendro = c(0,4,2,0))
dev.off()







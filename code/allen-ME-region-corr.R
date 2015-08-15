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
                       , function(x) c(brainRegionV, as.character(x$structure_acronym)))
brainRegionV <- unlist(brainRegionV)
brainRegionV <- as.factor(brainRegionV)
# Boxplots of the ME expression by brain region for each ME
sizeGrWindow(12,12)
pdf("../analysis/3b_BW_ME_region_corr.pdf", height=12, width=12)
par(mfrow = c(4,4))
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
for(i in 1:16) {
  boxplot(ME~brainRegionV, data=MEbrainRegionLDF[[i]], main=MEnames[[i]])
}
dev.off()


markerMEcor <- cor(bwModules$MEs, t(arrayDataSubsetAvgProbesDF[c("ALDH1A1", "TH", "SLC18A2", "KCNJ6", "CALB1"), ]))

markerMEcor <- cor(bwModules$MEs, t(arrayDataSubsetAvgProbesDF[c("ALDH1A1", "TH", "SLC18A2", "CACNA1D", "CALB1", "CALB2", "KCNJ6", "LMX1A", "FOXA2", "NR4A2", "ALDH1A1"), ]))

labeledHeatmap(markerMEcor
               , colnames(markerMEcor)
               , rownames(markerMEcor)
               , colorLabels = FALSE
               , colors=greenWhiteRed(50)
               , setStdMargins = TRUE
               , textMatrix = signif(markerMEcor,2)
               , cex.text = 0.5
               , cex.lab = 0.5
               , zlim = c(-1,1)
               , main = "Gene-module correlation")




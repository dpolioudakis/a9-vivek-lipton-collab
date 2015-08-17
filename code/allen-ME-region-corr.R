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


# Correlate ME expression for each sample to marker gene expression
# CALB1 and CALB2 are anti-markers
markerMEcor <- cor(bwModules[[3]]$MEs, t(arrayDataSubsetAvgProbesDF[
  c("ALDH1A1", "TH", "SLC18A2", "KCNJ6", "CALB1"), ]))

# Correlate with all markers listed in Vivek's work
markerMEcor <- cor(bwModules[[3]]$MEs, t(arrayDataSubsetAvgProbesDF[
  c("ALDH1A1", "TH", "SLC18A2", "CACNA1D", "CALB1"
    , "CALB2", "KCNJ6", "LMX1A", "FOXA2", "NR4A2", "ALDH1A1"), ]))

pdf("../analysis/3b_MEs_markers_corr.pdf", height=8, width=8)
labeledHeatmap(markerMEcor
               , colnames(markerMEcor)
               , rownames(markerMEcor)
               , colorLabels = FALSE
               , colors=greenWhiteRed(50)
               , setStdMargins = TRUE
               , textMatrix = signif(markerMEcor,2)
               , cex.text = 1
               , cex.lab = 1
               , zlim = c(-1,1)
               , main = "Gene to module eigenegene correlation")
dev.off()




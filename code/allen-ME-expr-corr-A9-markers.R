print("#######################################################################")
print("Starting allen-ME-expr-corr-A9-markers.R script...")
sessionInfo()

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/allen_BW_modules.rda")
load("../processed_data/array_data_subset_avg_probes.rda")

# bwModules is list of modules from 3 different ME merge cut heights
blockwiseMEs <- moduleEigengenes(exprData, bwModules[[3]]$colors)$eigengenes

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
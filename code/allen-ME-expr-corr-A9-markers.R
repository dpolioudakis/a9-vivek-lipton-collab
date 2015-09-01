# Correlate A9 marker genes to all module groups produced by each run of WGCNA
# blockwise module function with different parameters

print("#######################################################################")
print("Starting allen-ME-expr-corr-A9-markers.R script...")
sessionInfo()

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/allen_BW_modules.rda")
load("../processed_data/array_data_subset_avg_probes.rda")

# bwModulesLL is list of modules from different blockwiseModules parameters used

# Correlate ME expression for each sample to marker gene expression
# CALB1 and CALB2 are anti-markers
markerMEcor <- lapply(bwModulesLL, function(x)
  cor(x$MEs, t(arrayDataSubsetAvgProbesDF[
  c("ALDH1A1", "TH", "SLC18A2", "KCNJ6", "CALB1"), ])))

# Correlate with all markers listed in Vivek's work
markerMEcor <- lapply(bwModulesLL, function(x)
  cor(x$MEs, t(arrayDataSubsetAvgProbesDF[
  c("ALDH1A1", "TH", "SLC18A2", "CACNA1D", "CALB1"
    , "CALB2", "KCNJ6", "LMX1A", "FOXA2", "NR4A2", "ALDH1A1"), ])))

# # May need to run this code if ArgsbwModulesLLDF was not saved in
# # allen_BW_modules.rda
# softPower = c(5,7,9)
# minModSize <- c(30, 50, 100)
# deepSplit <- c(2)
# # Cut height of 0.25, corresponds to a correlation of 0.75, to merge ME:
# MEmergeCutHeight <- c(0.2, 0.25, 0.3)
# maxBlockSize <- c(12000)
# ArgsbwModulesLLDF <- expand.grid(softPower = softPower
#                                  , minModSize = minModSize
#                                  , deepSplit = deepSplit
#                                  , MEmergeCutHeight = MEmergeCutHeight
#                                  , maxBlockSize = maxBlockSize)

pdf("../analysis/Allen_MEs_markers_corr.pdf", height=10, width=10)
for (i in 1:length(bwModulesLL)) {
  labeledHeatmap(markerMEcor[[i]]
                 , xLabels = colnames(markerMEcor[[i]])
                 , yLabels = rownames(markerMEcor[[i]])
                 , ySymbols = rownames(markerMEcor[[i]])
                 , colorLabels = FALSE
                 , colors = blueWhiteRed(50)
                 , setStdMargins = TRUE
                 , textMatrix = round(markerMEcor[[i]],2)
                 , cex.text = 0.9
                 , cex.lab = 0.9
                 , zlim = c(-1,1)
                 , main = paste(colnames(ArgsbwModulesLLDF)
                                , c(ArgsbwModulesLLDF[i,]), collapse=" "))
}
dev.off()
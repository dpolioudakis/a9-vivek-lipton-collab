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

pdf("../analysis/3b_MEs_markers_corr-2015-08-20.pdf", height=8, width=8)
for (i in 1:length(blockwiseMEs)) {
  labeledHeatmap(markerMEcor[[i]]
                 , xLabels = colnames(markerMEcor[[i]])
                 , yLabels = rownames(markerMEcor[[i]])
                 , ySymbols = rownames(markerMEcor[[i]])
                 , colorLabels = FALSE
                 , colors = greenWhiteRed(50)
                 , setStdMargins = TRUE
                 , textMatrix = signif(markerMEcor[[i]],2)
                 , cex.text = 1
                 , cex.lab = 1
                 , zlim = c(-1,1)
                 , main = paste(colnames(ArgsbwModulesLLDF)
                                , c(ArgsbwModulesLLDF[i,]), collapse=" "))
}
dev.off()
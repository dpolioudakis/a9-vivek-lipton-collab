# Cluster samples and construct modules

# Workflow
#   1a-allen-subset-to-basal-ganglia.R
#   1b-allen-combine-robes.R
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
expr.data <- t(array.data.subset.avg.probes.df)
colnames(expr.data) <- expr.data[1, ]
expr.data <- expr.data[-1, ]
# # Selecting 5000 genes with highest expression values (avg across samples)
# expr.data.top.5000 <- expr.data[,rank(-colMeans(expr.data))<=5000]
# expr.data= expr.data.top.5000

soft.power = 5
# Biweight midcorrelation is considered to be a good alternative to Pearson
# correlation since it is more robust to outliers.
adjacency = adjacency(expr.data, power= soft.power, corFnc= "bicor")

TOM = TOMsimilarity(adjacency)
diss.TOM = 1-TOM

save(adjacency, TOM, diss.TOM
     , file="../processed_data/allen_adjacency_TOM.rda")
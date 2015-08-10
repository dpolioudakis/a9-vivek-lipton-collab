# Cluster samples and construct modules

print("#######################################################################")
print("Starting allen-construct-network-modules.R script...")
sessionInfo()
 
library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/array.data.subset.avg.probes.rda")

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
gene.tree = hclust(as.dist(diss.TOM), method = "average")

# Module identification using hybrid tree cut:
# Function to construct modules, also merges modules based on ME
     # Args: minimum module size, cut height to merge ME, deep split
Make_modules <- function (min.mod.size, deep.split) {
     print("Treecut arguments:")
     # print(c("min.mod.size"=min.mod.size,"cut.height"=ME.merge.cut.height, "deep.split"=ds))
     print(c(min.mod.size, deep.split))
     tree = cutreeHybrid(dendro= gene.tree, pamRespectsDendro= FALSE
                         , minClusterSize= min.mod.size, cutHeight= 0.967
                         , deepSplit= deep.split, distM= as.matrix(diss.TOM))
     print("Table of genes per module:")
     print(table(tree$labels))
     tree$labels
}

# Merge modules based on ME function
Merge_modules_ME <- function (gene.module.color, ME.merge.cut.height) {
     # Call an automatic merging function
     # merged: The merged module colors
     # Cut height of 0.25, corresponds to a correlation of 0.75, to merge ME:
     merged <- mergeCloseModules(exprData= expr.data, colors= gene.module.color,
                                 cutHeight= ME.merge.cut.height)
     labels2colors(merged$colors)
}

# Test different parameters for constructing and merging modules
# Define arguments to test for cutreeHybrid
min.mod.sizes <- c(30,100,160)
deep.splits <- c(2,4)
ME.merge.cut.heights <- c(0.1,0.2,0.25)
modules.colors <- NULL
module.labels <- NULL
for (min.mod.size in min.mod.sizes) {
     for (deep.split in deep.splits) {
          # Test multiple cutreeHybrid parameters
          module <- Make_modules(min.mod.size, deep.split)
          for (ME.merge.cut.height in ME.merge.cut.heights) {
               # Test ME merge cut heights
               modules.colors <- cbind(modules.colors,
                                      Merge_modules_ME(module
                                                      , ME.merge.cut.height))
               # Make label from parameters used to make each module
               module.labels <- c(module.labels, paste(
                    "MMS=",min.mod.size
                    , " \nDS=",deep.split
                    , " \nMEcor=",ME.merge.cut.height
                    ))
          }
     }
}

sizeGrWindow(25,20)
pdf("../analysis/Dendro_test_module_parameters.pdf",height=25,width=20)
plotDendroAndColors(gene.tree
                    , modules.colors
                    , groupLabels=module.labels
                    , addGuide=TRUE
                    , dendroLabels=FALSE
                    , main="Dendrogram With Different Module Cutting Parameters")
dev.off()

save(expr.data, gene.tree, modules.colors, module.labels,
     file="../processed_data/allen_modules.rda")




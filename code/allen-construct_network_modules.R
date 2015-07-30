sessionInfo()
 
library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

print("Starting wgcna-allen.R script...")
################################################################################

# Choose soft-thresholding power

load("../processed_data/subset.sn.array.data.rda")

sn.expr.data <- t(subset.sn.array.data)
colnames(sn.expr.data) <- sn.expr.data[1, ]
sn.expr.data <- sn.expr.data[-1, ]
#2.b.1 Choosing the soft-thresholding power: analysis of network topology

pdf("../analysis/1.1_power_top_5000_expr.pdf", height=10, width=18)
# Choose a set of soft-thresholding powers
powers = c(1:30)

# Call the network topology analysis function
sft = pickSoftThreshold(
  sn.expr.data, powerVector= powers, verbose= 5, blockSize= 5000, corFnc= "bicor")

sn.expr.data.top.5000 <- sn.expr.data[,rank(-colMeans(sn.expr.data))<=5000]
sn.expr.data.random.5000 <- sn.expr.data[ , sample(ncol(sn.expr.data), 5000)]

sft.fxn <- function(expr.data) {
  sft = pickSoftThreshold(
    expr.data, powerVector= powers, verbose= 5, corFnc="bicor"
  )
}
sft.random.5000 <- sft.fxn(sn.expr.data.random.5000)
sft.top.5000 <- sft.fxn(sn.expr.data.top.5000)

sft.plots.fxn <- function(sft, output.file.path, plot.title) {
  # Plot the results:
  pdf(output.file.path, height=10, width=18)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
       , xlab="Soft Threshold (power)"
       , ylab="Scale Free Topology module Fit,signed R^2",type="n"
       , main = paste(plot.title))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  abline(h=0.80,col="blue")
  abline(h=0.70,col="orange")
  abline(h=0.60,col="green")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  dev.off()
}
sft.plots.fxn(sft.top.5000
              , "../analysis/1.1_power_top_5000_expr.pdf"
              , "Scale independence: Top 5000 most expressed probes")
sft.plots.fxn(sft.random.5000
              , "../analysis/1.1_power_random_5000.pdf"
              , "Scale independence: 5000 random probes")
sft.plots.fxn(sft, "../analysis/1.1_power.pdf", "Scale independence")
################################################################################

# Cluster samples and construct modules

load("../processed_data/subset.sn.array.data.rda")

# Transpose expression data and assign column names from column 1
sn.expr.data <- t(subset.sn.array.data)
colnames(sn.expr.data) <- sn.expr.data[1, ]
sn.expr.data <- sn.expr.data[-1, ]
# Selecting 5000 genes with highest expression values (avg across samples)
sn.expr.data.top.5000 <- sn.expr.data[,rank(-colMeans(sn.expr.data))<=5000]
expr.data= sn.expr.data.top.5000

soft.power = 4
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
     file="../processed_data/allen.sn.modules.rda")
################################################################################




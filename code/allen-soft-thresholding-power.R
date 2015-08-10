# Make graphs to choose soft-thresholding power

print("#######################################################################")
print("Starting allen-soft-thresholding-power.R script...")
sessionInfo()

library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/array_data_subset.avg.probes.rda")

expr.data <- t(array.data.subset.avg.probes.df)
colnames(expr.data) <- expr.data[1, ]
expr.data <- expr.data[-1, ]
#2.b.1 Choosing the soft-thresholding power: analysis of network topology

pdf("../analysis/1.1_power_top_5000_expr.pdf", height=10, width=18)
# Choose a set of soft-thresholding powers
powers = c(1:30)

# Call the network topology analysis function
sft = pickSoftThreshold(
  expr.data, powerVector= powers, verbose= 5, blockSize= 5000, corFnc= "bicor")

expr.data.top.5000 <- expr.data[,rank(-colMeans(expr.data))<=5000]
expr.data.random.5000 <- expr.data[ , sample(ncol(expr.data), 5000)]

sft.fxn <- function(expr.data) {
  sft = pickSoftThreshold(
    expr.data, powerVector= powers, verbose= 5, corFnc="bicor"
  )
}
sft.random.5000 <- sft.fxn(expr.data.random.5000)
sft.top.5000 <- sft.fxn(expr.data.top.5000)

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

print("Ending allen-soft-thresholding-power.R script...")
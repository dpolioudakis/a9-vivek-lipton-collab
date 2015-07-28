sessionInfo()
 
library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

print("Starting wgcna-allen.R script...")
################################################################################

# Cluster samples and compare to meta data

load("../processed_data/subset.sn.array.data.rda")

# Transpose expression data and assign column names from column 1
sn.expr.data <- t(subset.sn.array.data)
colnames(sn.expr.data) <- sn.expr.data[1, ]
sn.expr.data <- sn.expr.data[-1, ]

# Selecting 5000 genes with highest expression values (avg across samples)
sn.expr.data.top.5000 <- sn.expr.data[,rank(-colMeans(sn.expr.data))<=5000]


expr.data= sn.expr.data.top.5000

# Allen brain IDs derived from Allan brain folder names
brain.ids <- list(
     "178236545",
     "178238266",
     "178238316",
     "178238359",
     "178238373",
     "178238387"
)

# Add column of brain IDs to meta data for each brain
datTraits <- mapply(function(brain.meta.data, brain.id) cbind(brain.meta.data, brain.id)
                    , subset.sn.meta.data, brain.ids, SIMPLIFY= FALSE)
# Combine list of meta data for each brain into one data frame
datTraits <- do.call(rbind, datTraits)

# Cluster samples
sampleTree2 = hclust(dist(expr.data), method = "average")

# Convert traits to a color representation: for numeric traits white means low, 
# red means high, grey means missing entry
traitColors= data.frame(labels2colors(datTraits[,1:7])
                        , numbers2colors(datTraits[,8:13])
                        , labels2colors(datTraits[,14]))
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# Select specific traits
# Convert traits to a color representation: for numeric traits white means low, 
# red means high, grey means missing entry
traitColors= data.frame(labels2colors(datTraits[,c(2,5,14)])
                        , numbers2colors(datTraits[,c(3,8:10)]))
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits)[c(2,5,14,3,8:10)],
                    main = "Sample dendrogram and trait heatmap")
################################################################################

# Choose soft-thresholding power

load("../processed_data/subset.sn.array.data.rda")

sn.expr.data <- t(subset.sn.array.data)
colnames(sn.expr.data) <- sn.expr.data[1, ]
sn.expr.data <- sn.expr.data[-1, ]
#2.b.1 Choosing the soft-thresholding power: analysis of network topology

pdf("1.1_power_top_5000_expr.pdf", height=10, width=18)
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

sft.plots.fxn <- function(sft, output.file.name, plot.title) {
  # Plot the results:
  pdf(output.file.name, height=10, width=18)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
       , xlab="Soft Threshold (power)"
       , ylab="Scale Free Topology Model Fit,signed R^2",type="n"
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
sft.plots.fxn(sft.top.5000, "1.1_power_top_5000_expr.pdf", "Scale independence: Top 5000 most expressed probes")
sft.plots.fxn(sft.random.5000, "1.1_power_random_5000.pdf", "Scale independence: 5000 random probes")
sft.plots.fxn(sft, "1.1_power.pdf", "Scale independence")
################################################################################

expr.data= sn.expr.data.top.5000

soft.power = 12;
# Biweight midcorrelation is considered to be a good alternative to Pearson
# correlation since it is more robust to outliers.
adjacency = adjacency(expr.data, power= soft.power, corFnc= "bicor")

TOM = TOMsimilarity(adjacency)
diss.TOM = 1-TOM
gene.tree = hclust(as.dist(diss.TOM), method = "average")
sizeGrWindow(12,9)
plot(gene.tree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
min.module.size = 30
# Module identification using dynamic tree cut:
dynamic.mods = cutreeDynamic(dendro = gene.tree, distM= diss.TOM
                             , deepSplit= 2, pamRespectsDendro= FALSE
                             , minClusterSize= min.module.size)
table(dynamic.mods)

# Convert numeric lables into colors
dynamic.colors = labels2colors(dynamic.mods)
table(dynamic.colors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(gene.tree, dynamic.colors, "Dynamic Tree Cut"
                    , dendroLabels= FALSE, hang= 0.03
                    , addGuide= TRUE, guideHang= 0.05
                    , main= "Gene dendrogram and module colors")

# Calculate eigengenes
ME.list = moduleEigengenes(expr.data, colors= dynamic.colors)
MEs = ME.list$eigengenes
# Calculate dissimilarity of module eigengenes
ME.diss = 1-cor(MEs);
# Cluster module eigengenes
ME.tree = hclust(as.dist(ME.diss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(ME.tree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Cut height of 0.25, corresponding to a correlation of 0.75, to merge ME:
ME.diss.thres = 0.25
# Plot the cut line into the dendrogram
abline(h= ME.diss.thres, col= "red")
# Call an auTOMatic merging function
merge = mergeCloseModules(expr.data, dynamic.colors, cutHeight= ME.diss.thres, verbose = 3)
# The merged module colors
merged.colors = merge$colors
# Eigengenes of the new merged modules:
merged.MEs = merge$newMEs

# Plot gene dendogram again with original and merged module colors underneath
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(gene.tree, cbind(dynamic.colors, merged.colors)
                    , c("Dynamic Tree Cut", "Merged dynamic")
                    , dendroLabels= FALSE, hang= 0.03
                    , addGuide= TRUE, guideHang= 0.05)

####

# Cluster samples and compare to meta data

load("../processed_data/subset.sn.array.data.rda")

# Transpose expression data and assign column names from column 1
sn.expr.data <- t(subset.sn.array.data)
colnames(sn.expr.data) <- sn.expr.data[1, ]
sn.expr.data <- sn.expr.data[-1, ]

# Selecting 5000 genes with highest expression values (avg across samples)
sn.expr.data.top.5000 <- sn.expr.data[,rank(-colMeans(sn.expr.data))<=5000]


expr.data= sn.expr.data.top.5000

# Allen brain IDs derived from Allan brain folder names
brain.ids <- list(
     "178236545",
     "178238266",
     "178238316",
     "178238359",
     "178238373",
     "178238387"
)

# Add column of brain IDs to meta data for each brain
datTraits <- mapply(function(brain.meta.data, brain.id) cbind(brain.meta.data, brain.id)
                    , subset.sn.meta.data, brain.ids, SIMPLIFY= FALSE)
# Combine list of meta data for each brain into one data frame
datTraits <- do.call(rbind, datTraits)

traitmat <- datTraits[,c(2,5,14,3)]
geneSigs=matrix(NA,nrow=4,ncol=ncol(expr.data)) # create a vector to hold the data


for(i in 1:ncol(geneSigs)) {
     
     exprvec= as.numeric(expr.data[,i]) # get the expression vector for ith gene
     slab.num= sqrt(max(summary(lm(exprvec~as.factor(traitmat[,1])))$adj.r.squared,0))
     structure.acroynm= sqrt(max(summary(lm(exprvec~as.factor(traitmat[,2])))$adj.r.squared,0))
     brain.id= sqrt(max(summary(lm(exprvec~as.factor(traitmat[,3])))$adj.r.squared,0))
     well.id=bicor(traitmat[,4],exprvec)# calculate r correlation value for numeric variables
     # Well IDs are numbers, ie: "160535191 160535175 160091869 160091634"

     geneSigs[, i]=c(slab.num, structure.acroynm, brain.id, well.id)
}

geneSigs[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100),lim=c(0,1)) # For categorical
geneSigs[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=FALSE,centered=FALSE,blueWhiteRed(100),lim=c(0,1)) # For categorical
geneSigs[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=FALSE,centered=FALSE,blueWhiteRed(100),lim=c(0,1)) # For categorical
geneSigs[4,] =numbers2colors(as.numeric(geneSigs[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
rownames(geneSigs)=c("slab.num","structure.acroynm","brain.id", "well.id")

####

mColorh <- mLabelh <- colorLabels <- NULL  
for (minModSize in c(40,100,160)) {
     for (dthresh in c(0.1,0.2,0.25)) {
          for (ds in c(2,4)) {
               print("Trying parameters:")
               print(c(minModSize,dthresh,ds))
               tree = cutreeHybrid(dendro = gene.tree, pamStage=FALSE,
                                   minClusterSize = minModSize, cutHeight = 0.9999,
                                   deepSplit = ds, distM = as.matrix(dissTOM))
               
               merged <- mergeCloseModules(exprData = expr.data,colors = tree$labels,
                                           cutHeight = dthresh)
               mColorh <- cbind(mColorh,labels2colors(merged$colors))
               mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
          }
     }
}

mColorh1=cbind(mColorh,geneSigs[1,],geneSigs[2,],geneSigs[3,],geneSigs[4,])
mLabelh1=c(mLabelh,rownames(geneSigs))
plotDendroAndColors(gene.tree,mColorh,groupLabels=mLabelh,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different
Module Cutting Parameters")

pdf("Signed_New_Dendro1.pdf",height=25,width=20)
plotDendroAndColors(gene.tree,mColorh1,groupLabels=mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different Module Cutting Parameters")
dev.off()






--------------------------


###################### WGCNA TOM##############

softPower = 20;
#Biweight midcorrelation is considered to be a good alternative to Pearson correlation since it is more robust to outliers.
adjacency = adjacency(expr.data, power = softPower, type = "signed",corFnc="bicor");

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

gene.tree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
save(TOM,dissTOM,gene.tree,adjacency,softPower,file="TOM_Ctx.rda")



### Relating dendrogram with traits
datTraits<- targets[,c(2,6,7)]

traitmat=as.data.frame(cbind(as.numeric(factor(datTraits[,1])),as.numeric(datTraits[,2]),as.numeric(datTraits[,3])))# convert categorical variables in factor and numeric as numeric

marker=c("ENSG00000180176","ENSG00000157542","ENSG00000165646","ENSG00000157388","ENSG00000162761","ENSG00000125798","ENSG00000153234","ENSG00000165092") ##TYH, KCNJ7,VMAT2,CACNA1D,LMX1A,FOXA2,NURR1,ALDH1
anti.marker=c("ENSG00000104327","ENSG00000172137","ENSG00000165588","ENSG00000165655","ENSG00000148680") ##CALB1/2, OTX2, NOLZ1, HTR7


rownames(traitmat)=rownames(datTraits)
colnames(traitmat)=c("Cell-Batch","RIN","Ratio260.280")

geneSigs=matrix(NA,nrow=3,ncol=ncol(expr.data)) # create a vector to hold the data

for(i in 1:ncol(geneSigs)) {
     
     exprvec=as.numeric(expr.data[,i]) # get the expression vector for ith gene
     batchr=sqrt(max(summary(lm(exprvec~as.factor(traitmat[,1])))$adj.r.squared,0))
     rinr=bicor(traitmat[,2],exprvec)# calculate r correlation value for numeric variables
     ratior=cor(traitmat[,3],exprvec)
     
     geneSigs[,i]=c(batchr,rinr,ratior)
     
     cat('Done for gene...',i,'\n')
}

geneSigs[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1)) # For categorical
geneSigs[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1)) 

geneSigs[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1)) 


rownames(geneSigs)=c("Cell-type.Batch","RIN","Ratio260.280")



mColorh <- mLabelh <- colorLabels <- NULL  
for (minModSize in c(40,100,160)) {
     for (dthresh in c(0.1,0.2,0.25)) {
          for (ds in c(2,4)) {
               print("Trying parameters:")
               print(c(minModSize,dthresh,ds))
               tree = cutreeHybrid(dendro = gene.tree, pamStage=FALSE,
                                   minClusterSize = minModSize, cutHeight = 0.9999,
                                   deepSplit = ds, distM = as.matrix(dissTOM))
               
               merged <- mergeCloseModules(exprData = expr.data,colors = tree$labels,
                                           cutHeight = dthresh)
               mColorh <- cbind(mColorh,labels2colors(merged$colors))
               mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
          }
     }
}

mColorh1=cbind(mColorh,geneSigs[1,],geneSigs[2,],geneSigs[3,])
mLabelh1=c(mLabelh,rownames(geneSigs))

pdf("Signed_New_Dendro1.pdf",height=25,width=20)
plotDendroAndColors(gene.tree,mColorh1,groupLabels=mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different Module Cutting Parameters")
dev.off()



print("End of wgcna-allen.R script...")



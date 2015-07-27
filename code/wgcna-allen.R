sessionInfo()
 
library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

print("Starting wgcna-allen.R script...")
################################################################################

#Load Allen Brain microarray data and corresponding metadata
print("Loading Allen Brain microarray data...")

#Functions to load data
load.array.data <- function (bio.rep.folder) {
     array.data.path <- paste(
          "../raw_data/allen_data/", bio.rep.folder, "/MicroarrayExpression.csv"
          , sep= ""
          )
     read.csv(array.data.path, header=FALSE)
}
load.meta.data <- function (bio.rep.folder) {
     meta.data.path <- paste(
          "../raw_data/allen_data/", bio.rep.folder, "/SampleAnnot.csv"
          , sep= ""
          )
     read.csv(meta.data.path)
}

#Sub folders with Allen Brain span data, each folder is 1 brain
array.data.folder.names <- list(
     "178236545-2015-07-15",
     "178238266-2015-07-15",
     "178238316-2015-07-15",
     "178238359-2015-07-15",
     "178238373-2015-07-15",
     "178238387-2015-07-15"
)

#Load the Allen Brain microarray data
#Columns correspond to samples and samples are listed in the same order in the
#meta and array data
exp.array.data <- lapply(array.data.folder.names, load.array.data)
meta.array.data <- lapply(array.data.folder.names, load.meta.data)
               
# save(array.data, file="../processed.data/all.allen.array.data.rda")
################################################################################

#Subset the data to retain only the substantia nigra samples
print("Subsetting data to select the substantia nigra samples...")

# Function to subset the meta data
subset.sn.meta.data.fxn <- function (bio.rep.meta.data) {
     subset(bio.rep.meta.data, grepl("substantia nigra", structure_name))
}

#Function to subset the array expression data
subset.array.data.fxn <- function (bio.rep.exp.array.data, bio.rep.meta.data) {
     subset.array.data <- bio.rep.exp.array.data[ , grepl(
          "substantia nigra", bio.rep.meta.data$structure_name)]
}

#Subset the meta data
subset.sn.meta.data <- lapply(meta.array.data, subset.sn.meta.data.fxn)
#Meta data substantia nigra number of samples and number of features
print("Meta data substantia nigra number of samples and number of features:")
sapply(subset.sn.meta.data, dim)
#Total number of substantia nigra samples
print("Total number of substantia nigra samples:")
sum(sapply(subset.sn.meta.data, dim)[1, ])

#Subset the array expression data
subset.sn.array.data <- do.call(cbind,
     mapply(subset.array.data.fxn, exp.array.data, meta.array.data))
#Add the probe ID numbers as column 1
subset.sn.array.data <- cbind(probe=exp.array.data[[1]][,1], subset.sn.array.data)
print("Dimensions of subset expression data, columns are samples:")
print("(row 1 is Probe ID)")
dim(subset.sn.array.data)

save(subset.sn.array.data, subset.sn.meta.data,
     file="../processed_data/subset.sn.array.data.rda")
################################################################################

load("../processed_data/subset.sn.array.data.rda")

sn.expr.data <- t(subset.sn.array.data)
colnames(sn.expr.data) <- sn.expr.data[1, ]
sn.expr.data <- sn.expr.data[-1, ]

sn.expr.data.top.5000 <- sn.expr.data[,rank(-colMeans(sn.expr.data))<=5000]

datExpr= sn.expr.data.top.5000

brain.ids <- list(
     "178236545",
     "178238266",
     "178238316",
     "178238359",
     "178238373",
     "178238387"
)

datTraits <- mapply(function(brain.meta.data, brain.id) cbind(brain.meta.data, brain.id)
                    , subset.sn.meta.data, brain.ids, SIMPLIFY= FALSE)
datTraits <- do.call(rbind, datTraits)

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry

traitColors= data.frame(labels2colors(datTraits[,1:7])
                        , numbers2colors(datTraits[,8:13])
                        , labels2colors(datTraits[,14]))
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


traitColors= data.frame(labels2colors(datTraits[,c(2,5,14)])
                        , numbers2colors(datTraits[,c(3,8:10)]))

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits)[c(2,5,14,3,8:10)],
                    main = "Sample dendrogram and trait heatmap")
################################################################################

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
datExpr= sn.expr.data.top.5000

softPower = 12;
#Biweight midcorrelation is considered to be a good alternative to Pearson correlation since it is more robust to outliers.
adjacency = adjacency(datExpr, power= softPower, corFnc= "bicor")

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM




--------------------------
     
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits)[2,5,14,3,8:10],
                    main = "Sample dendrogram and trait heatmap")





###################### WGCNA TOM##############

softPower = 20;
#Biweight midcorrelation is considered to be a good alternative to Pearson correlation since it is more robust to outliers.
adjacency = adjacency(datExpr, power = softPower, type = "signed",corFnc="bicor");

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
save(TOM,dissTOM,geneTree,adjacency,softPower,file="TOM_Ctx.rda")



### Relating dendrogram with traits
datTraits<- targets[,c(2,6,7)]

traitmat=as.data.frame(cbind(as.numeric(factor(datTraits[,1])),as.numeric(datTraits[,2]),as.numeric(datTraits[,3])))# convert categorical variables in factor and numeric as numeric

marker=c("ENSG00000180176","ENSG00000157542","ENSG00000165646","ENSG00000157388","ENSG00000162761","ENSG00000125798","ENSG00000153234","ENSG00000165092") ##TYH, KCNJ7,VMAT2,CACNA1D,LMX1A,FOXA2,NURR1,ALDH1
anti.marker=c("ENSG00000104327","ENSG00000172137","ENSG00000165588","ENSG00000165655","ENSG00000148680") ##CALB1/2, OTX2, NOLZ1, HTR7


rownames(traitmat)=rownames(datTraits)
colnames(traitmat)=c("Cell-Batch","RIN","Ratio260.280")

geneSigs=matrix(NA,nrow=3,ncol=ncol(datExpr)) # create a vector to hold the data

for(i in 1:ncol(geneSigs)) {
     
     exprvec=as.numeric(datExpr[,i]) # get the expression vector for ith gene
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
               tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
                                   minClusterSize = minModSize, cutHeight = 0.9999,
                                   deepSplit = ds, distM = as.matrix(dissTOM))
               
               merged <- mergeCloseModules(exprData = datExpr,colors = tree$labels,
                                           cutHeight = dthresh)
               mColorh <- cbind(mColorh,labels2colors(merged$colors))
               mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
          }
     }
}

mColorh1=cbind(mColorh,geneSigs[1,],geneSigs[2,],geneSigs[3,])
mLabelh1=c(mLabelh,rownames(geneSigs))

pdf("Signed_New_Dendro1.pdf",height=25,width=20)
plotDendroAndColors(geneTree,mColorh1,groupLabels=mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different Module Cutting Parameters")
dev.off()



print("End of wgcna-allen.R script...")



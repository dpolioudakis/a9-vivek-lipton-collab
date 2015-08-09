sessionInfo()

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

################################################################################

# Compare modules to meta data

load("../processed_data/subset.sn.array.data.rda")
load("../processed_data/allen.sn.modules.rda")

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
trait.data <- mapply(function(brain.meta.data, brain.id)
     cbind(brain.meta.data, brain.id)
     , subset.sn.meta.data, brain.ids, SIMPLIFY= FALSE)
# Combine list of meta data for each brain into one data frame
trait.data <- do.call(rbind, trait.data)

# Select specific traits
traitmat <- trait.data[,c(2,5,14,3)]

gene.sigs=matrix(NA,nrow=4,ncol=ncol(expr.data)) # create a vector to hold the data
for(i in 1:ncol(gene.sigs)) {
     
     exprvec= as.numeric(expr.data[,i]) # get the expression vector for ith gene
     slab.num= sqrt(max(summary(lm(exprvec~as.factor(traitmat[,1])))$adj.r.squared,0))
     structure.acroynm= sqrt(max(summary(lm(exprvec~as.factor(traitmat[,2])))$adj.r.squared,0))
     brain.id= sqrt(max(summary(lm(exprvec~as.factor(traitmat[,3])))$adj.r.squared,0))
     well.id=bicor(traitmat[,4],exprvec)# calculate r correlation value for numeric variables
     # Well IDs are numbers, ie: "160535191 160535175 160091869 160091634"
     
     gene.sigs[, i]=c(slab.num, structure.acroynm, brain.id, well.id)
}


gene.sigs[1,] =numbers2colors(as.numeric(gene.sigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100),lim=c(0,1)) # For categorical
gene.sigs[2,] =numbers2colors(as.numeric(gene.sigs[2,]),signed=FALSE,centered=FALSE,blueWhiteRed(100),lim=c(0,1)) # For categorical
gene.sigs[3,] =numbers2colors(as.numeric(gene.sigs[3,]),signed=FALSE,centered=FALSE,blueWhiteRed(100),lim=c(0,1)) # For categorical
gene.sigs[4,] =numbers2colors(as.numeric(gene.sigs[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
rownames(gene.sigs)=c("slab.num","structure.acroynm","brain.id", "well.id")

# Parameters used to construct module [,3]: DS= 2, MMS= 30, MEcor= 0.25
modules.traits=data.frame(modules.colors[,3]
                          , gene.sigs[1,]
                          , gene.sigs[2,]
                          , gene.sigs[3,]
                          , gene.sigs[4,])
modules.traits.labels=c(module.labels[3],rownames(gene.sigs))

pdf("../analysis/dendro_modules_traits_corr.pdf",height=25,width=20)
plotDendroAndColors(gene.tree,modules.traits,groupLabels=modules.traits.labels,addGuide=TRUE,dendroLabels=FALSE,main="Dendro and traits correlation")
dev.off()

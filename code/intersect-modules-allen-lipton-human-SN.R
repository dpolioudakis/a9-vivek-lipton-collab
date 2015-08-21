# Intersect Vivek's Lipton substantia nigra sample A9 modules and Allen brain A9
# modules

library(WGCNA)
library(reshape2)

load("../Vivek_WGCNA_Lipton_A9_SN/Vivek_SN_mod_colors.rda")
load("../processed_data/allen_BW_modules.rda")

# Split by model into lists or dataframes of genes in each model
allenGenesModels <- split(colnames(exprData), bwModulesLL[[3]]$colors)
humanSNgenesModels <- split(geneInfo[,1:2], geneInfo$`Initially Assigned Module Color`)
# Dataframe of all combinations of models from Allen and Lipton human Substantia
# Nigra
modelCombos <- expand.grid(names(allenGenesModels), names(humanSNgenesModels))

# Intersection genes for every combination of models
allenHumanSNintersect <- apply(modelCombos, 1
  , function(x) length(intersect(allenGenesModels[[x[1]]]
                                 , humanSNgenesModels[[x[2]]]$GeneSymbol)))
# Add combination of models labels
allenHumanSNintersect <- cbind(modelCombos, allenHumanSNintersect)
# Reshape into matrix
allenHumanSNintersect <- dcast(allenHumanSNintersect
                               , Var1~Var2, value.var="allenHumanSNintersect")

# Filter for intersection of genes in Allen and Lipton human substantia nigra
totalIntersectGenesProfiled <- intersect(geneInfo$GeneSymbol
                                         , colnames(exprData))
# Filter each model
allenGenesModels <- lapply(allenGenesModels, function(x)
                             intersect(x, totalIntersectGenesProfiled))
humanSNgenesModels <- lapply(humanSNgenesModels, function(x)
  intersect(x$GeneSymbol, totalIntersectGenesProfiled))

# Intersection genes for every combination of models
allenHumanSNintersect <- apply(modelCombos, 1, function(x)
  length(intersect(allenGenesModels[[x[1]]], humanSNgenesModels[[x[2]]])))
allenHumanSNintersect <- cbind(modelCombos, allenHumanSNintersect)
allenHumanSNintersect <- dcast(allenHumanSNintersect
                               , Var1~Var2, value.var="allenHumanSNintersect")

phyper(q=3, m=12, n=7147-12, k=211, lower.tail=FALSE, log.p=FALSE)
length(allenGenesModels[[29]])
length(humanSNgenesModels$skyblue)
length(totalIntersectGenesProfiled)

sapply(allenHumanSNintersect[29,], function(x)
         phyper(
           q=(x-1) #size of overlap - 1
         , m=length(allenGenesModels[[29]]) # genes in Allen module
         # total number of genes shared between Lipton and Allen - genes in
         # Lipton human SN model
         , n=length(totalIntersectGenesProfiled)-length(humanSNgenesModels$skyblue)
         , k=length(humanSNgenesModels$skyblue)  # genes in Lipton human SN model
         , lower.tail=FALSE, log.p=FALSE))

pvals <- apply(allenHumanSNintersect, 1, function(x)
  #print(length(allenGenesModels[[x[1]]])))
  phyper(
    q=as.numeric(x[3])-1 #size of overlap - 1
    , m=length(allenGenesModels[[x[1]]]) # genes in Allen module
    # total number of genes shared between Lipton and Allen - genes in
    # Lipton human SN model
    , n=length(totalIntersectGenesProfiled)-length(humanSNgenesModels[[x[2]]])
    , k=length(humanSNgenesModels[[x[2]]])  # genes in Lipton human SN model
    , lower.tail=FALSE, log.p=FALSE))
pvals <- sapply(pvals, function(x) round(x,2))

allenHumanSNintersectPvals <- cbind(modelCombos, pvals)
allenHumanSNintersectPvals <- dcast(allenHumanSNintersectPvals
                               , Var1~Var2, value.var="pvals")
allenHumanSNintersectPvals <- data.matrix(allenHumanSNintersectPvals)
rownames(allenHumanSNintersectPvals) <- allenHumanSNintersectPvals[ ,1]
allenHumanSNintersectPvals <- allenHumanSNintersectPvals[ ,-1]


allenHumanSNintersect[29,]
allenHumanSNintersect$skyblue

allenHumanSNintersect
data.class(allenHumanSNintersect)
allenHumanSNintersect <- data.matrix(allenHumanSNintersect)
data.class(allenHumanSNintersect[1,1])

labeledHeatmap(allenHumanSNintersect, xLabels=NULL)
labeledHeatmap(allenHumanSNintersectPvals, xLabels=NULL)

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
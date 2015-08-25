# Intersect Vivek's Lipton substantia nigra sample A9 modules and Allen brain A9
# modules

library(WGCNA)
library(reshape2)
library(biomaRt)

load("../Vivek_WGCNA_Lipton_A9_SN/Vivek_SN_mod_colors.rda")
load("../processed_data/allen_BW_modules.rda")

# Variable for minModSize parameter used in blockwiseModules WGCNA function
# to record in output graph titles
minModSize <- "30"

# bwModulesLL is list of modules from different blockwiseModules parameters used
# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 11

# Split by model into lists or dataframes of genes in each model
allenGenesModels <- split(colnames(exprData)
                          , bwModulesLL[[modNetworkToUse]]$colors)
humanSNgenesModelsLL <- split(geneInfo[,1:2]
                              , geneInfo$`Initially Assigned Module Color`)
# Dataframe of all combinations of models from Allen and Lipton human Substantia
# Nigra
moduleCombos <- expand.grid(names(allenGenesModels), names(humanSNgenesModelsLL))

# Intersection genes for every combination of models
allenHumanSNintersect <- apply(moduleCombos, 1
  , function(x) length(intersect(allenGenesModels[[x[1]]]
                                 , humanSNgenesModelsLL[[x[2]]])))
# Add combination of models labels
allenHumanSNintersect <- cbind(moduleCombos, allenHumanSNintersect)
# Reshape into matrix
allenHumanSNintersect <- dcast(allenHumanSNintersect
                               , Var1~Var2, value.var="allenHumanSNintersect")

# Filter for intersection of all genes in Allen and Lipton human substantia nigra
totalIntersectGenesProfiled <- intersect(geneInfo$GeneSymbol
                                         , colnames(exprData))
# Filter each model in Allen and Lipton for all genes profiled that intersect
allenGenesModels <- lapply(allenGenesModels, function(x)
                             intersect(x, totalIntersectGenesProfiled))
humanSNgenesModelsLL <- lapply(humanSNgenesModelsLL, function(x)
  intersect(x$GeneSymbol, totalIntersectGenesProfiled))

# Intersection genes for every combination of models
allenHumanSNintersect <- apply(moduleCombos, 1, function(x)
  length(intersect(allenGenesModels[[x[1]]], humanSNgenesModelsLL[[x[2]]])))
allenHumanSNintersect <- cbind(moduleCombos, allenHumanSNintersect)
allenHumanSNintersect <- dcast(allenHumanSNintersect
                               , Var1~Var2, value.var="allenHumanSNintersect")

# Move row names from column 1 to row names
rownames(allenHumanSNintersect) <- allenHumanSNintersect[ ,1]
allenHumanSNintersect <- allenHumanSNintersect[ ,-1]
# Initialize empty matrix to hold pvalues
# Set dimension labels as module and number of genes in module
allenModels <- split(colnames(exprData)
                          , bwModulesLL[[modNetworkToUse]]$colors)
allenLabels <- paste(names(allenGenesModels), sapply(allenModels, length))
liptonLabels <- paste(names(humanSNgenesModelsLL), sapply(humanSNgenesModelsLL, length))
pVals <- matrix( , nrow = nrow(allenHumanSNintersect)
                 , ncol = ncol(allenHumanSNintersect)
                 , dimnames = list(allenLabels, liptonLabels))
paste(names(allenGenesModels), sapply(allenGenesModels, length))

# Calc pbinomal of interesection for each model in Lipton human substantia nigra
# with each model in Allen
for (j in 1:(ncol(allenHumanSNintersect))) {
  for (i in 1:nrow(allenHumanSNintersect)) {
#     pVals[i,j] <- phyper(
#         q = (as.numeric(allenHumanSNintersect[i,j])-1)
#       , m = length(allenGenesModels[[i]]) # genes in Allen module
#       , n = length(totalIntersectGenesProfiled)-length(humanSNgenesModelsLL[[j]])
#       , k = length(humanSNgenesModelsLL[[j]])
#       , lower.tail = FALSE, log.p = FALSE)
    pVals[i,j] <- phyper(
      q = (as.numeric(allenHumanSNintersect[i,j])-1)
      , m = length(allenGenesModels[[i]]) # genes in Allen module
      , n = length(totalIntersectGenesProfiled)-length(allenGenesModels[[i]])
      , k = length(humanSNgenesModelsLL[[j]])
      , lower.tail = FALSE, log.p = FALSE)
  }
}
# Format number of genes that intersect with respective pbinom
textMatrix <- paste(as.matrix(allenHumanSNintersect), "\n", round(pVals,3), sep = "")
# Plot pbinoms with number of genes that intersect included on label
pdf(paste("../analysis/intersection_Allen_human_SN-minModSize", minModSize, ".pdf"
          , sep = ""), height=12, width=12)
labeledHeatmap(pVals
               , yLabels = rownames(pVals)
               , xLabels = colnames(pVals)
               , textMatrix = textMatrix
               , setStdMargins = TRUE
               , cex.text = 0.5
               , cex.lab = 1
               , zlim = c(0,1)
               , main = "pbinom Lipton Human Substantia nigra and Allen"
               )
dev.off()

# List of genes that is intersection of each Allen Lipton Human SN module
# combination.
# List is in order of module combinations listed in moduleCombos
intModsGenesLL <- apply(moduleCombos, 1, function(x)
  intersect(allenGenesModels[[x[1]]], humanSNgenesModelsLL[[x[2]]]))

save(intModsGenesLL, moduleCombos
   , file = paste("../processed_data/allen_human_SN_modules_intersect_modsize"
                  , minModSize, ".rda", sep = ""))
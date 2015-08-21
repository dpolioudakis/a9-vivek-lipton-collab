# Intersect Vivek's Lipton substantia nigra sample A9 modules and Allen brain A9
# modules

library(WGCNA)
library(reshape2)

load("../Vivek_WGCNA_Lipton_A9_SN/Vivek_SN_mod_colors.rda")
load("../processed_data/allen_BW_modules.rda")

# Split by model into lists or dataframes of genes in each model
allenGenesModels <- split(colnames(exprData), bwModulesLL[[18]]$colors)
humanSNgenesModelsLL <- split(geneInfo[,1:2]
                              , geneInfo$`Initially Assigned Module Color`)
# Dataframe of all combinations of models from Allen and Lipton human Substantia
# Nigra
modelCombos <- expand.grid(names(allenGenesModels), names(humanSNgenesModelsLL))

# Intersection genes for every combination of models
allenHumanSNintersect <- apply(modelCombos, 1
  , function(x) length(intersect(allenGenesModels[[x[1]]]
                                 , humanSNgenesModelsLL[[x[2]]])))
# Add combination of models labels
allenHumanSNintersect <- cbind(modelCombos, allenHumanSNintersect)
# Reshape into matrix
allenHumanSNintersect <- dcast(allenHumanSNintersect
                               , Var1~Var2, value.var="allenHumanSNintersect")

# Filter for intersection of all genes in Allen and Lipton human substantia nigra
totalIntersectGenesProfiled <- intersect(geneInfo$GeneSymbol
                                         , colnames(exprData))
# Filter each model in Allen and Lipton for all genes that intersect
allenGenesModels <- lapply(allenGenesModels, function(x)
                             intersect(x, totalIntersectGenesProfiled))
humanSNgenesModelsLL <- lapply(humanSNgenesModelsLL, function(x)
  intersect(x$GeneSymbol, totalIntersectGenesProfiled))

# Intersection genes for every combination of models
allenHumanSNintersect <- apply(modelCombos, 1, function(x)
  length(intersect(allenGenesModels[[x[1]]], humanSNgenesModelsLL[[x[2]]])))
allenHumanSNintersect <- cbind(modelCombos, allenHumanSNintersect)
allenHumanSNintersect <- dcast(allenHumanSNintersect
                               , Var1~Var2, value.var="allenHumanSNintersect")

# Move row names from column 1 to row names
rownames(allenHumanSNintersect) <- allenHumanSNintersect[ ,1]
allenHumanSNintersect <- allenHumanSNintersect[ ,-1]
# Initialize empty matrix to hold pvalues
pVals <- matrix( , nrow = nrow(allenHumanSNintersect)
                 , ncol = ncol(allenHumanSNintersect)
                 , dimnames = list(rownames(allenHumanSNintersect)
                                   , colnames(allenHumanSNintersect)))
# Calc pbinomal of interesection for each model in Lipton human substantia nigra
# with each model in Allen
for (j in 1:(ncol(allenHumanSNintersect))) {
  for (i in 1:nrow(allenHumanSNintersect)) {
    pVals[i,j] <- phyper(
        q = (as.numeric(allenHumanSNintersect[i,j])-1)
      , m = length(allenGenesModels[[i]]) # genes in Allen module
      , n = length(totalIntersectGenesProfiled)-length(humanSNgenesModelsLL[[j]])
      , k = length(humanSNgenesModelsLL[[j]])
      , lower.tail = FALSE, log.p = FALSE)
  }
}
# Format number of genes that intersect with respective pbinom
textMatrix <- paste(as.matrix(allenHumanSNintersect), "\n", round(pVals,3), sep = "")
# Plot pbinoms with number of genes that intersect included on label
pdf(paste("../analysis/intersection_Allen_human_SN-minModSize100-", Sys.Date(), ".pdf"
          , sep = ""), height=8, width=8)
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
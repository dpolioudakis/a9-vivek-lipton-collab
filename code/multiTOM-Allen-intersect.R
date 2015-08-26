# Intersect Vivek's Lipton substantia nigra sample A9 modules and Allen brain A9
# modules

library(reshape2)
library(biomaRt)

load("../processed_data/allen_BW_modules.rda")

multiTOMdataDF <- read.csv(
  "../processed_data/multiTOM_allen_neighbors1.csv")
multiTOMdataDF <- read.csv(
  "../processed_data/multiTOM_allen_neighbors100_ALDH1A1_TH_SLC18A2_KCNJ6.csv")
multiTOMdataDF <- read.csv(
  "../processed_data/multiTOM_allen_neighbors30_ALDH1A1_TH_SLC18A2_KCNJ6.csv")
multiTOMdataDF <- read.csv(
  "../processed_data/multiTOM_allen_neighbors30_CALB1.csv")
multiTOMdataDF <- read.csv(
  "../processed_data/multiTOM_allen_neighbors100_CALB1.csv")

# Variable for minModSize parameter used in blockwiseModules WGCNA function
# to record in output graph titles
minModSize <- "30"

# bwModulesLL is list of modules from different blockwiseModules parameters used
# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 11

# Split by model into lists or dataframes of genes in each model
allenGenesModelsLL <- split(colnames(exprData)
                          , bwModulesLL[[modNetworkToUse]]$colors)

multiTOMgenesV <- multiTOMdataDF[ ,1]

intModulesLL <- sapply(allenGenesModelsLL
                   , function(AllenGenes) intersect(AllenGenes, multiTOMgenesV))

intAllenMultiTOMDF <- data.frame(sapply(allenGenesModelsLL, length))
intAllenMultiTOMDF <- cbind(intAllenMultiTOMDF
                            , data.frame(sapply(intModulesLL, length)))
colnames(intAllenMultiTOMDF) <- c("Allen", "Intersection")

write.csv(intAllenMultiTOMDF
          , file = "../analysis/MultiTOM_Allen_intersect_minMod30_ALDH1A1_TH_SLC18A2_KCNJ6_neighbors100.csv")


# bwModulesLL is list of modules from different blockwiseModules parameters used
# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 18

# Split by model into lists or dataframes of genes in each model
allenGenesModelsLL <- split(colnames(exprData)
                            , bwModulesLL[[modNetworkToUse]]$colors)

multiTOMgenesV <- multiTOMdataDF[ ,1]

intModulesLL <- sapply(allenGenesModelsLL
                   , function(AllenGenes) intersect(AllenGenes, multiTOMgenesV))

intAllenMultiTOMDF <- data.frame(sapply(allenGenesModelsLL, length))
intAllenMultiTOMDF <- cbind(intAllenMultiTOMDF
                            , data.frame(sapply(intModulesLL, length)))
colnames(intAllenMultiTOMDF) <- c("Allen", "Intersection")

write.csv(intAllenMultiTOMDF
          , file = "../analysis/MultiTOM_Allen_intersect_minMod100_ALDH1A1_TH_SLC18A2_KCNJ6_neighbors100.csv")




# bwModulesLL is list of modules from different blockwiseModules parameters used
# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 11

# Split by model into lists or dataframes of genes in each model
allenGenesModelsLL <- split(colnames(exprData)
                            , bwModulesLL[[modNetworkToUse]]$colors)

multiTOMgenesV <- multiTOMdataDF[ ,1]

intModulesLL <- sapply(allenGenesModelsLL
                       , function(AllenGenes) intersect(AllenGenes, multiTOMgenesV))

intAllenMultiTOMDF <- data.frame(sapply(allenGenesModelsLL, length))
intAllenMultiTOMDF <- cbind(intAllenMultiTOMDF
                            , data.frame(sapply(intModulesLL, length)))
colnames(intAllenMultiTOMDF) <- c("Allen", "Intersection")

write.csv(intAllenMultiTOMDF
          , file = "../analysis/MultiTOM_Allen_intersect_minMod30_neighbors100_CALB1.csv")

# bwModulesLL is list of modules from different blockwiseModules parameters used
# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 11

# Split by model into lists or dataframes of genes in each model
allenGenesModelsLL <- split(colnames(exprData)
                            , bwModulesLL[[modNetworkToUse]]$colors)

multiTOMgenesV <- multiTOMdataDF[ ,1]

intModulesLL <- sapply(allenGenesModelsLL
                       , function(AllenGenes) intersect(AllenGenes, multiTOMgenesV))

intAllenMultiTOMDF <- data.frame(sapply(allenGenesModelsLL, length))
intAllenMultiTOMDF <- cbind(intAllenMultiTOMDF
                            , data.frame(sapply(intModulesLL, length)))
colnames(intAllenMultiTOMDF) <- c("Allen", "Intersection")

write.csv(intAllenMultiTOMDF
          , file = "../analysis/MultiTOM_Allen_intersect_minMod30_neighbors100_CALB1.csv")

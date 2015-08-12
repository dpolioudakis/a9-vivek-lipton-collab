print("#######################################################################")
print("Starting allen-combine-probes.R script...")
sessionInfo()

library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

# Load dataframe of combined allen brain data
load("../processed_data/array_data_subset_rda")

# Load probe and gene name data
probesDataDF <- read.csv(
  "../raw_data/allen_data/178236545-2015-07-15/Probes.csv")

# Add gene symbol and entrez to array data data frame
probesArrayDataDF <- merge(arrayDataSubsetDF, probesDataDF[ ,c(1,4,6)]
                              , by.x="probe", by.y="probe_id")

# Format data for collapseRows fxn
probesArrayDataDF <- merge(probesDataDF[ ,c(1,4,6)], arrayDataSubsetDF 
                              , by.x="probe_id", by.y="probe")
rownames(probesArrayDataDF) <- probesArrayDataDF$probe_id
rowGroups <- probesArrayDataDF$gene_symbol
probesArrayDataDF <- probesArrayDataDF[ ,-c(1:3)]
rowID <- rownames(probesArrayDataDF)

# Average probe expression
# Outputs list, 1st object is the collapsed array data
collapseObjectLDF <- collapseRows(probesArrayDataDF
                                    , rowGroup=rowGroups
                                    , rowID=rowID
                                    , method="Average")
arrayDataSubsetAvgProbesDF <- collapseObjectLDF[[1]]

save(arrayDataSubsetAvgProbesDF, metaDataSubsetLDF,
     file="../processed_data/array_data_subset_avg_probes.rda")

print("End of allen-combine-probes.R script...")

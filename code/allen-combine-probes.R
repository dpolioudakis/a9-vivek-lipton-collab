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
load("../processed_data/array.data.subset.df.rda")

# Load probe and gene name data
probes.data.df <- read.csv("../raw_data/allen_data/178236545-2015-07-15/Probes.csv")

# Add gene symbol and entrez to array data data frame
probes.array.data.df <- merge(array.data.subset.df, probes.data.df[ ,c(1,4,6)]
                              , by.x="probe", by.y="probe_id")

# Format data for collapseRows fxn
probes.array.data.df <- merge(probes.data.df[ ,c(1,4,6)], array.data.subset.df 
                              , by.x="probe_id", by.y="probe")
rownames(probes.array.data.df) <- probes.array.data.df$probe_id
row.groups <- probes.array.data.df$gene_symbol
probes.array.data.df <- probes.array.data.df[ ,-c(1:3)]
row.ID <- rownames(probes.array.data.df)

# Average probe expression
# Outputs list, 1st object is the collapsed array data
collapse.object.ldf <- collapseRows(probes.array.data.df
                                    , rowGroup=row.groups
                                    , rowID=row.ID
                                    , method="Average")
array.data.subset.avg.probes.df <- collapse.object.ldf[[1]]

save(array.data.subset.avg.probes.df, meta.data.subset.ldf,
     file="../processed_data/array_data.subset.avg.probes.rda")

print("End of allen-combine-probes.R script...")

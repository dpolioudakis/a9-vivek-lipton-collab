print("#######################################################################")
print("Starting allen-cluster-compare-meta-data.R script...")
sessionInfo()

# Cluster samples and compare to meta data

load("../processed_data/array.data.subset.avg.probes.rda")

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
trait.data <- mapply(function(brain.meta.data, brain.id)
                        cbind(brain.meta.data, brain.id)
                        , subset.sn.meta.data, brain.ids, SIMPLIFY= FALSE)
# Combine list of meta data for each brain into one data frame
trait.data <- do.call(rbind, trait.data)

# Cluster samples
sampleTree2 = hclust(dist(expr.data), method = "average")

# Convert traits to a color representation: for numeric traits white means low, 
# red means high, grey means missing entry
traitColors= data.frame(labels2colors(trait.data[,1:7])
                        , numbers2colors(trait.data[,8:13])
                        , labels2colors(trait.data[,14]))
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(trait.data),
                    main = "Sample dendrogram and trait heatmap")

# Select specific traits
# Convert traits to a color representation: for numeric traits white means low, 
# red means high, grey means missing entry
traitColors= data.frame(labels2colors(trait.data[,c(2,5,14)])
                        , numbers2colors(trait.data[,c(3,8:10)]))
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(trait.data)[c(2,5,14,3,8:10)],
                    main = "Sample dendrogram and trait heatmap")

print("Ending allen-cluster-compare-meta-data.R script...")
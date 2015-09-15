# CACNA1D expression in Allen brain regions

# Input:
#   Allen brain atlas expression data with probe expression averaged

print("#######################################################################")
print("Starting allen-CACNA1D-expr-by-brain-region.R script...")
sessionInfo()

library(reshape2)
library(ggplot2)

load("../processed_data/array_data_subset_avg_probes.rda")
exprDataDF <- arrayDataSubsetAvgProbesDF

# Each list of metaDataSubsetLDF is a DF of metadata corresponding to each brain
# Brains and samples are listed in the order of the observations in
# arrayDataSubsetAvgProbesDF
# Make vector of structure acronyms in order of arrayDataSubsetAvgProbesDF
# columns
brainRegionV <- NULL
brainRegionV <- lapply(metaDataSubsetLDF
                       , function(x) c(brainRegionV, as.character(x$structure_acronym)))
brainRegionV <- unlist(brainRegionV)
brainRegionV <- as.factor(brainRegionV)

# Make column names of exprData brain region acronym
colnames(exprDataDF) <- brainRegionV

# Melt for ggplot2
exprDataDF <- melt(exprDataDF)
colnames(exprDataDF) <- c("gene", "region", "expression")
# Subset to CACNA1D expression
cacna1dExpr <- exprDataDF[exprDataDF$gene %in% "CACNA1D", ]

ggplot(cacna1dExpr, aes(x = as.factor(region), y = expression)) +
  geom_boxplot(aes(fill = region)) +
  guides(fill = FALSE) +
  ylab("Expression\n(Normalized Intensity)") +
  xlab("Brain Region") +
  labs(title = paste(
    "allen-CACNA1D-expr-by-brain-region.R"
    , "\nCACNA1D expression by brain region in Allen data"
    , sep = "")) +
  theme_grey(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(color = "black")) +
  theme(plot.title = element_text(size = rel(0.6)))
  theme(aspect.ratio = 4/4)
ggsave(file = "../analysis/CACNA1D expression brain region Allen data.pdf"
       , width = 5, height = 6)

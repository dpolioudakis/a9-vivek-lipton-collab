# CACNA1D expression by sample type

# Input:
#   Expression data for Lipton iPSC A9 cultures and human substantia nigra
#   samples, Vivek human cortex, Yuan iPSC neurons

# Output:
#   Boxplot of CACNA1D expression by sample type

sessionInfo()

library(reshape2)
library(ggplot2)

# Parameters
readDepthFilt <- 5
graphSubTitle <- paste("\nVivek Cortex, Yuan's iPSC neuron, Lipton's A9 and SN"
                       , "\nCQN normalized: GC, gene length, quantile"
                       , "\nRIN regressed out"
                       , "\nread depth filter: ", readDepthFilt
                       , "\nCACNA1D-expr-by-sample-type.R", sep = "")
outPathInfo <- paste("Cortex-iPSCneuron-LiptonA9sN_RDF", readDepthFilt
                     , "_CQN-geneLength-GC-quantile_OutlierRemoved"
                     , sep = "")

# Input file paths
load("../processed_data/HTSeqUnion_Gene_A9-SN-cortex-iPSCneuron_RDF5_CQN-geneLength-GC-quantile_OutlierRemoved_regRIN.rda")
exprDataDF <- as.data.frame(exprRegM)
metaDataDF <- metaDatDF

# Output file paths
outBoxPlotExpr <- "../analysis/CACNA1D expression by sample type.pdf"

# Subset expression to CACNA1D
cacna1dExpr <- exprDataDF[row.names(exprDataDF) %in% "ENSG00000157388", ]
# Add column with sample type labels
cacna1dExpr <- data.frame(expression = t(cacna1dExpr[1, ]), type = metaDataDF$Type)

ggplot(cacna1dExpr, aes(x = as.factor(type), y = ENSG00000157388)) +
  geom_boxplot(aes(fill = type)) +
  guides(fill = FALSE) +
  ylab("Expression\n(Normalized Intensity)") +
  xlab("Brain Region") +
  labs(title = paste("CACNA1D expression by sample type"
                     , graphSubTitle, sep = "")) +
  scale_x_discrete(labels = c("Cortex", "(2) High MEF2C", "(7) Low MEF2C"
                              , "iPSC neurons", "Human Substantia Nigra")) +
  theme_grey(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(color = "black")) +
  theme(plot.title = element_text(size = rel(0.6))) +
theme(aspect.ratio = 4/4)
ggsave(file = outBoxPlotExpr, width = 5, height = 6)

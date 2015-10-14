# Boxplot of A9 marker genes comparing Vivek's cortex samples,
# Yuan's iPSC neurons, and Lipton's A9 cultures and human substantia nigra
# samples

# Inputs:
#   CQN GC and read length normalized and read depth filtered HTseq expression
#   values

# Outputs:
#   Boxplot of log2 CQN normalized expression for each marker gene

################################################################################
sessionInfo()

library(reshape2)
library(ggplot2)
library(biomaRt)

# Parameters
outpathInfo <- " cortex-humanSN-iPSCa9 readDF5 CQN-GC-geneLength-quantile regRIN"
graphSubTitle <- paste("\nVivek Cortex, Lipton iPSC A9 and human SN"                       
                       , "\nread depth filter: 5"
                       , "\nCQN GC, gene length, quantile"
                       , "\nregressed out RIN"
                       , "\nA9-marker-expr-A9-SN-cortex-iPSCneuron.R", sep = "")
# For properly labeling samples
type <- as.factor(c(rep("Human_Cortex", 90)
                    , rep("2_High_MEF2C", 30)
                    , rep("Human_Substantia_Nigra", 50)
                    , rep("7_Low_MEF2C", 30)
                    , rep("Human_Substantia_Nigra", 40)
                    , rep("iPSC_Neuron", 80)))

# Input file paths (choose 1)
# Load CQN normalized expression values for hESC A9 and human substantia nigra
# samples
load("../processed_data/HTSeqUnion_Gene_A9-SN-cortex-iPSCneuron_RDF5_CQN-geneLength-GC-quantile_OutlierRemoved_regRIN.rda")
exprDataDF <- as.data.frame(exprRegM)

# Output file paths and variables
outpathMarks <- paste(
  "../analysis/A9 marker expression boxplot"
  , outpathInfo, ".pdf", sep="")

# Other variables
# Data frame of Ensembl IDs, gene symbols, and status as marker or anti-marker
markersDF <- data.frame(
  ensembl = c(
    "ENSG00000180176",
    "ENSG00000157542",
    "ENSG00000165646",
    "ENSG00000157388",
    "ENSG00000162761",
    "ENSG00000125798",
    "ENSG00000153234",
    "ENSG00000165092",
    
    "ENSG00000104327",
    "ENSG00000172137",
    "ENSG00000165588",
    "ENSG00000165655",
    "ENSG00000148680"
  )
  , gene = c(
    #markers
    "TH", "KCNJ7","VMAT2","CACNA1D","LMX1A","FOXA2","NURR1","ALDH1",
    #anti markers
    "CALB1", "CALB2", "OTX2", "NOLZ1", "HTR7"
    #"NOLZ1"(also known as ZNF503)
  )
  , marker = c(rep("marker",8),rep("anti-marker",5)))
################################################################################

# Subset expression data by marker genes present
markersExprDF <- exprDataDF[(rownames(exprDataDF) %in% markersDF[,1]),]

# Combine data frames to add gene symbol and marker information
markersExprDF <- merge(markersExprDF, markersDF
                       , by.x="row.names", by.y="ensembl")

# Reshape for ggplot2 using Reshape2 package
markersExprDF <- melt(markersExprDF, value.name="expression"
                      , variable.name = "sample"
                      , id = c("gene","marker", "Row.names"))

# Add column with sample type label
markersExprDF <- cbind(markersExprDF
                       , type = type)

# Plot marker genes expression
ggplot(markersExprDF, aes(x = gene, y = expression)) + 
  geom_boxplot(aes(fill = type)) +
  scale_fill_discrete(name = "Sample Type"
                       , labels = c("(2) High MEF2C", "(7) Low MEF2C"
                                    , "Human Cortex", "Human Substantia Nigra"
                                    , "iPSC neuron")) +
  xlab("Gene") +
  ylab("log2 Normalized expression") +
  labs(title = paste(
    "A9 markers expression", graphSubTitle, sep = "")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
ggsave(file = outpathMarks, height = 8)

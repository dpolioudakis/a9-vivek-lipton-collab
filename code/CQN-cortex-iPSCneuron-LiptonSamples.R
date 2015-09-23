# CQN normalization on Vivek's human cortex samples (controls from ASD study),
# Yuan's iPSC derived neurons, Lipton's iPSC A9 cultures, and Lipton's human
# substantia nigra samples

# Inputs:
#   HTseq
#     Samples were aligned with Tophat, FPKMs calculated with HTseq

# Outputs:
#   Histogram post CQN normalization of GC and gene length correlations
#   Table of normalized gene expression

# Load data
  # Subset Vivek cortex data to desired samples
# Merge into data frame
# Option for read depth filter
# CQN
# Diagnostics

library(cqn)
library(ggplot2)
library(reshape2)
library(WGCNA)

# Parameters
readDepthFilt <- 5
graphSubTitle <- paste("\nVivek Cortex, Yuan's iPSC neuron, Lipton's A9 and SN"
                       , "\nread depth filter: ", readDepthFilt
                       , "\nCQN-cortex-IPSCneuron-LiptonSamples.R", sep = "")
outPathInfo <- paste("Cortex-iPSCneuron-LiptonA9sN_RDF", readDepthFilt
                     , "_CQN-geneLength-GC-quantile_OutlierRemoved"
                     , sep = "")

# Load inputs
# Vivek cortex
load("../raw_data/human_cortex_from_vivek_ASD_project/AllData_N263_HTSC_filteredExpr_3-13-2014.Rdata")
cortexDatDF <- datExpr.HTSC.wholegene.CTX
cortexMetaDF <- read.csv("../raw_data/human_cortex_from_vivek_ASD_project/human_cortex_metadata_vivek.csv")

# Yuan iPSC neurons
neuronDatDF <- read.table("../raw_data/yuan_iPSC_neuron/Processed_data_htseqcount_83samples.txt", header = TRUE, row.names = 1)
neuronMetaDF <- read.table("../raw_data/yuan_iPSC_neuron/sample.info_CIRMGageLab_combined_109samples_112213.txt")

# Lipton iPSC A9 and human substantia nigra
a9sNdatDF <- read.csv("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnionGene.csv")
a9sNmetaDF <- read.csv("../Vivek_WGCNA_Lipton_A9_SN/metadata_SN.csv")

# For CQN normalization:
load("../Vivek_WGCNA_Lipton_A9_SN/GC18unionAnno.Rdata")
gc18unionAnno <- gc18unionAnno

# Output file paths and variables
outpathCQNhist <- paste("../analysis/CQN QC histogram "
  , outPathInfo, ".pdf", sep="")
outpathExprBoxplot <- paste("../analysis/CQN-post boxplot expression "
                            , outPathInfo, ".pdf", sep="")
outpathData <- paste("../processed_data/HTSeqUnion_Gene_"
                     , outPathInfo, ".rda", sep = "")
################################################################################

# Prepare input data

# Subset Vivek's data to 10 cortex control (non ASD) samples to use for further
# analysis
cortexDatDF <- cortexDatDF[ ,colnames(cortexDatDF) %in% cortexMetaDF$Sample.Name]

# Subset Yuan's data to iPSC neuron samples
neuronDatDF <- neuronDatDF[ ,46:53]
# Subset Yuan's metadata to iPSC neuron samples in expression data and
# reorganize
neuronMetaDF <- neuronMetaDF[neuronMetaDF$Disease.Status == "ctrl" 
                             & neuronMetaDF$cell.type == "Neuron", ]
neuronMetaDF <- rbind(
    subset(neuronMetaDF, patient.id == "cent" & clone.ID == "3-3")
  , subset(neuronMetaDF, patient.id == "cent" & clone.ID == "3-6")
  , subset(neuronMetaDF, patient.id == "chap" & clone.ID == "1-2")
  , subset(neuronMetaDF, patient.id == "chap" & clone.ID == "1-5")
  , subset(neuronMetaDF, patient.id == "clay" & clone.ID == "1-2")
  , subset(neuronMetaDF, patient.id == "clue" & clone.ID == "4-7")
  , subset(neuronMetaDF, patient.id == "cove" & clone.ID == "3-7")
  , subset(neuronMetaDF, patient.id == "cove" & clone.ID == "3-1")
)

# Format Lipton's A9 and SN dataframe
a9sNdatDF <- data.frame(a9sNdatDF[ ,-1], row.names = a9sNdatDF[ ,1])
colnames(a9sNdatDF) <- a9sNmetaDF$SampleID
# Remove empty rows of Lipton's metadata
a9sNmetaDF <- a9sNmetaDF[1:16, ]

# Merge Vivek's cortex data and Lipton's A9 and SN data
exprDatDF <- merge(cortexDatDF, a9sNdatDF, by = "row.names")
exprDatDF <- data.frame(exprDatDF[ ,-1], row.names = exprDatDF[ ,1])
# Remove .## suffix from ensembl IDs
rownames(exprDatDF) <- substring(rownames(exprDatDF), 1, 15)
# Merge in Yuan's iPSC neuron data
exprDatDF <- merge(exprDatDF, neuronDatDF, by = "row.names")
exprDatDF <- data.frame(exprDatDF[ ,-1], row.names = exprDatDF[ ,1])

# Remove .## suffix from ensembl IDs
rownames(gc18unionAnno) <- substring(rownames(gc18unionAnno), 1, 15)
################################################################################

FilterDepth <- function (exrDatDF) {
  # Filter genes - Keep: Genes with more than X counts in at least 80% of samples
  passV <- apply(exprDatDF, 1, quantile, 0.8) > readDepthFilt
  exprDatDF <- exprDatDF[passV,]
  exrDatDF
}
exprDatDF <- FilterDepth(exprDatDF)


# Use CQN to normalize for GC content and gene length
RunCQN <- function (exprDatDF) {
  # Remove genes not in GC and gene length annotation file
  keepV <- intersect(rownames(gc18unionAnno), rownames(exprDatDF))
  geneAnno <- gc18unionAnno[match(keepV, rownames(gc18unionAnno)), ]
  # Set genes with length 0 to length 1 - why? - from Vivek's code
  geneAnno[geneAnno[,1] == 0] <- 1
  exprDatDF <- exprDatDF[match(keepV, rownames(exprDatDF)), ]
  
  # Run CQN with specified depths and no quantile normalization
  cqnDat <- cqn(exprDatDF, lengths = as.numeric(geneAnno[,1])
                , x = as.numeric(geneAnno[,2]), lengthMethod = c("smooth")
                , sqn = TRUE)
  # Get the log2(Normalized FPKM) values
  cqnDat <- cqnDat$y + cqnDat$offset
  cqnDat
}
cqnDat <- RunCQN(exprDatDF)

# Check correlation of expression level to GC content and gene length pre and
# post CQN normalization
CheckCQNnorm <- function (preCQN, postCQN) {
  # Remove genes not in GC and gene length annotation file
  keepV <- intersect(rownames(gc18unionAnno), rownames(exprDatDF))
  geneAnno <- gc18unionAnno[match(keepV, rownames(gc18unionAnno)), ]
  # Set genes with length 0 to length 1 - why? - from Vivek's code
  geneAnno[geneAnno[,1] == 0] <- 1
  keepgenes <- intersect(rownames(preCQN),rownames(postCQN))
  preCQN <- preCQN[match(keepgenes,rownames(preCQN)),]
  postCQN <- postCQN[match(keepgenes,rownames(postCQN)),]
  geneAnno1 <- geneAnno[match(keepgenes,rownames(geneAnno)),]
  
  qcCorrCQNm <- matrix(NA,nrow=ncol(preCQN),ncol=4)
  colnames(qcCorrCQNm) <- c("preNormGCcorr", "preNormLengthCorr"
                        ,"postNormGCcorr", "postNormLengthCorr")
  for (i in 1:nrow(qcCorrCQNm)) {
    qcCorrCQNm[i,1] <- cor(preCQN[,i], geneAnno1[,2], method="spearman")
    qcCorrCQNm[i,2] <- cor(preCQN[,i], geneAnno1[,1], method="spearman")
    qcCorrCQNm[i,3] <- cor(postCQN[,i], geneAnno1[,2], method="spearman")
    qcCorrCQNm[i,4] <- cor(postCQN[,i], geneAnno1[,1], method="spearman")
  }
  qcCorrCQNm
}
qcCorrCQNm <- CheckCQNnorm(exprDatDF, cqnDat)
apply(qcCorrCQNm, 2, quantile)
qcCorrCQNm <- data.frame(qcCorrCQNm)
qcCorrCQNm <- melt(qcCorrCQNm)
colnames(qcCorrCQNm) <- c("CorrType", "Corr")

# Histogram of spearman's rho pre and post normalization for GC and gene length
ggplot(qcCorrCQNm, aes(x = Corr)) +
  facet_wrap(~CorrType, nrow = 2) +
  geom_histogram(binwidth = 0.01) +
  ylab("Counts") +
  xlab("Spearman's rho across samples") +
  labs(title = paste("Histogram: Spearman's rho across samples - pre and post CQN"
                     , graphSubTitle, sep = "")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = outpathCQNhist, height = 6)

# Boxplot of log2 expression post CQN
meltCQNdat <- melt(cqnDat)
ggplot(meltCQNdat, aes(y = value, x = Var2)) +
  geom_boxplot() +
  ylab("log2 Expression") +
  xlab("Samples") +
  labs(title = paste("Expression: Post CQN GC and gene length normalization"
                     , graphSubTitle, sep = "")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(file = outpathExprBoxplot, height = 6)

# Outlier Removal based on connectivity
sdOut <- 2
normAdj <- (0.5 + 0.5 * bicor(cqnDat)^2)

## Calculate connectivity
netSummary <- fundamentalNetworkConcepts(normAdj)
ku <- netSummary$Connectivity
zKu <- ku - (mean(ku)) / sqrt(var(ku))
# Declare as outliers those samples which are more than sdOut sd above the mean
# connectivity based on the chosen measure - this is basically filtering based
# on clustering, could use different type of clustering, Horvath found this was
# more robust
outliers <- (zKu > mean(zKu) + sdOut * sd(zKu))|(zKu < mean(zKu) - sdOut * sd(zKu))
print(paste("There are ", sum(outliers)
            , " outliers samples based on a bicor distance sample network"
            , "connectivity standard deviation above ", sdOut, sep = ""))
print(colnames(cqnDat)[outliers])
print(table(outliers))
cqnDatAll <- cqnDat
metaDatA9sNalldF <- a9sNmetaDF
cqnDat <- cqnDatAll[ ,!outliers]
# a9sNmetaDF <- metaDatA9sNalldF[!outliers, ]

save(cqnDat, a9sNmetaDF, cortexMetaDF, neuronMetaDF, file = outpathData)

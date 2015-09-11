# Boxplot of normalized expression by sample for Lipton hESC A9 cultures and
# human substantia nigra samples

library(reshape2)
library(ggplot2)
library(biomaRt)

load("../processed_data/allen_BW_modules.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells_5.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN_5.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells.rda")
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_humanSN.rda")

# variable for read depth filter to record in output graph titles
readDepthFilt <- "5"

# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
modNetworkToUse <- 11
modsToUse <- c("plum1", "grey60", "brown", "red", "cyan", "yellowgreen"
               , "sienna3", "royalblue")

# Data frame of A9 samples and human SN samples, each column is a sample
exprA9sNdF <- merge(datExpr.HTSC.A9, datExpr.HTSC.SN
                    , by.x = "row.names", by.y = "row.names")
exprA9sNdF <- data.frame(exprA9sNdF[ ,-1], row.names = exprA9sNdF[ ,1])

# Boxplot of normalized expression by sample
exprA9 <- melt(exprA9sNdF)
colnames(exprA9) <- c("sample", "expression")
ggplot(exprA9, aes(x = sample, y = expression)) +
  geom_boxplot() +
  labs(title = paste(
  "a9-human-substantia-nigra-expr-by-sample.R"
  , "\nLipton A9 and SN samples CQN normalized expression"
  , "\nread depth filter: ", readDepthFilt)
  , sep = "") +
  theme_grey(base_size = 14) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste(
  "../analysis/Lipton A9 and SN samples expression boxplot"
  , readDepthFilt, ".pdf", sep=""))
# Code adapted from Vivek's Cufflinks&HTseqCounts_LiptonSN.R

# Regress out rin from Vivek's Cortex and Lipton's A9 and human substantia
# nigra data and compare to principal components pre and post regression

# Inputs
#   CQN GC and read length normalized and read depth filtered HTseq expression
#   values Metadata with RIN values for each sample

# Outputs
#   Expression values with RIN and Total Reads regressed out
#   Correlation plots of principal compenents and technical confounders before
#     and after regression

sessionInfo()

library(boot)
library(WGCNA)
library(ggplot2)
library(reshape2)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

# Parameters
readDepthFilt <- 5
graphSubTitle <- paste("\nVivek Cortex, Yuan's iPSC neuron, Lipton's A9 and SN"
                       , "\nCQN normalized: GC, gene length, quantile"
                       , "\nRIN regressed out"
                       , "\nread depth filter: ", readDepthFilt
                       , "\nA9-SN-cortex-iPSCneuron-regression.R", sep = "")
outPathInfo <- paste("A9-SN-cortex-iPSCneuron_RDF", readDepthFilt
                     , "_CQN-geneLength-GC-quantile_OutlierRemoved"
                     , "_regRIN", sep = "")

# Load CQN normalized expression values for hESC A9 and human substantia nigra
# samples
load("../processed_data/HTSeqUnion_Gene_Cortex-iPSCneuron-LiptonA9sN_RDF5_CQN-geneLength-GC-quantile_OutlierRemoved.rda")
exprDatM <- cqnDat
a9sNmetaDF <- a9sNmetaDF 
cortexMetaDF <- cortexMetaDF
neuronMetaDF <- neuronMetaDF

# Output file paths
outPathCorr <- paste("../analysis/corr PCA confounders "
                     , outPathInfo, ".pdf", sep = "")
outPathReg <- paste("../processed_data/HTSeqUnion_Gene_"
                    , outPathInfo, ".rda", sep = "")
outPathExprBoxplot <- paste(
  "../analysis/Post confounders regression boxplot expression "
  , outPathInfo, ".pdf", sep = "")
outPathCorrReg <- paste("../analysis/corr PCA confounders regressed "
                        , outPathInfo, ".pdf", sep = "")
################################################################################

# Subset to desired metadata and combine metadata

# Adjust iPSC A9 culture names to make subsetting easier
a9sNmetaDF$SampleID <- as.character(a9sNmetaDF$SampleID)
a9sNmetaDF$SampleID[c(1:11)] <- paste(
  "X", a9sNmetaDF$SampleID[c(1:11)], sep = "")

# Subset and combine
subA9sNmetaDF <- a9sNmetaDF[ ,c("SampleID", "RIN.RQI", "Total_Reads"
                                , "Percent.Aligned", "Type")]
subA9sNmetaDF$Percent.Aligned <- subA9sNmetaDF$Total_Reads * subA9sNmetaDF$Percent.Aligned 
subCortexMetaDF <- cbind(
  cortexMetaDF[ ,c("Sample.Name", "RIN", "TotalReads.picard"
                   , "Aligned.Reads.picard")], Type = "cortex")
subNeuronMetaDF <- cbind(colnames(exprDatM)[25:32]
  , neuronMetaDF[ ,c("RNA.RIN", "sepDepth.PF1.indexMatch", "UniqueMapped")]
  , Type = "iPSCneuron")
mDatColNames <- c("SampleID", "RIN", "Total_Reads", "Aligned_Reads", "Type")
metaDatDF <- rbind(setNames(subA9sNmetaDF, mDatColNames)
                   , setNames(subCortexMetaDF, mDatColNames)
                   , setNames(subNeuronMetaDF, mDatColNames))
################################################################################

# PC Analysis

exprDataDF <- as.data.frame(t(exprDatM))
metaDatDF <- metaDatDF[metaDatDF$SampleID %in% row.names(exprDataDF) , ]

# TYH, KCNJ7,VMAT2,CACNA1D,LMX1A,FOXA2,NURR1,ALDH1
markers <- c("ENSG00000180176", "ENSG00000157542", "ENSG00000165646"
             , "ENSG00000157388", "ENSG00000162761", "ENSG00000125798"
             , "ENSG00000153234", "ENSG00000165092")
# CALB1/2, OTX2, NOLZ1,
antiMarkers <- c("ENSG00000104327", "ENSG00000172137", "ENSG00000165588"
                 , "ENSG00000165655", "ENSG00000148680")

markExprDF <- exprDataDF[ ,match(intersect(markers, colnames(exprDataDF))
                                 , colnames(exprDataDF))]
antiExprDF <- exprDataDF[ ,match(intersect(antiMarkers, colnames(exprDataDF))
                                 , colnames(exprDataDF))]

# Module Eigengene of markers and anti-markers
markersME <- as.numeric(as.matrix(
  moduleEigengenes(markExprDF
                   , colors = rep("red", ncol(markExprDF)), nPC=1)$eigengenes))
antiMarkME <- as.numeric(as.matrix(
  moduleEigengenes(antiExprDF
                   , colors = rep("red", ncol(antiExprDF)), nPC=1)$eigengenes))


# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(exprDatM), scale=F))
# Run PCA
pCdat <- prcomp(meanCenteredM, center=F);
topPCs <- pCdat$rotation[,1:5];
# Calculate variance explained by each PC
varExp <- (pCdat$sdev)^2 / sum(pCdat$sdev^2)
topVar <- varExp[1:5]
colnames(topPCs) <- paste("Unregressed\n", colnames(topPCs)
                         , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")

pairsDat <- data.frame(type = as.factor(metaDatDF$Type)
                     , totalReads = as.numeric(metaDatDF$Total_Reads)
                     , alignedReads = as.numeric(metaDatDF$Aligned_Reads)
                     , RIN = as.numeric(metaDatDF$RIN)
                     , A9.markers = markersME
                     , A9.antiMarkers = antiMarkME)

cond <- labels2colors(metaDatDF$Type)  ## colors

# Useful function for comparing multivariate data
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (class(x) == "numeric" & class(y) == "numeric") {
    r <- abs(cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$adj.r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pdf(outPathCorr, height = 20, width = 24)
pairs(cbind(pairsDat, topPCs), col = cond, pch = 19
      , upper.panel = panel.cor
      , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")
dev.off()
################################################################################

###For categorical data use this
# age=as.numeric(targets[7:16, ]$Age)
# sex <- as.numeric(factor(targets[7:16, ]$Sex))-1
# pmi=as.numeric(targets[7:16, ]$PMI..Autolysis.time..hr..)
# rin=as.numeric(targets[7:16, ]$RIN.RQI)
# regVars <- as.data.frame(cbind(age,sex,pmi,rin))

regVars <- data.frame(RIN = as.numeric(metaDatDF$RIN))

# 1) lme using covariates on unnormalized data
options(stringsAsFactors=FALSE)
boot <- FALSE
numboot <- 1000
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows bootstrap function to select samples
  fit <- lm(formula, data=d)
  return(coef(fit))
}  

# 2) lmboot for 106 samples + additional sample set (mega analysis - adding covariates in the model)
# Get the covariate data

# Regress out confounding variables
RegConf <- function (exprDatM) {
  exprDatRegM <- matrix(NA, nrow = nrow(exprDatM), ncol = ncol(exprDatM))
  rownames(exprDatRegM) <- rownames(exprDatM)
  colnames(exprDatRegM) <- colnames(exprDatM)
  ## change it to ncol(regVars)+1 when condition has 2 levels
  coefmat <- matrix(NA, nrow = nrow(exprDatM), ncol = ncol(regVars) + 1)
  if (boot==TRUE) {
    set.seed(8675309)
    for (i in 1:nrow(exprDatM)) {
      if (i%%1000 == 0) {print(i)}
      thisexp <- as.numeric(exprDatM[i, ])
      bs.results <- boot(data = data.frame(thisexp, regVars), statistic = bs,
                         R = numboot, formula = thisexp~. + age + sex + PMI)
      ## get the median - we can sometimes get NA values here... so let's exclude
      ## these - old code #bs.stats <- apply(bs.results$t,2,median) 
      bs.stats <- rep(NA, ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
      for (n in 1:ncol(bs.results$t)) {
        bs.stats[n] <- median(na.omit(bs.results$t[ ,n]))
      }
      coefmat[i,] <- bs.stats
      exprDatRegM[i,] <- thisexp - bs.stats[2]*regVars[,"age"] - bs.stats[3]*regVars[,"sex"] - bs.stats[4]*regVars[,"PMI"]
      # cat('Done for Gene',i,'\n')
    }
  } else {
    for (i in 1:nrow(exprDatM)) {
      if (i%%1000 == 0) {print(paste("Done for", i, "genes..."))}
      linMod <- lm(as.numeric(exprDatM[i, ]) ~ ., data = regVars)
      # The full data - the undesired covariates
      exprDatRegM[i,] <- coef(linMod)[1] + linMod$residuals
      # lmmod1 <- lm(as.numeric(exprDatM[i, ]) ~ condition + age + sex + pmi, data = regVars)
      ##datpred <- predict(object=lmmod1,newdata=regVars)
      # coef <- coef(lmmod1)
      # coefmat[i,] <- coef
      # The full data - the undesired covariates
      # exprDatRegM[i,] <- coef[1] + coef[2]*regVars[,"condition"] + lmmod1$residuals
      ## Also equivalent to <- thisexp - coef*var expression above
      # cat('Done for Genes',i,'\n')
    }
  }
  exprDatRegM
}
exprRegM <- RegConf(exprDatM)
quantile(exprDatM[ ,7],c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(exprRegM[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))
save(exprRegM, exprDatM, metaDatDF, file = outPathReg)

# Boxplot of log2 expression post CQN
meltRegExpr <- as.data.frame(melt(exprRegM))
ggplot(meltRegExpr, aes(y = value, x = Var2)) +
  geom_boxplot() +
  ylab("log2 Expression") +
  xlab("Samples") +
  labs(title = paste("Expression: Post Regression, CQN GC and gene length normalization"
                     , graphSubTitle, sep = "")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(file = outPathExprBoxplot)
################################################################################

#PC Analysis - Technical covariates - Regressed and Unregressed

# PC Analysis

# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(exprRegM), scale=F))
# Run PCA
pCdat <- prcomp(meanCenteredM, center=F);
topPCsReg <- pCdat$rotation[,1:5];
# Calculate variance explained by each PC
varExp <- (pCdat$sdev)^2 / sum(pCdat$sdev^2)
topVar <- varExp[1:5]
colnames(topPCsReg) <- paste("Regressed\n", colnames(topPCsReg)
                          , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")

pairsDat <- data.frame(type = as.factor(metaDatDF$Type)
                       , totalReads = as.numeric(metaDatDF$Total_Reads)
                       , AlignedReads = as.numeric(metaDatDF$Aligned_Reads)
                       , RIN = as.numeric(metaDatDF$RIN))

pdf(outPathCorrReg, height = 20, width = 24)
pairs(cbind(pairsDat, topPCs, topPCsReg), col = cond, pch = 19
      , upper.panel = panel.cor
      , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")
dev.off()

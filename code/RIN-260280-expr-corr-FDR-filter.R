# Filter genes based on correlation to RIN

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

# Load CQN normalized expression values for hESC A9 and human substantia nigra
# samples
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5.rda")
exprDataM <- datExpr.HTSC.A9SN[ ,1:6]
metadataDF <- targets.A9SN[1:6, ]
exprDataM <- datExpr.HTSC.A9SN
metadataDF <- targets.A9SN

# Filter genes based on RIN correlation FDR
corRIN <- cor(t(exprDataM), metadataDF$RIN)
corPvalRIN <- corPvalueStudent(corRIN, ncol(exprDataM))
corFDRrIN <- p.adjust(corPvalRIN, method = "fdr", n = length(corPvalRIN))
table(corFDRrIN > 0.1)
exprFDRfiltM <- exprDataM[corFDRrIN > 0.1, ]

# Filter genes based on Ratio 260/280 correlation FDR
cor260280 <- cor(t(exprFDRfiltM), metadataDF$Ratio_260.280)
corPval260280 <- corPvalueStudent(cor260280, ncol(exprFDRfiltM))
corFDr260280 <- p.adjust(corPval260280, method = "fdr", n = length(corPval260280))
table(corFDr260280 > 0.3)
exprFDRfiltM <- exprFDRfiltM[corFDr260280 > 0.3, ]

save(exprFDRfiltM, metadataDF
     , file = "../processed_data/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_FDrRIn260280.rda")





thisdat.HTSC <- t(scale(t(exprDataM),scale=F)) ## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC.HTSC <- prcomp(thisdat.HTSC,center=F);
topPC <- PC.HTSC$rotation[,1:5];
varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
topvar <- varexp[1:5]
colnames(topPC) <- paste("Unfiltered\n",colnames(topPC)," (",signif(100*topvar[1:5],2),"%)",sep="")

thisdat.HTSC <- t(scale(t(exprFDRfiltM),scale=F)) ## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC.HTSC <- prcomp(thisdat.HTSC,center=F);
topPCreg <- PC.HTSC$rotation[,1:5];
varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
topvar <- varexp[1:5]
colnames(topPCreg) <- paste("FDR < 0.01\n",colnames(topPCreg)," (",signif(100*topvar[1:5],2),"%)",sep="")

# pairsdat <- data.frame(totalReads=as.numeric(targets$Total_Reads),percentAligned=as.numeric(targets$Percent.Aligned),RIN=as.numeric(targets$RIN.RQI),A9.marker=markerME,A10.marker=antiMarkME)
pairsdat <- data.frame(  type = as.factor(metadataDF$Type)
                         , ratio260280 = as.numeric(metadataDF$Ratio_260.280)
                         , totalReads = as.numeric(metadataDF$Total_Reads)
                         , percentAligned = as.numeric(metadataDF$Percent.Aligned)
                         , RIN = as.numeric(metadataDF$RIN.RQI))

cond <- labels2colors(metadataDF$Type)  ## colors

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (class(x) == "numeric" & class(y) == "numeric") {
    r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$adj.r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(cbind(pairsdat, topPCreg), col = cond, pch = 19
      , upper.panel = panel.cor
      , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")
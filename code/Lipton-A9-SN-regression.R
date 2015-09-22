# Code adapted from Vivek's Cufflinks&HTseqCounts_LiptonSN.R

# Regress out rin and 260/280 from Lipton's A9 and human substantia
# nigra data and compare to principal components pre and post regression

# Adopted from Vivek's code

library(boot)
library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

# Load CQN normalized expression values for hESC A9 and human substantia nigra
# samples
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5.rda")
exprDataM <- datExpr.HTSC.A9SN
metadataDF <- targets.A9SN

###PC Analysis

exprDataDF <- as.data.frame(t(exprDataM))

marker=c("ENSG00000180176","ENSG00000157542","ENSG00000165646","ENSG00000157388","ENSG00000162761","ENSG00000125798","ENSG00000153234","ENSG00000165092") ##TYH, KCNJ7,VMAT2,CACNA1D,LMX1A,FOXA2,NURR1,ALDH1
anti.marker=c("ENSG00000104327","ENSG00000172137","ENSG00000165588","ENSG00000165655","ENSG00000148680") ##CALB1/2, OTX2, NOLZ1,

markExprDF <- exprDataDF[,match(intersect(marker,colnames(exprDataDF)),colnames(exprDataDF))]
antiExprDF <- exprDataDF[,match(intersect(anti.marker,colnames(exprDataDF)),colnames(exprDataDF))]

##Module Eigengene of marker and anti-marker
markerME <- as.numeric(as.matrix(
  moduleEigengenes(markExprDF
                   , colors = rep("red",ncol(markExprDF)), nPC=1)$eigengenes))
antiMarkME <- as.numeric(as.matrix(
  moduleEigengenes(antiExprDF
                   , colors = rep("red",ncol(antiExprDF)), nPC=1)$eigengenes))

pdf("corr_PCA_techs_A9_SN_RDF5.pdf",height=20,width=24)

thisdat.HTSC <- t(scale(t(exprDataM),scale=F)) ## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC.HTSC <- prcomp(thisdat.HTSC,center=F);
topPC <- PC.HTSC$rotation[,1:5];
varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
topvar <- varexp[1:5]
colnames(topPC) <- paste("Unregressed\n",colnames(topPC)," (",signif(100*topvar[1:5],2),"%)",sep="")

# pairsdat <- data.frame(totalReads=as.numeric(targets$Total_Reads),percentAligned=as.numeric(targets$Percent.Aligned),RIN=as.numeric(targets$RIN.RQI),A9.marker=markerME,A10.marker=antiMarkME)
pairsdat <- data.frame(  type = as.factor(metadataDF$Type)
                         , ratio260280 = as.numeric(metadataDF$Ratio_260.280)
                         , totalReads = as.numeric(metadataDF$Total_Reads)
                         , percentAligned = as.numeric(metadataDF$Percent.Aligned)
                         , RIN = as.numeric(metadataDF$RIN.RQI)
                         , A9.marker = markerME
                         , A9.anti.marker = antiMarkME)

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

pairs(cbind(pairsdat, topPC), col = cond, pch = 19
      , upper.panel = panel.cor
      , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")

dev.off()

### 1) lme using covariates on unnormalized data
options(stringsAsFactors=FALSE)
boot <- TRUE
numboot <- 20
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows bootstrap function to select samples
  fit <- lm(formula, data=d)
  return(coef(fit))
}  

### 2) lmboot for 106 samples + additional sample set (mega analysis - adding covariates in the model)
## Get the covariate data


###For categorical data use this
# age=as.numeric(targets[7:16, ]$Age)
# sex <- as.numeric(factor(targets[7:16, ]$Sex))-1
# pmi=as.numeric(targets[7:16, ]$PMI..Autolysis.time..hr..)
# rin=as.numeric(targets[7:16, ]$RIN.RQI)
# regvars <- as.data.frame(cbind(age,sex,pmi,rin))

regvars <- data.frame(  age = as.numeric(metadataDF[7:16, ]$Age)
                        , sex = as.numeric(factor(metadataDF[7:16, ]$Sex))-1
                        , PMI = as.numeric(metadataDF[7:16, ]$PMI..Autolysis.time..hr..))

## Run the regression
RegConf1 <- function (exprDataM) {
  exprDataRegM <- matrix(NA,nrow=nrow(exprDataM),ncol=ncol(exprDataM))
  rownames(exprDataRegM) <- rownames(exprDataM)
  colnames(exprDataRegM) <- colnames(exprDataM)
  coefmat <- matrix(NA,nrow=nrow(exprDataM),ncol=ncol(regvars)+1)## change it to ncol(regvars)+1 when condition has 2 levels
  if (boot==TRUE) {
    set.seed(8675309)
    for (i in 1:nrow(exprDataM)) {
      if (i%%1000 == 0) {print(i)}
      thisexp <- as.numeric(exprDataM[i,])
      bs.results <- boot(data=data.frame(thisexp,regvars), statistic=bs,
                         R=numboot, formula=thisexp~. +age+sex+PMI)
      ## get the median - we can sometimes get NA values here... so let's exclude
      ## these - old code #bs.stats <- apply(bs.results$t,2,median) 
      bs.stats <- rep(NA,ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
      for (n in 1:ncol(bs.results$t)) {
        bs.stats[n] <- median(na.omit(bs.results$t[,n]))
      }
      coefmat[i,] <- bs.stats
      exprDataRegM[i,] <- thisexp - bs.stats[2]*regvars[,"age"] - bs.stats[3]*regvars[,"sex"] - bs.stats[4]*regvars[,"PMI"]
      # cat('Done for Gene',i,'\n')
    }
  } else {
    for (i in 1:nrow(exprDataM)) {
      if (i%%1000 == 0) {print(i)}
      lmmod1 <- lm(as.numeric(exprDataM[i,])~condition+age+sex+pmi,data=regvars)
      ##datpred <- predict(object=lmmod1,newdata=regvars)
      coef <- coef(lmmod1)
      coefmat[i,] <- coef
      exprDataRegM[i,] <- coef[1] + coef[2]*regvars[,"condition"] + lmmod1$residuals ## The full data - the undesired covariates
      ## Also equivalent to <- thisexp - coef*var expression above
      cat('Done for Genes',i,'\n')
    }
  }
  exprDataRegM
}
exprRegSnM <- RegConf(exprDataM[ ,7:16])
save(exprRegSnM, exprDataM, metadataDF, datExpr.HTSC.A9SN
     , file = "../processed_data/HTSeqUnion_Exon_CQN_OutlierRemoved_SN_RDF5_regAgeSexPMI.rda")
quantile(exprDataM[ ,7],c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(exprRegSnM[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))

exprA9regSN <- cbind(exprDataM[ ,1:6], exprRegSnM)

regvars <- data.frame(  RIN = as.numeric(metadataDF$RIN.RQI)
                      , ratio260280 = as.numeric(metadataDF$Ratio_260.280))
regvars <- data.frame(  RIN = as.numeric(metadataDF$RIN.RQI))

## Run the regression
RegConf2 <- function (exprDataM) {
  exprDataRegM <- matrix(NA,nrow=nrow(exprDataM),ncol=ncol(exprDataM))
  rownames(exprDataRegM) <- rownames(exprDataM)
  colnames(exprDataRegM) <- colnames(exprDataM)
  coefmat <- matrix(NA,nrow=nrow(exprDataM),ncol=ncol(regvars)+1)## change it to ncol(regvars)+1 when condition has 2 levels
  if (boot==TRUE) {
    set.seed(8675309)
    for (i in 1:nrow(exprDataM)) {
      if (i%%1000 == 0) {print(i)}
      thisexp <- as.numeric(exprDataM[i,])
      bs.results <- boot(data=data.frame(thisexp,regvars), statistic=bs,
                         R=numboot, formula=thisexp~. +ratio260280+RIN)
      ## get the median - we can sometimes get NA values here... so let's exclude
      ## these - old code #bs.stats <- apply(bs.results$t,2,median) 
      bs.stats <- rep(NA,ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
      for (n in 1:ncol(bs.results$t)) {
        bs.stats[n] <- median(na.omit(bs.results$t[,n]))
      }
      coefmat[i,] <- bs.stats
      exprDataRegM[i,] <- thisexp - bs.stats[2]*regvars[,"RIN"] - bs.stats[3]*regvars[,"ratio260280"]
      # cat('Done for Gene',i,'\n')
    }
  } else {
    for (i in 1:nrow(exprDataM)) {
      if (i%%1000 == 0) {print(i)}
      lmmod1 <- lm(as.numeric(exprDataM[i,])~condition+age+sex+pmi,data=regvars)
      ##datpred <- predict(object=lmmod1,newdata=regvars)
      coef <- coef(lmmod1)
      coefmat[i,] <- coef
      exprDataRegM[i,] <- coef[1] + coef[2]*regvars[,"condition"] + lmmod1$residuals ## The full data - the undesired covariates
      ## Also equivalent to <- thisexp - coef*var expression above
      cat('Done for Genes',i,'\n')
    }
  }
  exprDataRegM
}
exprReg <- RegConf2(exprA9regSN)
exprDataRegM <- exprReg
quantile(exprDataM[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(exprReg[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))

save(exprDataRegM, exprDataM, metadataDF, datExpr.HTSC.A9SN
     , file = "../processed_data/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_regRIN.rda")
save(exprDataRegM, exprDataM, metadataDF, datExpr.HTSC.A9SN
     , file = "../processed_data/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_regAgeSexPMiRIn260280.rda")
load(file = "../processed_data/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5_regRIN260280.rda")

quantile(exprDataM[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(exprDataRegM[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))
print("#######################################################################")


#PC Analysis - Technical covariates - Regressed and Unregressed

pdf("../analysis/PCA_A9_SN_RDF5_regressed_RIN.pdf",height=20,width=24)
pdf("../analysis/corr_PCA_techs_A9_SN_RDF5_regressed_RIN_260280_age_sex_PMI.pdf",height=20,width=24)

thisdat.HTSC <- t(scale(t(exprDataM),scale=F)) ## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC.HTSC <- prcomp(thisdat.HTSC,center=F);
topPC <- PC.HTSC$rotation[,1:5];
varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
topvar <- varexp[1:5]
colnames(topPC) <- paste("Unregressed\n",colnames(topPC)," (",signif(100*topvar[1:5],2),"%)",sep="")

thisdat.HTSC <- t(scale(t(exprDataRegM),scale=F)) ## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC.HTSC <- prcomp(thisdat.HTSC,center=F);
topPCreg <- PC.HTSC$rotation[,1:5];
varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
topvar <- varexp[1:5]
colnames(topPCreg) <- paste("Regressed\n",colnames(topPCreg)," (",signif(100*topvar[1:5],2),"%)",sep="")

# pairsdat <- data.frame(totalReads=as.numeric(targets$Total_Reads),percentAligned=as.numeric(targets$Percent.Aligned),RIN=as.numeric(targets$RIN.RQI),A9.marker=markerME,A10.marker=antiMarkME)
pairsdat <- data.frame(  type = as.factor(metadataDF$Type)
                         , ratio260280 = as.numeric(metadataDF$Ratio_260.280)
                         , totalReads = as.numeric(metadataDF$Total_Reads)
                         , percentAligned = as.numeric(metadataDF$Percent.Aligned)
                         , RIN = as.numeric(metadataDF$RIN.RQI))

cond <- labels2colors(metadataDF$Type)  ## colors

pairs(cbind(pairsdat, topPC, topPCreg), col = cond, pch = 19
      , upper.panel = panel.cor
      , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")

dev.off()


#PC Analysis - Marker MEs - Regressed

pdf("../analysis/corr_PCA_techs_markers_A9_SN_RDF5_regressed_RIN_260280_age_sex_PMI.pdf",height=20,width=24)

exprDataDF <- as.data.frame(t(exprDataRegM))

marker=c("ENSG00000180176","ENSG00000157542","ENSG00000165646","ENSG00000157388","ENSG00000162761","ENSG00000125798","ENSG00000153234","ENSG00000165092") ##TYH, KCNJ7,VMAT2,CACNA1D,LMX1A,FOXA2,NURR1,ALDH1
anti.marker=c("ENSG00000104327","ENSG00000172137","ENSG00000165588","ENSG00000165655","ENSG00000148680") ##CALB1/2, OTX2, NOLZ1,

markExprRegDF <- exprDataDF[,match(intersect(marker,colnames(exprDataDF)),colnames(exprDataDF))]
antiExprRegDF <- exprDataDF[,match(intersect(anti.marker,colnames(exprDataDF)),colnames(exprDataDF))]

##Module Eigengene of marker and anti-marker
markerME <- as.numeric(as.matrix(
  moduleEigengenes(markExprRegDF
                   , colors = rep("red",ncol(markExprRegDF)), nPC=1)$eigengenes))
antiMarkME <- as.numeric(as.matrix(
  moduleEigengenes(antiExprRegDF
                   , colors = rep("red",ncol(antiExprRegDF)), nPC=1)$eigengenes))

# pairsdat <- data.frame(totalReads=as.numeric(targets$Total_Reads),percentAligned=as.numeric(targets$Percent.Aligned),RIN=as.numeric(targets$RIN.RQI),A9.marker=markerME,A10.marker=antiMarkME)
pairsdat <- data.frame(  type = as.factor(metadataDF$Type)
                         , ratio260280 = as.numeric(metadataDF$Ratio_260.280)
                         , totalReads = as.numeric(metadataDF$Total_Reads)
                         , percentAligned = as.numeric(metadataDF$Percent.Aligned)
                         , RIN = as.numeric(metadataDF$RIN.RQI)
                         , A9.marker = markerME
                         , A9.anti.marker = antiMarkME)

cond <- labels2colors(metadataDF$Type)  ## colors

pairs(cbind(pairsdat, topPCreg), col = cond, pch = 19
      , upper.panel = panel.cor
      , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")
dev.off()
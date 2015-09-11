# Correlation of PCA of Lipton hESC A9 cultures and human substantia nigra
# samples data to RNAseq technical measures and A9 markers

###PC Analysis

inputFile <- "../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9_SN_RDF5.rda"
load(inputFile)
exprDataDF <- datExpr.HTSC.A9SN
metadataDF <- targets.A9SN
load("../processed_data/allen_BW_modules.rda")
# 11 corresponds to softPower 7, minModSize 30, deepSplit 2,
# MEmergeCutHeight 0.25, maxBlockSize 12000
genesColorsDF <- data.frame(colnames(exprData)
                            , bwModulesLL[[11]]$colors)
colnames(genesColorsDF) <- c("gene", "module")

print("#######################################################################")

metadataDF <- metadataDF[ ,c("Type", "RIN.RQI", "Ratio_260.280", "Total_Reads", "Percent.Aligned")]

exprDataDF <- as.data.frame(t(exprDataDF))

# Marker genes

marker <- c( "ENSG00000180176","ENSG00000157542","ENSG00000165646"
            ,"ENSG00000157388","ENSG00000162761","ENSG00000125798"
            ,"ENSG00000153234","ENSG00000165092") ##TYH, KCNJ7,VMAT2,CACNA1D,LMX1A,FOXA2,NURR1,ALDH1
antiMarker <- c( "ENSG00000104327","ENSG00000172137","ENSG00000165588"
                ,"ENSG00000165655","ENSG00000148680") ##CALB1/2, OTX2, NOLZ1,
 
exprDataDF.marker=exprDataDF[,match(intersect(marker,colnames(exprDataDF)),colnames(exprDataDF))]
exprDataDF.antiMarker=exprDataDF[,match(intersect(antiMarker,colnames(exprDataDF)),colnames(exprDataDF))]

##Module Eigengene of marker and anti-marker
ME.marker=as.numeric(as.matrix(moduleEigengenes(exprDataDF.marker, colors = rep("red",ncol(exprDataDF.marker)), nPC=1)$eigengenes))
ME.antiMarker=as.numeric(as.matrix(moduleEigengenes(exprDataDF.antiMarker, colors = rep("red",ncol(exprDataDF.antiMarker)), nPC=1)$eigengenes))

# Marker module genes

AddEnsembl <- function (geneList) {
  geneListDF <- data.frame(geneList)
  # bioMart manual:
  # http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  moduleEnsemblDF <- getBM(  attributes = c("ensembl_gene_id", "hgnc_symbol")
                             , filters = "hgnc_symbol"
                             , values = geneListDF
                             , mart = ensembl)
  moduleEnsemblDF
}

# Add Ensembl gene ID to gene symbol and module color data frame
genesEnsemblDF <- AddEnsembl(colnames(exprData))
ensemblColorsDF <- merge(genesEnsemblDF, genesColorsDF
                         , by.x = "hgnc_symbol", by.y = "gene")
# Merge with Lipton expression values
modsExprA9sNdF <- merge(ensemblColorsDF, t(exprDataDF)
                        , by.x = "ensembl_gene_id", by.y = "row.names")

moduleME <- moduleEigengenes(t(modsExprA9sNdF[ ,4:ncol(modsExprA9sNdF)])
                   , colors = modsExprA9sNdF$module[ ,drop=TRUE], nPC=1)$eigengenes
print("#######################################################################")

# Panel correlation plots of PCA and A9 and SN
pdf("PCA_metadata_corr_A9_SN_RDF5_CQNtogether.pdf",height=20,width=24)
thisdat.HTSC <- t(scale(exprDataDF,scale=F)) ## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC.HTSC <- prcomp(thisdat.HTSC,center=F);
topPC <- PC.HTSC$rotation[,1:5];
varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
topvar <- varexp[1:5]
colnames(topPC) <- paste(colnames(topPC)," (",signif(100*topvar[1:5],2),"%)",sep="")

pairsdat <- cbind(
          data.frame( type = as.factor(metadataDF$Type)
                    , ratio260280 = as.numeric(metadataDF$Ratio_260.280)
                    , totalReads = as.numeric(metadataDF$Total_Reads)
                    , percentAligned = as.numeric(metadataDF$Percent.Aligned)
                    , RIN = as.numeric(metadataDF$RIN.RQI)
                    , A9.marker = ME.marker
                    , A9.antiMarker = ME.antiMarker)
                  , plum1 = moduleME$MEplum1)

cond=labels2colors(metadataDF$Type)  ## colors

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs",method="spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (class(x) == "numeric" & class(y) == "numeric") {
    r <- abs(cor(x, y,use="pairwise.complete.obs",method="spearman"))
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
print("#######################################################################")

# Panel correlation plots of PCA and A9 and SN

pdf("PCA_metadata_corr_A9_RDF5_CQNtogether.pdf",height=20,width=24)
thisdat.HTSC <- t(scale(exprDataDF[1:6, ],scale=F)) ## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC.HTSC <- prcomp(thisdat.HTSC,center=F);
topPC <- PC.HTSC$rotation[,1:5];
varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
topvar <- varexp[1:5]
colnames(topPC) <- paste("Regressed\n",colnames(topPC)," (",signif(100*topvar[1:5],2),"%)",sep="")

pairsdat <- cbind(
  data.frame( type = as.factor(metadataDF$Type)
              , ratio260280 = as.numeric(metadataDF$Ratio_260.280)
              , totalReads = as.numeric(metadataDF$Total_Reads)
              , percentAligned = as.numeric(metadataDF$Percent.Aligned)
              , RIN = as.numeric(metadataDF$RIN.RQI)
              , A9.marker = ME.marker
              , A9.antiMarker = ME.antiMarker)
  , plum1 = moduleME$MEplum1)
pairsdat <- pairsdat[1:6, ]
moduleME <- moduleME[1:6, ]

cond=labels2colors(metadataDF$Type)  ## colors

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs",method="spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (class(x) == "numeric" & class(y) == "numeric") {
    r <- abs(cor(x, y,use="pairwise.complete.obs",method="spearman"))
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

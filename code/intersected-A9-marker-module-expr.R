# Expression Ratio of intersected Allen and Lipton A9 marker modules in Lipton
# hESC A9 cells

library(biomaRt)

load("../processed_data/allen_BW_modules.rda")
load("../processed_data/array_data_subset_avg_probes.rda")
# Vivek normalized RNAseq FPKMs from hESC derived A9 neuronal cultures
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells.rda")
# Loads intModsGenesLL and moduleCombos
# intModsGenesLL is lists of genes that is intersection of each Allen Lipton
# Human SN module combination.
# List is in order of module combinations listed in moduleCombos
load("../processed_data/allen_human_SN_modules_intersect_modsize100.rda")
minModSize <- 100
modCombos <- rbind(
    c("brown", "cyan")
  , c("brown", "darkred")
  , c("brown", "royalblue")
  , c("purple", "darkgreen")
  , c("purple", "white")
  , c("greenyellow", "magenta")
  , c("greenyellow", "orange")
)

load("../processed_data/allen_human_SN_modules_intersect_modsize30.rda")
minModSize <- 30
modCombos <- rbind(
    c("plum1", "skyblue")
  , c("brown", "turquoise")
  , c("brown", "yellow")
  , c("cyan", "turquoise")
  , c("cyan", "yellow")
)

SelectMod <- function (allenMod, liptonMod) {
  intModsGenesLL[moduleCombos[ ,1] == allenMod
                 & moduleCombos[ ,2] == liptonMod]
}

AddEnsembl <- function (geneList) {
  moduleGenes <- data.frame(geneList)
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
                             , values = moduleGenes
                             , mart = ensembl)
  moduleEnsemblDF
}

SubsetMarkerModInA9 <- function (moduleEnsemblDF) {
  # Subset Vivek produced expression data by marker module genes present
  markersExprDF <- datExpr.HTSC.A9[
    (rownames(datExpr.HTSC.A9) %in% moduleEnsemblDF[,1]), ]
  # Combine data frames to add gene symbol and module information
  markersExprDF <- merge(markersExprDF, moduleEnsemblDF
                         , by.x="row.names", by.y=1)
  # Changing column name like this is not working
  colnames(markersExprDF[1]) <- "ensembl"
  # Reshape for ggplot2 using Reshape2 package
  markersExprDF <- melt(markersExprDF
                        , value.name="expression", variable.name="sample")
  # Add column with treatment group label
  numGenesInTreatmentGroup <- (nrow(markersExprDF) / 2)
  markersExprDF <- cbind(markersExprDF
                         , bio.rep = as.factor(
                           c(rep(2,numGenesInTreatmentGroup)
                             , rep(7,numGenesInTreatmentGroup))))
  markersExprDF
}

geneListsLL <- apply(modCombos, 1, function(modCombo) SelectMod(modCombo[1], modCombo[2]))
geneListsLDF <- lapply(geneListsLL, AddEnsembl)
# Add gene expression (normalized FPKM) values to gene list
# Gene list is the intersection of Allen module and Lipton Human SN module
intModExprLDF <- lapply(geneListsLDF, SubsetMarkerModInA9)
names(intModExprLDF) <- apply(modCombos, 1
                          , function(modCombo) paste(modCombo[1], modCombo[2]))

# Calculate ratio of expression in high MEF2C samples versus low MEF2C samples
# for each gene in each module
# List (modules) of lists (each gene in that module)
ratioExprLL <- lapply(intModExprLDF, 
                      function(x) {
                        markExpr <- dcast(x, Row.names~sample, value.var="expression")
                        row.names(markExpr) <- markExpr[ ,1]
                        markExpr <- markExpr[ ,-1]
                        apply(markExpr, 1, function(x) (sum(x[1:3]) / sum(x[4:6])))
                      }
)
# Mean ratio of expression in each module
sapply(ratioExprLL, mean)
# Paired T-test of high MEF2C versus low MEF2C expression (not ratios) for each
# module
sapply(intModExprLDF, function(module) {
  highMEF2C <- module[module$bio.rep=="2", ]$expression
  lowMEF2C <-  module[module$bio.rep=="7", ]$expression
  round(t.test(highMEF2C, lowMEF2C, paired = TRUE)$p.value, 5)
})

# Boxplot of expression ratios for each module
ratioExprDF <- melt(ratioExprLL)
colnames(ratioExprDF) <- c("ratio.expr", "module")
# Preserve module order in boxplot
ratioExprDF$module <- factor(ratioExprDF$module
                             , levels = as.character(unique(ratioExprDF$module)))
ggplot(data = ratioExprDF, aes(x=module, y=ratio.expr)) + 
  geom_boxplot() +
  # geom_boxplot(aes(fill=module)) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(title = paste(
    "intersected-A9-marker-module-expr.R\nExpression Ratio of intersected Allen"
    , "and Lipton A9 marker modules in Lipton hESC A9 cells\nAllen minModSize: "
    , minModSize), sep = "") +
  ylab("Expression Ratio in high MEF2C versus low MEF2C") +
  xlab("Modules Intersected") +
  theme_grey(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  ggsave(file = paste(
    "../analysis/Intersected Modules Expression Ratio"
    , " ModSize", minModSize, ".pdf", sep=""))  #  "-", Sys.Date(),
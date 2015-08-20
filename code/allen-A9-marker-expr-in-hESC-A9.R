# Expression of A9 module 28 genes in hESC derived A9 neuronal cultures

print("#######################################################################")
print("Starting allen-A9-marker-expr-in-hESC-A9.R script...")
sessionInfo()

library(WGCNA)
library(reshape2)
library(ggplot2)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/allen_BW_modules.rda")
load("../processed_data/array_data_subset_avg_probes.rda")
# Vivek normalized RNAseq FPKMs from hESC derived A9 neuronal cultures
load("../Vivek_WGCNA_Lipton_A9_SN/HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells.rda")

# bwModulesLL is list of modules from different blockwiseModules parameters used
modulesToUse <- 15
blockwiseMEs <- moduleEigengenes(exprData
                                 , bwModulesLL[[modulesToUse]]$colors)$eigengenes

# Write table of gene names in module 28 (ME correlated with A9 markers)
geneModuleMembership <- as.data.frame(cor(exprData, blockwiseMEs, use = "p"))
module=28
moduleGenes <- bwModulesLL[[modulesToUse]]$colors==module
geneModuleMembership$ME28[moduleGenes]
row.names(geneModuleMembership[moduleGenes, ])
module28Genes <- data.frame(row.names(geneModuleMembership[moduleGenes, ]))
write.table(module28Genes, "../analysis/SN_module_genes.txt")

# Data frame of Ensembl IDs, gene symbols, and status as marker or anti-marker
# Start by making list of Enzembl IDs and gene symbols
module28Ensembl <- c("ENSG00000165092"=	"ALDH1A1"
                     ,"ENSG00000164512"=	"ANKRD55"
                     ,"ENSG00000186897"=	"C1QL4"
                     ,"ENSG00000205856"=	"C22orf42"
                     ,"ENSG00000004948"=	"CALCR"
                     ,"ENSG00000075275"=	"CELSR1"
                     ,"ENSG00000147432"=	"CHRNB3"
                     ,"ENSG00000167600"=	"CYP2S1"
                     ,"ENSG00000132437"=	"DDC"
                     ,"ENSG00000185559"=	"DLK1"
                     ,"ENSG00000108001"=	"EBF3"
                     ,"ENSG00000163064"=	"EN1"
                     ,"ENSG00000129514"=	"FOXA1"
                     ,"ENSG00000102287"=	"GABRE"
                     ,"ENSG00000131979"=	"GCH1"
                     ,"ENSG00000087460"=	"GNAS"
                     ,"ENSG00000137252"=	"HCRTR2"
                     ,"ENSG00000150361"=	"KLHL1"
                     ,"ENSG00000259974"=	"LINC00261"
                     ,"ENSG00000136944"=	"LMX1B"
                     ,"ENSG00000167419"=	"LPO"
                     ,"ENSG00000109805"=	"NCAPG"
                     ,"ENSG00000138653"=	"NDST4"
                     ,"ENSG00000168743"=	"NPNT"
                     ,"ENSG00000065320"=	"NTN1"
                     ,"ENSG00000101188"=	"NTSR1"
                     ,"ENSG00000183395"=	"PMCH"
                     ,"ENSG00000136546"=	"SCN7A"
                     ,"ENSG00000169432"=	"SCN9A"
                     ,"ENSG00000115884"=	"SDC1"
                     ,"ENSG00000145248"=	"SLC10A4"
                     ,"ENSG00000036565"=	"SLC18A1"
                     ,"ENSG00000165646"=	"SLC18A2"
                     ,"ENSG00000142319"= "SLC6A3"
                     ,"ENSG00000276996"=	"SLC6A3"
                     ,"ENSG00000159167"=	"STC1"
                     ,"ENSG00000169836"=	"TACR3"
                     ,"ENSG00000160180"=	"TFF3"
                     ,"ENSG00000180176"=	"TH"
                     ,"ENSG00000182223"=	"ZAR1")

# Convert to data frame
module28EnsemblDF <- data.frame(names(module28Ensembl)
                                , module28Ensembl, row.names = NULL)

# Subset Vivek produced expression data by marker genes present
markersExprDF <- datExpr.HTSC.A9[
  (rownames(datExpr.HTSC.A9) %in% module28EnsemblDF[,1]), ]

# Combine data frames to add gene symbol and marker information
markersExprDF <- merge(markersExprDF, module28EnsemblDF
                       , by.x="row.names", by.y=1)
# Changing column name like this is not working
colnames(markersExprDF[1]) <- "ensembl"


# Reshape for ggplot2 using Reshape2 package
markersExprDF <- melt(markersExprDF
                      , value.name="expression", variable.name="sample")

# Add column with biological replicate label
markersExprDF <- cbind(markersExprDF
                       , bio.rep = as.factor(c(rep(2,36), rep(7,36))))


# Plot expression marker genes high MEF2C and low MEF2C samples
ggplot(markersExprDF, aes(x=module28Ensembl, y=expression)) + 
  geom_boxplot(aes(fill=bio.rep)) +
  scale_fill_discrete(name= "Biological\nReplicate",
                      labels= c("2_(high MEF2C)", "7_(low MEF2C)")) +
  labs(title= "Allen derived A9 marker expression") +
  ylab("Expression (normalized FPKM)") +
  xlab("Gene") +
  theme_grey(base_size = 18)
ggsave(file= "A9 marker gene expression.pdf")

# Expression of A9 module 28 genes in hESC derived A9 neuronal cultures

print("#######################################################################")
print("Starting allen-A9-marker-expr-in-hESC-A9.R script...")
sessionInfo()

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../processed_data/allen_BW_modules.rda")
load("../processed_data/array_data_subset_avg_probes.rda")

# bwModules is list of modules from 3 different ME merge cut heights
blockwiseMEs <- moduleEigengenes(exprData, bwModules[[3]]$colors)$eigengenes


# Expression of A9 module 28 genes in hESC derived A9 neuronal cultures


geneModuleMembership <- as.data.frame(cor(exprData, blockwiseMEs, use = "p"))
module=28
moduleGenes <- bwModules[[3]]$colors==module
geneModuleMembership$ME28[moduleGenes]
row.names(geneModuleMembership[moduleGenes, ])
module28Genes <- data.frame(row.names(geneModuleMembership[moduleGenes, ]))
write.table(module28Genes, "../analysis/SN_module_genes.txt")



library(reshape2)
library(ggplot2)

setwd(".")
load("../HTSeqUnion_Exon_CQN_OutlierRemoved_A9cells.rda")
```

Data frame of Ensembl IDs, gene symbols, and status as marker or anti-marker
```{r marker_ensembl_labels}
markers.df <- data.frame(
  ensembl=
    c(
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
    ),
  gene=
    c(
      #markers
      "TH", "KCNJ7","VMAT2","CACNA1D","LMX1A","FOXA2","NURR1","ALDH1",
      #anti markers
      "CALB1", "CALB2", "OTX2", "NOLZ1", "HTR7"
      #"NOLZ1"(also known as ZNF503)
    ),
  marker.type=
    c(rep("marker",8),rep("anti-marker",5))
)
```

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



module28EnsemblDF <- data.frame(names(module28Ensembl), module28Ensembl, row.names = NULL)

Subset expression data by marker genes present
```{r subset_markers}
markers.exp.df <- datExpr.HTSC.A9[(rownames(datExpr.HTSC.A9) %in% module28EnsemblDF[,1]),]
```
Combine data frames to add gene symbol and marker information
```{r combine_marker_expression_dfs}
markers.exp.df <- merge(markers.exp.df, module28EnsemblDF, by.x="row.names", by.y=1)
colnames(markers.exp.df[1]) <- "ensembl"
colnames(markers.exp.df[7])
```
Reshape for ggplot2 using Reshape2 package
```{r reshape_marker_expression}
markers.exp.df <- melt(markers.exp.df, value.name="expression", variable.name="sample")
```
Add column with biological replicate label
```{r add_bio_rep_label}
markers.exp.df <- cbind(markers.exp.df, bio.rep = as.factor(c(rep(2,36), rep(7,36))))
```

Plot marker genes expression
```{r plot_markers}
ggplot(markers.exp.df, aes(x=module28Ensembl, y=expression)) + 
  geom_boxplot(aes(fill=bio.rep)) +
  scale_fill_discrete(name= "Biological\nReplicate",
                      labels= c("2_(high MEF2C)", "7_(low MEF2C)")) +
  labs(title= "Allen derived A9 marker expression") +
  ylab("Expression (normalized FPKM)") +
  xlab("Gene") +
  theme_grey(base_size = 18)
ggsave(file= "A9 marker gene expression.pdf")
```

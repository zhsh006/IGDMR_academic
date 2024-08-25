imt <- read.csv("./imprint_target.csv")
getwd()

head(imt)
### SCore即bindinggp
imt2 <- imt[imt$bindingp >= 0.9, ]
nrow(imt2)
nrow(imt)
head(imt2)
imt_targetscan <- imt[imt$TargetScan == 1, ]
imt_mirdb <- imt[imt$miRDB == 1, ]
nrow(imt_mirdb)
nrow(imt_targetscan)

# gene=imt_mirdb$genesymbol
length(gene)
gene <- imt2$genesymbol
gene2 <- unique(gene)
length(gene2)
length(gene2)
library(clusterProfiler)
gene3 <- bitr(gene2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Rn.eg.db")
kegg <- enrichKEGG(
    gene = gene3$ENTREZID,
    organism = "rno",
    pvalueCutoff = 0.5
)
head(kegg, 10)
load(file = "./mRNA成年/GSEA-mRNA/kegg.rat.sig.Rdata")
load(file = "./mRNA成年/GSEA-mRNA/kegg.rat.sigmet.Rdata")
getwd()

###########
kegg.rat.sig.df <- data.frame(
    TERM = rep(
        names(kegg.rat.sig),
        sapply(kegg.rat.sig, length)
    ),
    GENE = unlist(kegg.rat.sig)
)
head(kegg.rat.sig.df)
#####
kegg.rat.sigmet.df <- data.frame(
    TERM = rep(
        names(kegg.rat.sigmet),
        sapply(kegg.rat.sigmet, length)
    ),
    GENE = unlist(kegg.rat.sigmet)
)
head(kegg.rat.sigmet.df)

length(gene3)
kegg2 <- enricher(
    gene = gene3$ENTREZID,
    pvalueCutoff = 0.5,
    pAdjustMethod = "fdr",
    TERM2GENE = kegg.rat.sigmet.df
)
nrow(kegg2)
head(kegg2, 10)
head(kegg2, 30)
keggread <- setReadable(kegg2, OrgDb = "org.Rn.eg.db", keyType = "ENTREZID")
head(keggread, 10)
library(enrichplot)
library(circlize)
library(DOSE)
library(ggplot2)
kegg_dot <- dotplot(kegg3) # nolint: infix_spaces_linter.

library(clusterProfiler.dplyr)
kegg3 <- kegg2 %>% mutate(FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
kegg_dot <- dotplot(kegg3, showCategory = 10, x = "FoldEnrichment")
dot <- kegg_dot +
    scale_color_continuous(low = "#E64B35FF", high = "#91D1C2FF") +
    scale_size(range = c(2, 10)) +
    xlim(1.2, 1.45)
dot
head(kegg2)
keggread[4, ]



head(imt)

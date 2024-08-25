# library(BiocManager)
# library(devtools)
# devtools::install_github("datapplab/pathview")
# devtools::install_local(pkg = "pathview_1.38.0.tgz")
# install.packages("pathview_1.38.0.tgz", repos = NULL, type="source") 
# install.packages("KEGGgraph")
# BiocManager::install("KEGGgraph")
# install.packages("KEGGgraph_1.58.3.tgz", repos = NULL, type="source") 


####加载包————————————————————————————————————————————————————————————————————
library(edgeR)
library(pathview)
library(gage)
library(clusterProfiler)
library(enrichplot)
library(GseaVis)
library(dplyr)
####读入数据————————————————————————————————————————————————————————————————————
library(openxlsx)
dge <- read.xlsx(
    xlsxFile = "./diff_gene/all_dge_edit.xlsx",
    sheet = 1,
    startRow = 1,
    colNames = T,
    check.names = F
)
rawcounts=dge[,c(3:8)]
rownames(rawcounts)=dge$gene_id
######——————————————————————————————————————————————————————————————————————
######差异分析——————————————————————————————————————————————————————————————————————
group <- c(rep("Control", 3), rep("PPCE", 3))
design <- model.matrix(~group)

dgelist <- DGEList(counts = rawcounts, group = group)
y <- calcNormFactors(dgelist, method = "TMM")
nrow(y)
keep <- rowSums(cpm(y)>1) >= 1 #########筛选基因CPM大于1的
y=y[keep,,keep.lib.sizes=FALSE]

# #####TMM标准化表达量##——————————————##
# head(y)
# tmm <- t(t(y$counts) / y$samples$lib.size / y$samples$norm.factors) * 1000000 ## 获取标准化矩阵
# head(tmm)
# tmm["rno-miR-770-3p", ]
# write.csv(tmm, "TMM-nomalized.csv")
# tmm2 <- as.data.frame(tmm)
# head(tmm2)
# ##################————————————————————

y2 <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y2, design, robust = TRUE)
lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist$counts))
lrt2 <- lrt$table
head(lrt2)
colnames(lrt2)[6]=c("ENSEMBL")
lrt2$ENSEMBL=rownames(lrt2)
bt=bitr(lrt2$ENSEMBL,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Rn.eg.db",drop = TRUE)
exp_dif=inner_join(lrt2,bt,by = "ENSEMBL")
# bt2=bitr(lrt2$ENSEMBL,fromType="ENSEMBL",toType=c("ENTREZID","SYMBOL"),OrgDb="org.Rn.eg.db",drop = TRUE)
# exp_dif2=left_join(lrt2,bt2,by = "ENSEMBL")
# write.csv(exp_dif2,"exp.csv")

###################————————————————————————————————————————————————————————————————————

##########富集分析GSEA和ORA，GSEA先做genelist，ORA直接筛选logFC的差异基因进行富集
exp_dif_order=exp_dif[order(exp_dif$logFC,decreasing = T),]
head(exp_dif_order)
genelist=exp_dif_order$logFC
names(genelist)=as.character(exp_dif_order$ENTREZID)

###pathview查看一下通路中的基因上下调情况
pathview(gene.data = genelist,pathway.id = "rno04150",species = "rno")
#######

########GSEA分析，用特定的Geneset
library(gage)
####*********此处从gage包中提取特定的基因集。此基因集只有KEGG的信号通路，不包含代谢、疾病等通路
kegg.rat <- kegg.gsets(species = "rno", id.type = "kegg")
save(kegg.rat, file = "kegg.rat.all.Rdata")
head(kegg.rat)
kegg.rat.sig <- kegg.rat$kg.sets[kegg.rat$sig.idx]
kegg.rat.met <- kegg.rat$kg.sets[kegg.rat$met.idx]
kegg.rat.sigmet <- kegg.rat$kg.sets[kegg.rat$sigmet.idx]
kegg.rat.dise <- kegg.rat$kg.sets[kegg.rat$dise.idx]
save(kegg.rat.sig, file = "kegg.rat.sig.Rdata")
save(kegg.rat.met, file = "kegg.rat.met.Rdata")
save(kegg.rat.dise, file = "kegg.rat.dise.Rdata")
class(kegg.rat)
##########将特地基因集转换为dataframe
kegg.rat.sigmet.df <- data.frame(
    TERM = rep(
        names(kegg.rat.sigmet),
        sapply(kegg.rat.sigmet, length)
    ),
    GENE = unlist(kegg.rat.sigmet)
)
head(kegg.rat.sigmet.df)
###产看下PI3k-AKT信号通路
library(tidyverse)
kegg.rat.sigmet.df[str_detect(kegg.rat.sigmet.df$TERM,"PI3K"),]
kegg.rat.sigmet.df[str_detect(kegg.rat.sigmet.df$TERM,"mTOR"),]
############用特定数据集进行富集分析GSEA分析

clp_gsea=clusterProfiler::GSEA(
    geneList = genelist,
    minGSSize = 5,
    maxGSSize = 500,
    pvalueCutoff = 1,
    TERM2GENE = kegg.rat.sigmet.df
)

gseaNb(clp_gsea,geneSetID = "rno04151 PI3K-Akt signaling pathway",addPval = T)
gseaNb(clp_gsea,geneSetID = "rno04150 mTOR signaling pathway",addPval = T)

######################用特定数据集进行ORA富集分析
gene=names(genelist)[abs(genelist)>0.2] 
clp_ora=clusterProfiler::enricher(
    gene = gene,
    pvalueCutoff = 1,
    TERM2GENE = kegg.rat.sigmet.df
)
head(clp_ora)
barplot(clp_ora,showCategory=20)
####用在线KEGG数据集进行ORA富集
ekegg=enrichKEGG(
    gene=gene,
    organism = "rno",
    keyType = "kegg",
    pvalueCutoff = 0.5
)
head(ekegg)
dotplot(ekegg)

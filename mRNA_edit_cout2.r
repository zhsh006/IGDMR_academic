######## 找出Protein digestion and absorption和Fc gamma R-mediated 和B cell
getwd()
rm(list = ls())
library(openxlsx)
et <- read.xlsx("/Users/zhsh006/Desktop/父系IGDMR文章/测序/mRNA成年/edit_count/exp_join3.xlsx")
head(et)
rownames(et) <- et$SYMBOL
# rownames(et)=et$ENTREZID
nrow(et)



###################
library(edgeR)
group <- c(rep("Control", 3), rep("PPCE", 3))
design <- model.matrix(~group)
head(design)

dgelist <- DGEList(counts = et[, 8:13], group = group)
y <- calcNormFactors(dgelist, method = "TMM")
nrow(y)
keep <- rowSums(cpm(y) > 1) >= 1
y <- y[keep, , keep.lib.sizes = FALSE]

y2 <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y2, design, robust = TRUE)
lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist$counts))
head(lrt)
lrt2 <- lrt$table
head(lrt2)
#######################





library(dplyr)
library(tidyverse)
library(clusterProfiler)
lrt2$SYMBOL <- rownames(lrt2)
head(lrt2)
tg <- bitr(lrt2$SYMBOL, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Rn.eg.db")
head(tg)
exp_dif <- right_join(tg, lrt2, by = "SYMBOL")
head(exp_dif)
#############
rownames(exp_dif)
library(EnhancedVolcano)
lab_italics <- exp_dif$SYMBOL
selectLab_italics <- ora_sig$SYMBOL
p <- EnhancedVolcano(
    exp_dif,
    x = "logFC",
    y = "PValue",
    lab = lab_italics, # 基因名
    selectLab = selectLab_italics, #s筛选基因名
    pCutoff = 0.05, # p值截断值：水平线
    FCcutoff = 1, # FC截断值：垂直线
    cutoffLineWidth = 0.7, # 虚线粗细
    cutoffLineType = "twodash", # 线的类型
    xlim = c(-5, 5), # x轴起始
    ylim = c(0, -log10(10e-5)), # y轴起始
    pointSize = 3.8, # 点的大小
    labSize = 4, # 标签大小
    xlab = bquote(~ Log[2] ~ "fold change"), # 此行为默认x轴名
    ylab = bquote(~ -Log[10] ~ italic(p - value)), # 修改y轴名为斜体
    axisLabSize = 17, # 坐标轴字体大小
    title = "mRNA differential expression", # 修改标题
    titleLabSize = 18, # 标题大小
    subtitle = bquote(italic("Volcano plot")), # 修改子标题：斜体
    subtitleLabSize = 17, # 子标题大小
    legendLabSize = 13, # legend大小
    col = c("grey20", "darkgreen", "royalblue", "red3"), # legend颜色
    colAlpha = 0.7, # 透明度
    gridlines.major = FALSE, # 背景网格
    gridlines.minor = FALSE,
    drawConnectors = T,
    widthConnectors = 0.4,
    colConnectors = "#7E6148FF"
)
p
###############
head(exp_order)
nrow(exp_order)
head(et)
head(et[, c(1, 8:13)])
exp_join <- left_join(exp_dif, et[, c(1, 8:13)], by = "SYMBOL")
head(exp_join)
getwd()
write.xlsx(exp_join, "./exp_join3.xlsx")


exp_order <- exp_join[order(exp_join$logFC, decreasing = T), ]
head(exp_order)
exp_order[exp_order$SYMBOL == "Akt2", ]



genelt <- exp_order$logFC
names(genelt) <- as.character(exp_order$ENTREZID)
head(genelt)


downgen <- genelt[genelt < 0]
names(downgen)

upgen <- genelt[genelt > 0]
names(upgen)
updown <- data.frame(down = names(downgen[1:6000]), up = names(upgen[1:6000]), fix.empty.names = T)

head(updown)
bt1 <- bitr(updown$down, fromType = "ENTREZID", toType = "SYMBOL", org.Rn.eg.db)
head(bt1)
bt2 <- bitr(updown$up, fromType = "ENTREZID", toType = "SYMBOL", org.Rn.eg.db)
head(bt2)
write.xlsx(bt1, "./edit_count/down.xlsx")
write.xlsx(bt2, "./edit_count/up.xlsx")


##################


load(file = "../GSEA-mRNA/kegg.rat.sig.Rdata")
load(file = "../GSEA-mRNA/kegg.rat.sigmet.Rdata")
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
#####
head(lrt2)

write.csv(setReadable(gseK, org.Rn.eg.db, keyType = "ENTREZID"), "readablegsek.csv")
head(genelist)
gse_auto <- gseKEGG(
    gene = genelt,
    organism = "rno",
    pvalueCutoff = 1,
    nPermSimple = 10000
)



dotplot(gse_auto)
library(enrichplot)
upsetplot(gse_auto)
gseaplot2(gse_auto, geneSetID = 1)


genelist2 <- sort(genelist)

head(genelist2)

dotplot(gse_auto, showCategory = 15)
upsetplot(gseK2)

gseK <- clusterProfiler::GSEA(
    geneList = genelt,
    pvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    TERM2GENE = kegg.rat.sig.df,
    by = "DOSE",
    nPerm = 500
)
upsetplot(gseK2, showCategory = 20)
dotplot(gseK2, showCategory = 20)

head(gseK, 10)
dotplot(gseK, showCategory = 15)
barplot(gseK)
library(GseaVis)
gseaNb(gseK, geneSetID = "rno04151 PI3K-Akt signaling pathway", addPval = T)
gseaNb(gseK, geneSetID = "rno04150 mTOR signaling pathway", addPval = T)

#######################


gseK2 <- clusterProfiler::GSEA(
    geneList = genelt,
    pvalueCutoff = 0.5,
    pAdjustMethod = "fdr",
    TERM2GENE = kegg.rat.sigmet.df,
    by = "fgsea",
    nPerm = 500
)

write.csv(gseK, "./edit_count/gseK.csv")

head(gseK2)
dotplot(gseK2, showCategory = 10)
clusterProfiler::cnetplot(gseK2, foldChange = genelist)

gseK2.2 <- pairwise_termsim(gseK2)
p1 <- treeplot(gseK2.2)
p1
p2 <- treeplot(gseK2.2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels = "A")

dotplot(gseK2, showCategory = 20)
write.xlsx(gesK, )
head(kegg.rat.sig.df)

rawgene <- kegg.rat.sig.df[kegg.rat.sig.df$TERM == "rno04151 PI3K-Akt signaling pathway", ]
head(rawgene)
rawimp <- read.csv("PI3k-target-imprint.csv")
head(rawimp)
rawimp2 <- na.omit(rawimp)
head(rawimp2)
head(exp_join)

exp_join <- inner_join(exp_order, rawimp2, by = "ENTREZID")
rawimp2$ENTREZID <- as.character(rawimp2$ENTREZID)
head(exp_join)
write.xlsx(exp_join, "exp_join_imp_mir.xlsx")


head(exp_order)
exp_order$ENTREZID <- rownames(exp_order)
gsek <- setReadable(gseK, OrgDb = "org.Rn.eg.db", keyType = "ENTREZID")
write.xlsx(gsek, "today1.xlsx")
dotplot(gsek, showCategory = 15)

ex <- as.data.frame(gseK2)
head(ex)
gseup <- ex[ex$NES > 0, ]
gsedown <- ex[ex$NES < 0, ]
nrow(gseup)
nrow(gsedown)
head(gsedown, 15)
head(gseup, 15)
browseKEGG(gseup, "rno04150")


# library("pathview")
# hsa04110 <- pathview(
#     gene.data = geneList,
#     pathway.id = "hsa04110",
#     species = "hsa",
#     limit = list(gene = max(abs(geneList)), cpd = 1)
# )
library(enrichplot)
install.packages("ggnewscale")
gseK3 <- setReadable(gseK2, "org.Rn.eg.db", "ENTREZID")
cnetplot(gseK3, foldChange = genelist, circular = TRUE, colorEdge = TRUE)
dotplot(gseK2)

library(clusterProfiler)
install.packages("ggupset")
upsetplot(gseK2, showCategory = 20)
head(gseK2)
write.csv(gseK2,"gseK2.csv")

sort(gseK2,by="Pvalue")
write.csv(setReadable(gseK2,"org.Rn.eg.db",keyType = "ENTREZID"),"gsek22.csv")
getwd()
head(gseK2,10)
write.csv(gseK2,"gseK2.csv")


write.csv(setReadable(gseK2_up,OrgDb = "org.Rn.eg.db",keyType = "ENTREZID"),"gseK2_up.read.csv")
write.csv(setReadable(gseK2_down,OrgDb = "org.Rn.eg.db",keyType = "ENTREZID"),"gseK2_down.read.csv")


read=setReadable(gseK2_down,OrgDb = "org.Rn.eg.db",keyType = "ENTREZID")
head(read)

read$ID 
selecet_gene_pi3k = read %>% filter(ID == "rno04151 PI3K-Akt signaling pathway")
selecet_gene_pi3k@result[[11]]
pi3kgene_core=strsplit(selecet_gene_pi3k@result[[11]],split = "/")
selecet_gene_mtor = read %>% filter(ID == "rno04150 mTOR signaling pathway")
mtorgene_core=strsplit(selecet_gene_mtor@result[[11]],split = "/")


library(tidyverse)
mtorgene_core2=as.character(mtorgene_core[[1]])
pi3kgene_core2=as.character(pi3kgene_core[[1]])
core_gene=c(mtorgene_core2,pi3kgene_core2)
core_gene
length(core_gene)
for_fi=data.frame(SYMBOL=core_gene,set=1:length(core_gene))
head(for_fi)
head(exp_dif)
path_core_gene=merge(exp_dif,for_fi,by="SYMBOL")
path_core_gene
write.csv(path_core_gene,"path_core_gene.csv")
##################heatmap
path_core_gene2=path_core_gene[!duplicated(path_core_gene$SYMBOL),]

nrow(path_core_gene2)
rownames(path_core_gene2)=path_core_gene2$SYMBOL
head(path_core_gene2)
library(ComplexHeatmap)
df=path_core_gene2[,2:7]
head(df)
class(t(scale(t(df))))

??Heatmap
col_fun=colorRamp2(c(-2, 0, 2), c('#E15A29','#AEC7E8FF','#009844'))
Heatmap(t(scale(t(df))),
        cluster_rows = F,
        width=ncol(path_core_gene2)*unit(7, "mm"),
        height=nrow(path_core_gene2)*unit(13, "mm"),
        #column_km = 1,
        cluster_columns = T,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson",
        clustering_method_rows = "complete",
        column_order = sort(colnames(path_core_gene2[,2:7])),
        column_dend_reorder = T,
        column_split = factor(rep(c('Control', 'PPCE'), each=3)),
        right_annotation = NULL,
        col = col_fun,
        border = T,
        border_gp = gpar(col="black")
        )
miR_exp_forannotion[91:92,]
shedingx=miR_exp_forannotion$logFC
shedingx[1:92]
library(RColorBrewer)
library(circlize)
col_fun=colorRamp2(c(-2, 0, 2), c('#009844','#AEC7E8FF','#E15A29'))
ha2=rowAnnotation(logFC=shedingx[1:92],col = list(logFC=col_fun))
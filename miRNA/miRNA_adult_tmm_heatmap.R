head(hip_cmvskh)
head(miRNA_adult_edit)
colnames(hip_cmvskh)
colnames(成年edit)=colnames(hip_cmvskh)


library(tidyverse)
exp_today1 = full_join(成年edit,hip_cmvskh,by="sRNA.readcount")
head(exp_today1)
write.csv(exp_today1,"all_adult_count.csv")
nrow(hip_cmvskh)
nrow(miRNA_adult_edit)

tmm2=log(tmm)

head(tmm2)
tmm2[!is.infinite(tmm2$hip_KH_1.readcount.y),]
colnames(tmm2)=c("c1","c2","c3","p1","p2","p3")
tmm3[!is.infinite(tmm3),]
class(tmm2)
tmm3=as.tibble(tmm2)
head(tmm3)
rownames(tmm3)
tmm3[is.finite(tmm3[,1:6]),]
remove(tmm4)
tmm3[sapply(tmm3, is.infinite),]
rowSums(tmm3)

tmm3[mapply(is.infinite, tmm3)]


tmm4=tmm3[is.finite(rowSums(tmm3)),]
head(tmm4)

nrow(tmm4)



library(edgeR)

group=c(rep("Control",3),rep("PPCE",3))
design=model.matrix(~group)
head(design)
head(all_adult_count)
co1=all_adult_count[1:99,]
tail(co1)

dgelist=DGEList(counts = co1[,2:7],group = group)
y<-calcNormFactors(dgelist,method = 'TMM')
head(y)
tmm<-t(t(y$counts)/y$samples$lib.size/y$samples$norm.factors)*1000000 ##获取标准化矩阵
head(tmm)

rownames(tmm)=co1$miRNA_ID
tmm["rno-miR-770-3p",]
write.csv(tmm,"TMM-nomalized.csv")

tmm_today=tmm
tmm_today=log(tmm_today)
tmm_today2=tmm_today[is.finite(rowSums(tmm_today)),]
head(tmm_today2)


head(miR_exp)
Heatmap(scale_imprint_count_tmm,
        cluster_rows = T,
        #width=ncol(all)*unit(7, "mm"),
        #height=nrow(all)*unit(3.5, "mm"),
        #column_km = 1,
        cluster_columns = T,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson",
        clustering_method_rows = "complete",
        column_order = sort(colnames(scale_imprint_count_tmm)),
        column_dend_reorder = T,
        column_split = factor(rep(c('A', 'B'), each=3)),
        right_annotation = ha2
        )

pheatmap(log(tmm_today2),
         color=colorRampPalette(c('green','white','red'))(100),
         fontsize_row=9,
         scale='col',
         cluster_rows = T,
         cluster_cols = T,
         cellheight = 12,
         cellwidth = 20,
         display_numbers = T,
         border_color='black')


xxx=scale(log(tmm_today2))
xxx["rno-miR-377-5p",]

exp <- apply(tmm_today2,1,scale)
head(exp)
rownames(exp)=colnames(tmm_today2)
exp=t(exp)
exp["rno-miR-541-5p",]
tmm_today2["rno-miR-541-5p",]

head(all_adult_count)
nrow(all_adult_count)
imprint_edit <- read_excel("imprint_edit.xlsx")
head(imprint_edit)
imprint_edit_today = full_join(imprint_edit,all_adult_count,by="miRNA_ID")
write.csv(imprint_edit_today,"imprint_edit_today.csv")
head(imprint_edit_today)

raw=read.csv("miRNA_rawcount.csv")
head(raw)

dgelist2=DGEList(counts = raw[,2:7],group = group)
y_today<-calcNormFactors(dgelist2,method = 'TMM')
head(y)
raw_today<-t(t(y_today$counts)/y_today$samples$lib.size/y_today$samples$norm.factors)*1000000 ##获取标准化矩阵
head(raw_today)
rownames(raw_today)=raw$miRNA_ID
write.csv(raw_today,"raw_today_TMM.csv")


nrow(y$counts)


raw_today=log(raw_today)
raw_today2=raw_today[is.finite(rowSums(raw_today)),]
head(raw_today2)
raw_today3 <- apply(raw_today2,1,scale)
head(raw_today3)
rownames(raw_today3)=colnames(raw_today2)
raw_today3=t(raw_today3)
exp["rno-miR-541-5p",]
raw_today2["rno-miR-541-5p",]
nrow(raw_today3)




#####全部counts：raw
#####全部TMM：raw_today
#####印迹counts：imprint_edit_today前99行
#####从全部TMM中，提取imprint_TMM进行作图

head(raw)
head(raw_today)
raw_today27=raw_today
lie=raw$miRNA_ID
raw_today27=as.data.frame(raw_today27)
raw_today27=raw_today27%>%mutate(miRNA_ID=lie,.before = hip_Control_1)
class(raw_today27)
head(raw_today27)
nrow(raw_today27)

all_count_tmm=full_join(raw,raw_today27,by="miRNA_ID")
write.csv(all_count_tmm,"miRNA_adult_all_count_tmm.csv")
imprint_count_tmm=all_count_tmm[1:99,]
write.csv(imprint_count_tmm,"imprint_adult_count_tmm.csv")

rownames(imprint_count_tmm)=imprint_count_tmm[,1]
imprint_count_tmm2=imprint_count_tmm[,-1]
imprint_count_tmm3=imprint_count_tmm2[,7:12]
log_imprint_count_tmm=log(imprint_count_tmm3)
log_imprint_count_tmm=log_imprint_countß_tmm[is.finite(rowSums(log_imprint_count_tmm)),]
head(log_imprint_count_tmm)
scale_imprint_count_tmm <- apply(log_imprint_count_tmm,1,scale)
head(scale_imprint_count_tmm)
rownames(scale_imprint_count_tmm)=colnames(log_imprint_count_tmm)
head(scale_imprint_count_tmm)
scale_imprint_count_tmm=t(scale_imprint_count_tmm)
head(scale_imprint_count_tmm)



dgelist2=DGEList(counts = raw_today27[,2:7],group = group)
y_today<-calcNormFactors(dgelist2,method = 'TMM')
head(y_today)
nrow(y_today$counts)
y2 <- estimateDisp(y_today, design, robust=TRUE)
fit <- glmQLFit(y2, design, robust = TRUE)
lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist2$counts))
head(lrt)

lrt2<-lrt$table
miR_diff <- lrt2[order(lrt2$FDR, lrt2$logFC, decreasing = c(FALSE, TRUE)), ]
head(miR_diff)
nrow(miR_diff)
write.csv(miR_diff,"miR_diff_today_all.csv")
miR_diff2=miR_diff %>% mutate(miRNA_ID=rownames(miR_diff),.before = logFC)

View(miR_diff)
miR_exp=full_join(miR_diff2,all_count_tmm,by="miRNA_ID")
head(miR_exp)
write.csv(miR_exp,"miRNA_adult_exp.csv")

imprint_count_tmm3
imprint_count_tmm4=imprint_count_tmm3 %>% mutate(miRNA_ID=rownames(imprint_count_tmm3),.before = hip_Control_1.y)
scale_imprint_count_tmm2=as.data.frame(scale_imprint_count_tmm)
scale_imprint_count_tmm2=scale_imprint_count_tmm2 %>% mutate(miRNA_ID=rownames(scale_imprint_count_tmm2),.before = hip_Control_1.y)
miR_exp_forannotion=full_join(scale_imprint_count_tmm2,miR_diff2,by="miRNA_ID")
head(miR_exp_forannotion)
ha=rowAnnotation(foo = shedingx[1:92])

nrow(scale_imprint_count_tmm2)
miR_exp_forannotion[91:92,]
shedingx=miR_exp_forannotion$logFC
shedingx[1:92]
col_fun=colorRamp2(c(-2, 0, 2), c("green", "grey", "red"))
ha2=rowAnnotation(logFC=shedingx[1:92],col = list(logFC=col_fun))


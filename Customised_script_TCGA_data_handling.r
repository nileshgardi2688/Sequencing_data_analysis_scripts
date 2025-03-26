#R 4.3

options(download.file.method = 'wininet')
options(download.file.method = "libcurl")
library(tidyverse)
library(palmerpenguins)
library(ggstatsplot)
library(doParallel)
#library(microbenchmark)
library(UCSCXenaTools)
library(survival)
library(survminer)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggridges)
library(enrichplot)

#my_path<-setwd("D:/Backup_90046043/Work/TCGA/xena")

df_exp=read.table("L_H_data.txt",header=T)
df_exp<-as_tibble(df_exp)

df_exp_new =
df_exp %>%
mutate_if(is.numeric, function(x){round(2^x-1, digits = 0)})

#df_exp[1:5, 1:6]
#df_exp_new[1:5, 1:6]

# check if missing value exist and convert dataframe into matrix format
df_exp_new %>%
select_if(is.numeric) %>%
as.matrix %>% {.->>m_exp} %>%
as.vector %>%
is.na %>%
sum

row.names(m_exp) = df_exp_new$sample
m_exp[1:6, 1:5]

# identify rows with all zero
all_zero = map_lgl(1:nrow(m_exp), function(i){ sum(m_exp[i,] == 0) == ncol(m_exp)})
table(all_zero)
all_zero

# retain rows not all zero
m_exp_new = m_exp[-which(all_zero == TRUE),]
dim(m_exp)
dim(m_exp_new)

# prep sample annotation
sample_class = map_chr(colnames(m_exp_new), function(x){
# 14th digits of sample ID: tumor ('0'); adjacent normal ('1')
ifelse(str_sub(x, 1,1)=='L', yes = 'Low_IRAK3', no = 'High_IRAK3')
}) %>% factor()

table(sample_class)

# create a *DGElist* object from a table of counts and annotation
d0 = DGEList(counts = m_exp_new, group = sample_class)
# prep normalization factor
d0 = calcNormFactors(d0)

# remove lowly expressed genes (try different cutoff by your own based on the two figures below)
cutoff = 100 # expression cutoff
cut = which(apply(cpm(d0), 1, max) < cutoff)
d1 = d0[-cut,]
dim(d0)
dim(d1)

design = model.matrix(~sample_class)

d2 = voom(d1, design, plot = T, normalize = 'quantile')

# linearly fit to identify DEGs with corresponding significance
fit = lmFit(d2,design)
fit = eBayes(fit,robust=TRUE,trend=TRUE)
top.table = topTable(fit, n=Inf) %>% na.omit()

write.table(top.table,"top_table.txt",sep="\t")

# arrange DEG table by log Fold change (logFC) and attach Entrez ID column
DEG_table =
top.table %>%
rownames_to_column('gene') %>%
as_tibble() %>%
dplyr::select(gene, everything()) %>%
arrange(desc(logFC))

DEG_symbol2entrez = bitr(DEG_table$gene, fromType = 'SYMBOL', toType =
'ENTREZID', OrgDb = 'org.Hs.eg.db', drop = TRUE) %>% as_tibble()
DEG_table_new = inner_join(DEG_table, DEG_symbol2entrez, by =
c('gene'='SYMBOL'))

# DEG gene list for ORA
DEG_table_up = filter(DEG_table_new, adj.P.Val <0.05, logFC >= 1.5) %>%
arrange(desc(logFC))
DEG_table_dn = filter(DEG_table_new, adj.P.Val <0.05, logFC <= -1.5) %>%
arrange(desc(logFC))
geneList_ora_up = DEG_table_up$gene
geneList_ora_dn = DEG_table_dn$gene

#enrichment analysis: over-representation (ORA) and geneset enrichment (GSEA􀔼
#ORA Biological process (BP)
ora_up_go_BP = enrichGO(
gene = geneList_ora_up,
OrgDb = org.Hs.eg.db, # species: human
keyType="SYMBOL", # type of ORA input
ont = "BP", # biological process from GO
pAdjustMethod = "BH")

ora_dn_go_BP = enrichGO(
gene = geneList_ora_dn,
OrgDb = org.Hs.eg.db, # species: human
keyType = "SYMBOL", # type of ORA input
ont = "BP", # biological process from GO
pAdjustMethod = "BH")

ora_up_go_BP@result %>% as_tibble() %>% head()
ora_dn_go_BP@result %>% as_tibble() %>% head()
write.table(ora_up_go_BP,"ORA_UP_BP.txt",sep="\t")
write.table(ora_dn_go_BP,"ORA_DN_BP.txt",sep="\t")

barplot(ora_up_go_BP, showCategory=15,font.size=7) + ggtitle("barplot for GO BP(upregulated DEG)")
barplot(ora_dn_go_BP, showCategory=15,font.size=7) + ggtitle("barplot for GO BP(downregulated DEG)")
dotplot(ora_up_go_BP, showCategory=15,font.size=7) + ggtitle("dotplot for GO BP(upregulated DEG)")
dotplot(ora_dn_go_BP, showCategory=15,font.size=7) + ggtitle("dotplot for GO BP(downregulated DEG)")

#ORA Molecular Function (MF)
ora_up_go_MF = enrichGO(
gene = geneList_ora_up,
OrgDb = org.Hs.eg.db, # species: human
keyType="SYMBOL", # type of ORA input
ont = "MF", # biological process from GO
pAdjustMethod = "BH")

ora_dn_go_MF = enrichGO(
gene = geneList_ora_dn,
OrgDb = org.Hs.eg.db, # species: human
keyType = "SYMBOL", # type of ORA input
ont = "MF", # biological process from GO
pAdjustMethod = "BH")

ora_up_go_MF@result %>% as_tibble() %>% head()
ora_dn_go_MF@result %>% as_tibble() %>% head()
write.table(ora_up_go_MF,"ORA_UP_MF.txt",sep="\t")
write.table(ora_dn_go_MF,"ORA_DN_MF.txt",sep="\t")

barplot(ora_up_go_MF, showCategory=15,font.size=7) + ggtitle("barplot for GO MF(upregulated DEG)")
barplot(ora_dn_go_MF, showCategory=15,font.size=7) + ggtitle("barplot for GO MF(downregulated DEG)")
dotplot(ora_up_go_MF, showCategory=15,font.size=7) + ggtitle("dotplot for GO MF(upregulated DEG)")
dotplot(ora_dn_go_MF, showCategory=15,font.size=7) + ggtitle("dotplot for GO MF(downregulated DEG)")

#ORA Cellular Component (CC)
ora_up_go_CC = enrichGO(
gene = geneList_ora_up,
OrgDb = org.Hs.eg.db, # species: human
keyType="SYMBOL", # type of ORA input
ont = "CC", # biological process from GO
pAdjustMethod = "BH")

ora_dn_go_CC = enrichGO(
gene = geneList_ora_dn,
OrgDb = org.Hs.eg.db, # species: human
keyType = "SYMBOL", # type of ORA input
ont = "CC", # biological process from GO
pAdjustMethod = "BH")

ora_up_go_CC@result %>% as_tibble() %>% head()
ora_dn_go_CC@result %>% as_tibble() %>% head()
write.table(ora_up_go_CC,"ORA_UP_CC.txt",sep="\t")
write.table(ora_dn_go_CC,"ORA_DN_CC.txt",sep="\t")

barplot(ora_up_go_CC, showCategory=15,font.size=7) + ggtitle("barplot for GO CC(upregulated DEG)")
barplot(ora_dn_go_CC, showCategory=15,font.size=7) + ggtitle("barplot for GO CC(downregulated DEG)")
dotplot(ora_up_go_CC, showCategory=15,font.size=7) + ggtitle("dotplot for GO CC(upregulated DEG)")
dotplot(ora_dn_go_CC, showCategory=15,font.size=7) + ggtitle("dotplot for GO CC(downregulated DEG)")

# DEG gene list for GSEA
#geneList_gsea = DEG_table_new$logFC %>% setNames(DEG_table_new$ENTREZID)
geneList_gsea_up = DEG_table_up$logFC %>% setNames(DEG_table_up$ENTREZID)
geneList_gsea_dn = DEG_table_dn$logFC %>% setNames(DEG_table_dn$ENTREZID)
geneList_gsea = c(geneList_gsea_up, geneList_gsea_dn)
#GSEA
gsea_kegg = gseKEGG(
geneList_gsea,
organism = "hsa", # species: human
keyType = "ncbi-geneid", # type of GSEA input; symbol not accepted
minGSSize = 10,
maxGSSize = 500,
eps = 1e-10,
pvalueCutoff = 0.05,
pAdjustMethod = "BH",
verbose = TRUE,
use_internal_data = FALSE
)


gsea_kegg@result %>% as_tibble() %>% arrange(desc(NES)) %>% head()
write.table(gsea_kegg,"gsea_kegg.txt",sep="\t")

########################################Visualisation
#concept network plot
enrichplot::cnetplot(gsea_kegg, foldChange = geneList_gsea) +
scale_color_viridis_c() +
labs(color = 'Fold\nChange', title = 'Concept network') +
theme_void() +
theme(plot.title = element_text(face = 'bold', size = 24))

#pairwise_termsim
pairwise_termsim(gsea_kegg) %>% emapplot(layout="kk")

#ridge plot
ridgeplot(gsea_kegg,label_format = 40,showCategory = 20,decreasing=T) + labs(x = 'NES', title= 'ridge plot for GSEA') +
    theme_classic()


# file_gmt_symbol ="C:/Users/90046043/Downloads/c2.cp.v2023.2.Hs.symbols.gmt"
# file_gmt_entrez ="C:/Users/90046043/Downloads/c2.cp.v2023.2.Hs.entrez.gmt"
# gs_hallmark_symbol = read.gmt(file_gmt_symbol)
# gs_hallmark_entrez = read.gmt(file_gmt_entrez)
# head(gs_hallmark_symbol)

# ora_up_hallmark = enricher(geneList_ora_up, TERM2GENE = gs_hallmark_symbol)
# ora_dn_hallmark = enricher(geneList_ora_dn, TERM2GENE = gs_hallmark_symbol)

# barplot(ora_up_hallmark, showCategory=20) + ggtitle("barplot for hallmark (up-regulated DEG)")

# gsea_hallmark = GSEA(geneList_gsea, TERM2GENE = gs_hallmark_entrez)
# gsea_hallmark
# idx_2 = grep(gsea_hallmark@result$Description, pattern =
# 'HALLMARK_G2M_CHECKPOINT', ignore.case = T)
# gseaplot2(gsea_hallmark, geneSetID = idx_2, title =
# gsea_hallmark@result$Description[idx_2])


###################################Other visualisations of GSEA results
# library(DOSE)
# geneList_gsea_up = DEG_table_up$logFC %>% setNames(DEG_table_up$ENTREZID)
# geneList_gsea_dn = DEG_table_dn$logFC %>% setNames(DEG_table_dn$ENTREZID)
# geneList_gsea = c(geneList_gsea_up, geneList_gsea_dn)

# de <- names(geneList_gsea)
# edo <- enrichDGN(de)

# edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
# p1 <- cnetplot(edox, color.params = list(foldChange = geneList_gsea))
# p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList_gsea)
# p3 <- cnetplot(edox, foldChange=geneList_gsea, circular = TRUE, colorEdge = TRUE) 
# cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

#########################
de <- names(geneList_gsea)
gsea_kegg = gseKEGG(
geneList_gsea,
organism = "hsa", # species: human
keyType = "ncbi-geneid", # type of GSEA input; symbol not accepted
minGSSize = 30,
maxGSSize = 500,
eps = 1e-10,
pvalueCutoff = 0.05,
pAdjustMethod = "BH",
verbose = TRUE,
use_internal_data = FALSE
)

edox <- setReadable(gsea_kegg, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, color.params = list(foldChange = geneList_gsea),showCategory=10,node_label="category")
p2 <- cnetplot(edox, categorySize="pvalue", color.params = list(foldChange = geneList_gsea),showCategory=10,node_label="category")
p3 <- cnetplot(edox, color.params = list(foldChange = geneList_gsea), circular = TRUE, colorEdge = TRUE,showCategory=10,node_label="category") 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))


#p1 <- cnetplot(edox, node_label="category",showCategory=10, 
        cex_label_category = 1.2) 
#p2 <- cnetplot(edox, node_label="gene", showCategory=10,
        cex_label_gene = 0.8) 
p3 <- cnetplot(edox, node_label="all",showCategory=10,cex_label_gene = 0.5,cex_label_category = 1,color_category='firebrick',color.params = list(foldChange = geneList_gsea))
	
#p4 <- cnetplot(edox, node_label="none", showCategory=10,
        color_category='firebrick', 
        color_gene='steelblue') 
#cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


#p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList_gsea, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

emapplot(edox)

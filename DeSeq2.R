main_path<-setwd("C:/Users/Nilesh/Desktop/Salmon study/Progesterone/UCSC_knownGenes/")
samples <- read.table(file.path(main_path, "Sample.txt"), header = TRUE)
files <- file.path(main_path, samples$Patient,samples$run, "quant.sf")
names(files)<-paste(samples$Patient,samples$run,sep="_")
all(file.exists(files))
tx2gene <-read.table("C:/Users/Nilesh/Desktop/Salmon study/UCSCID_to_geneID.txt")
library(tximport)
library(readr)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene,dropInfReps = TRUE)
#countData<-txi$counts[apply(txi$counts==0,1,sum)<=(20*ncol(txi$counts)/100),]
#write.table(countData,"countData_nozero.txt",sep="\t")
##################################################################DEseq2 started 
library("DESeq2")
dds <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ run)
dds<-dds[apply(counts(dds)==0,1,sum)<=(20*ncol(counts(dds))/100),]
dds$run <- factor(dds$run, levels=c("A","B","C"))
dds <- DESeq(dds)
#res <- results(dds)
#res
res1<-results(dds, contrast=c("run","B","A"))
write.table(as.data.frame(res1), file="B_vs_A.txt",sep="\t")
res2<-results(dds, contrast=c("run","C","A"))
write.table(as.data.frame(res2), file="C_vs_A.txt",sep="\t")
res3<-results(dds, contrast=c("run","C","B"))
write.table(as.data.frame(res3), file="C_vs_B.txt",sep="\t")

#####################SVASEQ 
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ run, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))

n.sv <- num.sv(dat,mod, method = 'leek')
n.sv

svseq <- svaseq(dat, mod, mod0, n.sv = 3)

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
ddssva$SV4 <- svseq$sv[,4]
design(ddssva) <- ~ SV1 + SV2 + SV3 + SV4+ run

ddssva <- DESeq(ddssva) 
resSV <- results(ddssva)
#resOrdered <- res[order(res$padj),]
#summary(res)
#sum(res$padj < 0.1, na.rm=TRUE)
#res05 <- results(dds, alpha=0.05)############alpha correspong to adjusted p value cutoff
#summary(res05)

##########################################################Exploring and exporting results
plotMA(res, ylim=c(-2,2))##,alpha=0.1
#idx <- identify(res$baseMean, res$log2FoldChange)############Time consuming step
#rownames(res)[idx]
#plotCounts(dds, gene=which.min(res$padj), intgroup="run")


####################################Individual gene plots
d<-plotCounts(dds, gene="FOSB", intgroup="run",returnData=TRUE)
ggplot(d, aes(x=run, y=count,color = run )) +
    geom_point(position=position_jitter(0.2),shape=16) +
	scale_y_log10()+ geom_point(size = 3) + geom_boxplot() + ggtitle("FOSB")


######################Exporting results to txt files
write.table(as.data.frame(res), file="pre_vs_post.txt",sep="\t")

#############Multi-factor designs
colData(dds)
ddsMF <- dds
design(ddsMF) <- formula(~ Patient + run)
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)
head(resMF)


#####Data transformations and visualization
rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.fast <- vst(dds, blind=FALSE)

write.table(assay(rld),"rld.txt",sep="\t")
write.table(assay(vsd),"vsd.txt",sep="\t")
write.table(assay(vsd.fast),"vsd_fast.txt",sep="\t")
###########################Effects of transformations on the variance
ntd <- normTransform(dds)
library("vsn")
notAllZero <- (rowSums(counts(dds))>0)

png("ntd.png")
meanSdPlot(assay(ntd)[notAllZero,])
dev.off()

png("rld.png")
meanSdPlot(assay(rld[notAllZero,]))
dev.off()

png("vsd.png")
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()
################################Data quality assessment by sample clustering and visualization
library("pheatmap")
#select<-which((res$padj < 0.1 & res$log2FoldChange < -2) | (res$padj < 0.1 & res$log2FoldChange > 2))
select<-which((res1$pvalue < 0.05 & res1$log2FoldChange > 1.4) | (res1$pvalue < 0.05 & res1$log2FoldChange < -1.5))
mat<-assay(rld)[select, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)

#up<-which(res$padj < 0.05 & res$log2FoldChange > 2)
#down<-which(res$padj < 0.05 & res$log2FoldChange < -2)
#select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]################To select any genes for heatmap
#df <- as.data.frame(colData(dds)[,c("condition","type")])

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]################To select any genes for heatmap
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
	 
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)	 
	 
#####################################Heatmap of the sample-to-sample distances	 
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$run, rld$Patient, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(rld, intgroup=c("run", "Patient"))
plotPCA(rld, intgroup=c("run"))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)	


  
pcaData <- plotPCA(rld, intgroup=c("run", "Patient"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=run, shape=run)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
coord_fixed()



#Differential Gene Expression Analysis

#Please, open README file before running and install required libraries

#Load libraries
library("DESeq2")
library(RColorBrewer)
library(gplots)
library(genefilter)
library(ggplot2)
library(EnhancedVolcano)

# Set your working directory
filepath <- "~choose/your/own/directory"
getwd()
setwd(filepath)

#Create output folder
folder <- "output" 
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}

##Paste csv file name of your raw counts
countsName <-"your_raw_counts_file.csv"

#Import counts from csv file
countData <- read.csv(countsName, header = TRUE,row.names=1)
countData <- as.matrix(countData)
head(countData)

#Build metadata matrix (coldata)
#Modify it to match your data
(condition <- factor(c(rep("ctl", 8), rep("exp", 8))))
(coldata <- data.frame(row.names=colnames(countData), condition))

#Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
dds <- DESeq(dds)
dds

# Variance Stabilizing Transformations (VST) corrects for size and normalization factors. 
#The transformed data is on the log2 scale for large counts.
vsd <- vst(dds, blind=FALSE) #takes ~1 second to run but gives a different result
head(assay(vsd))
rld <- rlogTransformation(dds)  #takes ~30 s #alternative form: rlog(dds, blind=FALSE)
head(assay(rld))
#Visualization of rlog transformation
hist(assay(rld))
#Save normalization results 
write.csv(results(dds), file="output/DGE-results.csv")

#Calculate distances between samples
(mycols <- brewer.pal(8, "BuPu")[1:length(unique(condition))])
sampleDists <- as.matrix(dist(t(assay(rld))))
head(sampleDists)

#Create a heatmap image from Sample Distance Matrix
png("output/DGE-heatmap.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "cyan", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()


#Create a PCA plot function
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
#Run and save PCA plot
png("output/DGE-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, intgroup="condition",xlim=c(-75, 35))
dev.off()

#RESULTS
res <- results(dds,
               contrast = c("condition", "exp", "ctl"),
               alpha = 0.1,
               lfcThreshold = 0.32) #Adjust for multiple testing =1.25 foldchange
summary(res)
#Visualize MA plot to observe significant DE genes
plotMA(res,ylim=c(-8,8))

#Shrink results
res <- lfcShrink(dds, contrast = c("condition", "exp", "ctl"), res = res, type="normal")
summary(res)
#Visualize MA plot to observe shrinkage
plotMA(res,ylim=c(-8,8))
head(res)

#Create Data Frame with significant DGE (alpha=0.05)
table(res$padj<0.05)
resdata <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
#Save shrunk results to csv file
write.csv(resdata, file="output/DGE-shrunkresults.csv")
#Quick visualization of significant findings
hist(resdata$pvalue, breaks=50, col="blue")

#VOLCANO PLOT
#Using Enhanced Volcano plot package
png("output/DGE-volcanoplot.png", 1200, 1000, pointsize=20)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'DGE Volcano Plot',
                subtitle= 'Expressed and Unexpressed Samples',
                pCutoff = 10e-6, #to change p-value threshold
                FCcutoff = 1.5, #to change fold change cut-lines
                col=c('gray', 'cyan', 'green', 'red3'),
                shape = 1,
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                cutoffLineWidth = 0.8,
                hlineType = c('solid', 'longdash', 'dotdash', 'dotted'),
                hlineWidth = c(1.0, 1.5, 2.0, 2.5),
                colAlpha = 1,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 8)
dev.off()

#Dispersion plot
png("output/DGE-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()


#Now you have 4 images in your folder with:
#a heatmap, PCA, volcano plot and dispersions plot
#And some csv files with your results

#good luck

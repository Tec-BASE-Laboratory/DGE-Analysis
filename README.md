# DGE-Analysis

This is a Differential Gene Expression (DGE) Analysis with [__DESeq2__](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) using a reference genome.

### __FILE REQUIREMENTS__

- csv file with your raw counts (not-normalized) obtained from a previous Gene Expression Quantification.

- The metadata matrix should contain a column with all your sample names and another column with the corresponding treatment levels (min 2 levels).
Here, metadata matrix is built within the DE_Analysis code, so you'll need to adapt that part of your code to your data.


### __FILE OUTPUTS__

- Heatmap `DGE-heatmap.png`
- Principal Component Analysis `DGE-pca.png`
- Volcano plot `DGE-volcanoplot.png`
- Dispersions plot `DGE-dispersions.png`
- Normalized Results `DGE-results.csv`
- Shrunk Results `DGE-shrunkresults.csv`


### __PACKAGE REQUIREMENTS__

Install required libraries from Bioconductor

*Run this code if the packages are not already installed*

```
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("genefilter")
BiocManager::install("calibrate")
BiocManager::install('MASS')
BiocManager::install('MatrixGenerics')
BiocManager::install('matrixStats')
```

### __Load required libraries__

```
# Differential Expression Analysis
library("DESeq2")
# Heatmap
library(RColorBrewer)
library(gplots)
# PCA plot
library(genefilter)
# Volcano plot
library(ggplot2)
library(EnhancedVolcano)
```

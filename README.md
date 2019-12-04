# FinalProjectDelivery
Differential expression analysis of the genes of male and female patients diagnosed with brain cancer. 

## Downloading the Data: 
To conduct a similar analysis, download gene counts data in the form: htseqcounts
This data can be found for various cancer types on the Cancer Genome Atlas Program Website. Choose various parameters that you would like to include in your analysis (E.g. in this example, the analysis of different genes in males and females are analysed. In my search criteria, the genders were specified. 

https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga -> access tcga data -> gender of choice -> view files in repository -> RNA Seq Data -> HTSeqCounts 

Upon input of this search criteria, the first 20 files were selected for each gender and downloaded. The files were stored in a common and accessible file. The location of storage is important. When retrieving the data in our notebook, this location must be specified. 


## Abstract: 
RNA sequencing data can be used to analyse the differential expression of genes for a particular disease amongst two subsets of patients. In this analysis, the counts data for each gene is used as the input. To this data, pre-processing, exploratory data analysis and differential expression analysis is done. The input data is derived from The Cancer Genome Atlas Program Database in the htseq-counts format. Various packages are used for the computation. EdgeR is used to store the data in a retrievable and organized format. Limma is used to conduct linear modelling and the Glimma package is used to create interactive graph's that allow the end user to explore the data. Additionally, packages such as gplots and heatmap.plus are used to create some useful diagramatic representation of the data. 

## Let's Begin!
## Setting Up: 

The libraries required for analysis are loaded in. Glimma and Limma are used for computation and analysis of differential expression. EdgeR is used to organize the data into a dataframe. Homo.sapiens is a key created by Bioconductor that allows the information of GeneID to be compared to other methods of naming. 

```{r}
library(Glimma)
library(edgeR)
library(limma)
library(Homo.sapiens)
library(reticulate)
```
## Loading in the Data

The working directory is set to the location where the files are stored.
```{r}
setwd("/Users/ramiyasivakumar/")
getwd()
```
The files are stored into a collective location and they are read. The first 5 rows are printed as a result. We can see that the files contain the ENSEMBL Gene Id in one column and the counts in another. The ENSEMBL Gene Id is in the format: ENS(species)(object type)(identifier).(version). 
The files do not contain a header, resulting in the column heads being denoted as V1 and V2.
```{r}
files <- c("M1.tsv","M2.tsv","M3.tsv","M4.tsv","M5.tsv","M6.tsv","M7.tsv","M8.tsv","M9.tsv","M10.tsv","M11.tsv","M12.tsv","M13.tsv","M14.tsv","M15.tsv","M16.tsv","M17.tsv","M18.tsv","M19.tsv","M20.tsv","F1.tsv","F2.tsv","F3.tsv","F4.tsv","F5.tsv","F6.tsv","F7.tsv","F8.tsv","F9.tsv","F10.tsv","F11.tsv","F12.tsv","F13.tsv","F14.tsv","F15.tsv","F16.tsv","F17.tsv","F18.tsv","F19.tsv","F20.tsv")

read.delim(files[1],nrow=5,header=FALSE)
```
The contents of each tsv file is combined into one dataframe called as 'x'. The data is grouped according to the parameter of Male vs. Female. EdgeR provides a simplified method of doing this with a single command. 
```{r}
x <- readDGE(files, columns=c(1,2))
class(x)
```
The dimensions of the dataframe are then printed. Each file contains 60,487 genes and we have a sample size of 40 people (20 Males and 20 Females)
```{r}
dim(x)
```
```{r}
samplenames<-colnames(x)
sampleNames
colnames(x) <- samplenames
```
Organizing the data: The samples are grouped according to their gender.
```{r}
group <- c(rep("Male",20),rep("Female",20))
x$samples$group <- group
x$samples
```
As the "Homo.sapiens" package doesn't contain the ENSEMBL Gene ID in the format provided by the loaded data, the data must be converted using a substitution function. The digits following the decimal are replaced to look more similar to the key provided.
```{r}
library(gsubfn)
geneid <- rownames(x)
geneid <-gsub("\\.[0-9]*$","", geneid)
head(geneid)

```
The dataframe genes is used to store the gene-level information corresponding to the rows of genes in our datasets. The information required for this data frame is retreived from the 'Homo.sapiens' library. The GeneIDs are associated with the symbol and the chromosome number of each gene. 
```{r}
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENSEMBL")
head(genes)
```
The duplicated genes are removed from the dataset. 
```{r}
genes <- genes[!duplicated(genes$ENSEMBL),]
x$genes <- genes
x
```
## Data Pre-Processing 

The data is transformed from its raw form to more useable. The calculations of "cpm: counts per million" and "lcpm: log of counts per million" are carried out to do so.
```{r}
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
```
For further calculations the data is processed using the calculations below.
```{r}
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
```
The summary of the lcpm calculation is presented graphically.
```{r}
summary(lcpm)
```
Some genes are lowly expressed and do not contribute to the statistical analysis in a meaningful manner. To understand the ratio of these genes the following line is run: 
```{r}
table(rowSums(x$counts==0)==9)
```
EdgeR provides an easy method to remove these lowly expressed genes. 
```{r}
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
```
The filtered versus unfiltered data is presented in a graphical manner using the following lines of code. The package "RColorBrewer" allows for the creation of visually appealing graphs. 
```{r}
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
```
For statistical analysis, the normalization of data is greatly beneficial. edgeR offers a function to do so. The resultant values are stored in a dataset called norm.factors. 
```{r}
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
```
To understand this effect the data is manipulated in such a way that the counts of the first sample are reduced to 5% and the second is increased by 500%.
```{r}
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5
```
The plot below shows data that is normalized versus unnormalized. 
```{r}
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
```

```{r}
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
```
The MDS plot is an informative representation of the similarities and dissimilarities in a sample. The extent of the differences can give an idea of the expected results of the differential expression analysis. The plotMDS function of limma allows us to create the graph. The distances on the plot correspond to the action of different categories on this data. 
```{r}
library(limma)
library(Glimma)
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
group
levels(group) <-  brewer.pal(nlevels(group), "Set1")
col.group <- as.character(group)
col.group <- c("purple","orange")[group]
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
```
An interactive version of this graph can be launched using glimma.
```{r}
glMDSPlot(lcpm,groups = group)
```

## Differential Expression Analysis

We are now beginning the Differential expression analysis. The first step is to greate a design matrix with both the groups of data. 
```{r}
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design
```
Contrasts for the two groups are created and stored in a matrix. 
```{r}
contr.matrix <- makeContrasts(
   MaleVsFemale = Male-Female, 
   levels = colnames(design))
contr.matrix
```
Voom is a function within limma that allows the raw counts data to be transformed into the normalized data that accomadates the mean-variance relationship.
```{r}
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
```
```{r}
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
```
To have a sneak peak at the differential expression levels, the following table can be created. (Up- up-regulated genes)
```{r}
summary(decideTests(efit))
```
T-statistics allows for a more stringent definition of signigicance to be applied to the data. 
```{r}
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
```
```{r}
de.common <- which(dt[,1]!=0)
length(de.common)
```
```{r}
head(tfit$genes$SYMBOL[de.common], n=20)
```
```{r}
vennDiagram(dt[,1], circle.col=c("turquoise", "salmon"))
```
```{r}
write.fit(tfit, dt, file="results.txt")
```

```{r}

Male.vs.Female <- topTreat(tfit, coef=1, n=Inf)
head(Male.vs.Female)
```

To summarise the genes, the log fold change is plot against the log counts per million using the plotMD function.
```{r}
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))
```

Glimma allows the same information to be presented interactively. 
```{r}
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENSEMBL", counts=lcpm, groups=group, launch=TRUE)
```

A heatmap is created of the highest expressed genes. Heatmap.plus is a library that allows for the creation of heatmaps for larger datasets. 
```{r}
library(gplots)
library(heatmap.plus)
Male.vs.Female.topgenes <- Male.vs.Female$ENSEMBL[1:50]
i <- which(v$genes$ENSEMBL %in% Male.vs.Female.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
par(cex.main=0.8,mar=c(1,1,1,1))
heatmap.plus(lcpm[i,], col=bluered(20),cexRow=1,cexCol=0.2, margins = c(20,13), main="HeatMap")

```

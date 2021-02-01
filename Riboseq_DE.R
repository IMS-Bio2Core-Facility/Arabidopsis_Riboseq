#Upload library for DiffExpr in RNAseq sample Arabidopsis Thaliana

options(stringsAsFactors = F)
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")
#BiocManager::install("Rsamtools")
library(Rsamtools)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(org.At.tair.db)
library("GO.db")
library("GOstats")
library(GenomicAlignments)
library(xlsx)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(Biostrings)
library(GenomeInfoDb)
library(BSgenome)
library('GenomicFeatures')
library(GenomeInfoDb)
library(GenomicRanges)
library(Biostrings)
library(GenomicFeatures)
library(ggplot2)
library(grid)
library(systemPipeR)
library(S4Vectors)
library(GGally)
library(RiboseQC)
library(limma)
library(Rsubread)
library(edgeR)
library(biomaRt)
library(RColorBrewer)
library(DESeq2, quietly = TRUE)
library(ape, warn.conflicts = FALSE)
library(gplots)

# Create directories and Upload BAM files

WORKDIR="/usr/data/bgfs1/maloku/Riboseq_genome/bam_file/";
#dir.create(WORKDIR)
setwd(WORKDIR)
Results="/usr/data/bgfs1/maloku/Riboseq_genome/Results/"; 
#dir.create(Results)
BAMS <- list.files(WORKDIR, pattern = ".bam$",full.names = T)

# featureCounts and Normalization

#Quantifying read counts for each gene
fc <-featureCounts(files=BAMS, annot.ext="/usr/data/bgfs1/maloku/Riboseq_genome/Arabidopsis_thaliana.TAIR10.45.gtf",
                   useMetaFeatures = F, isPairedEnd = F, isGTFAnnotationFile = T, nthreads = 20, annot.inbuilt="TAIR10")


# Create a  matrix Count and Convert counts to dge object

dgList <- DGEList(fc$counts) 
dgList$genes <- fc$annotation[,c( "GeneID","Length"), drop=FALSE]

# Adding gene annotation

dgList$genes$Symbol <- mapIds(org.At.tair.db, rownames(dgList),keytype="TAIR", column="SYMBOL")

# Filtering 

countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) > 10)
dgList <-dgList[keep,]

# Normalisation with TMM

dgList <-calcNormFactors(dgList)

# Data exploration

plotMDS(dgList)
plotMD(dgList, column=1)
abline(h=0, col="red", lty=2, lwd=2)


# Setting up the  model

targets <-readTargets("/usr/data/bgfs1/maloku/Riboseq_genome/sample_name_ribo.csv",sep=",")
time <-targets$Time
mm <- model.matrix(~0+time)
rownames(mm) <-targets$Sample

# Dispersion estimation  

dgList <- estimateGLMCommonDisp(dgList, design = mm)
dgList <- estimateGLMTrendedDisp(dgList, design = mm)
dgList <- estimateGLMTagwiseDisp(dgList, design = mm)

plotBCV(dgList)


#y <- estimateDisp(dgList, mm, robust=TRUE) 
#y$common.dispersio
#plotBCV(y)

# Differential Expression

fit1 <- glmQLFit(dgList, mm)

# make a contrast 
T1vsT3 <-makeContrasts(timeT1-timeT3, levels=mm)
T1vsT5 <-makeContrasts(timeT1-timeT5, levels=mm)
T1vsT7 <-makeContrasts(timeT1-timeT7, levels=mm)
T1vsT9 <-makeContrasts(timeT1-timeT9, levels=mm)
T1vsT11 <-makeContrasts(timeT1-timeT11, levels=mm)
T1vsT13 <-makeContrasts(timeT1-timeT13, levels=mm)
T1vsT15 <-makeContrasts(timeT1-timeT15, levels=mm)
T1vsT17 <-makeContrasts(timeT1-timeT17, levels=mm)
T1vsT19 <-makeContrasts(timeT1-timeT19, levels=mm)
T1vsT21 <-makeContrasts(timeT1-timeT21, levels=mm)
T1vsT23 <-makeContrasts(timeT1-timeT23, levels=mm)


qlf.T1vsT3 <- glmQLFTest(fit1, contrast=T1vsT3)
topTags(qlf.T1vsT3)
qlf.T1vsT5 <- glmQLFTest(fit1, contrast=T1vsT5)
topTags(qlf.T1vsT5)
qlf.T1vsT7 <- glmQLFTest(fit1, contrast=T1vsT7)
topTags(qlf.T1vsT7)
qlf.T1vsT9 <- glmQLFTest(fit1, contrast=T1vsT9)
topTags(qlf.T1vsT9)
qlf.T1vsT11<- glmQLFTest(fit1, contrast=T1vsT11)
topTags(qlf.T1vsT11)
qlf.T1vsT13 <- glmQLFTest(fit1, contrast=T1vsT13)
topTags(qlf.T1vsT13)
qlf.T1vsT15 <- glmQLFTest(fit1, contrast=T1vsT15)
topTags(qlf.T1vsT15)
qlf.T1vsT17 <- glmQLFTest(fit1, contrast=T1vsT17)
topTags(qlf.T1vsT17)
qlf.T1vsT19 <- glmQLFTest(fit1, contrast=T1vsT19)
topTags(qlf.T1vsT19)
qlf.T1vsT21 <- glmQLFTest(fit1, contrast=T1vsT21)
topTags(qlf.T1vsT21)
qlf.T1vsT23 <- glmQLFTest(fit1, contrast=T1vsT23)
topTags(qlf.T1vsT23)

#filtering criteria for significance based on both p-value (or FDR corrected p-value) and Fold Chance (or log fold change)

FDR <- p.adjust(qlf.T1vsT3$table$PValue, method="BH")
sum(FDR < 0.05)

is.de_T1vsT3 <- decideTestsDGE(qlf.T1vsT3)
summary(is.de_T1vsT3)
is.de_T1vsT5 <- decideTestsDGE(qlf.T1vsT5)
summary(is.de_T1vsT5)
is.de_T1vsT7 <- decideTestsDGE(qlf.T1vsT7)
summary(is.de_T1vsT7)
is.de_T1vsT9 <- decideTestsDGE(qlf.T1vsT9)
summary(is.de_T1vsT9)
is.de_T1vsT11 <- decideTestsDGE(qlf.T1vsT11)
summary(is.de_T1vsT11)
is.de_T1vsT13 <- decideTestsDGE(qlf.T1vsT13)
summary(is.de_T1vsT13)
is.de_T1vsT15 <- decideTestsDGE(qlf.T1vsT15)
summary(is.de_T1vsT15)
is.de_T1vsT17 <- decideTestsDGE(qlf.T1vsT17)
summary(is.de_T1vsT17)
is.de_T1vsT19 <- decideTestsDGE(qlf.T1vsT19)
summary(is.de_T1vsT19)
is.de_T1vsT21 <- decideTestsDGE(qlf.T1vsT21)
summary(is.de_T1vsT21)
is.de_T1vsT23 <- decideTestsDGE(qlf.T1vsT23)
summary(is.de_T1vsT23)

#Plot log-fold change against log-counts per million, with DE genes highlighted

plotMD(qlf.T1vsT13, status=is.de_T1vsT3, values=c(1,-1), col=c("red","blue"), legend="topright")
abline(h=c(-1, 1), col="blue")

#Heat map clustering

#we first convert the read counts into log2-counts-per-million (logCPM) values
#reduce the variability of the logCPM values for genes with low counts.

logCPM <- cpm(dgList, log=TRUE, prior.count=3)
rownames(logCPM) <- dgList$genes$Symbol
colnames(logCPM)<- targets$Sample
o <- order(qlf$table$PValue)
logCPM <- logCPM[o[1:50],]
logCPM <- t(scale(t(logCPM)))
my_palette <- colorRampPalette(c("blue", "white", "red")) 
heatmap.2(logCPM,col=rev(my_palette(70)),trace="none", main="Top 30 most variable genes across samples for T1vsT13",scale="row")


dev.off()


# Gene Ontology
keg <- kegga(qlf, species="At")
topKEGG(keg, n=15, truncate=34)
topKEGG(keg, sort="up")

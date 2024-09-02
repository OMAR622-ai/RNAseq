#Installing packages (this takes > 1 hr, need to be done before the class)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")


#loading packages:
library(limma)
library(edgeR)
#loading data
#mobData_ <- read.delim("mobData.tsv",sep="\t",header=TRUE,row.names=1)
load("mobData.RData")
head(mobData)

#checking the null and NA reads
is.null(mobData)
is.na(mobData)
sum(is.na(mobData))

#to see the distribution of the data :
hist(as.matrix(mobData))
hist(log2(as.matrix(mobData)))

#help(mobData)
mobDataGroups <- c("MM", "MM", "WM", "WM", "WW", "WW")

#creating DGE object
d <- DGEList(counts=mobData,group=factor(mobDataGroups))
d

#filtering
dim(d)
d.full <- d # keep the old one in case we mess up
head(cpm(d))
apply(d$counts, 2, sum) # total gene counts per sample
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
d$samples$lib.size <- colSums(d$counts)
d$samples

#Normalizing
d <- calcNormFactors(d, method="TMM")
d

#Data exploration
plotMDS(d, col=as.numeric(d$samples$group))
legend("topleft", as.character(unique(d$samples$group)), col=1:3, pch=23)

#Estimating the Dispersion
#Classical Quantile-adjusted conditional maximum likelihood (qCML) => single factor
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)


#GLM estimates of dispersion => For general experiments (with multiple factors), edgeR uses the Cox-Reid profile-adjusted
#likelihood (CR) method in estimating dispersions [25]. The CR method is derived to overcome
#the limitations of the qCML method as mentioned above. It takes care of multiple factors by
#fitting generalized linear models (GLM) with a design matrix
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)


#differential expression
#qCML
et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
et13 <- exactTest(d1, pair=c(1,3)) # compare groups 1 and 3
et23 <- exactTest(d1, pair=c(2,3)) # compare groups 2 and 3
topTags(et12, n=10)

de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)
View(de1)
# differential expressed tags from the naive method in d1
de1tags12 <- rownames(d1)[as.logical(de1)]
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")

#GLM testing for differential expression
design.mat
fit <- glmQLFit(d2, design.mat)
# compare (group 1 - group 2) to 0:
lrt12 <- glmQLFTest(fit, contrast=c(1,-1,0))
lrt13 <- glmQLFTest(fit, contrast=c(1,0,-1))
lrt23 <- glmQLFTest(fit, contrast=c(0,1,-1))
topTags(lrt12, n=10)
de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05 )
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt12, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")



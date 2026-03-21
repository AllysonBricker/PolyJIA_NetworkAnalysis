library(GEOquery)
library(Biobase)
library(limma)

#Get GSE67596 dataset
JIA <- getGEO("GSE67596", GSEMatrix =TRUE, AnnotGPL = TRUE)[[1]]

#Ensure syntax of names is valid for R and view phenotype data
fvarLabels(JIA) <- make.names(fvarLabels(JIA))

#Display table of characteristic data 
table(pData(JIA)$characteristics_ch1, pData(JIA)$characteristics_ch1.1)

#characteristics_ch1 and characteristics_ch1.1 will need to be filtered
#characteristics_ch1 describes the cell types
#characteristics_ch1.1 describes control and disease groups

#Filter characteristics_ch1 column by PBMC
JIA <- JIA[, pData(JIA)$characteristics_ch1 == "cell type: PBMC" ]

#Filter out Pauciarticular JIA from characteristics_ch1.1 column
JIA <- JIA[, pData(JIA)$characteristics_ch1.1 != "disease state: Pauciarticular JIA" ]

#Display table of filtered data
table(pData(JIA)$characteristics_ch1, pData(JIA)$characteristics_ch1.1)

#Check normalization
summary(exprs(JIA))
boxplot(exprs(JIA),outline=FALSE)

#Log2 transformed values will be between 0 and 16, so JIA has not been log2 transformed
#Apply log2 transformation
exprs(JIA) <- log2(exprs(JIA))
summary(exprs(JIA))
boxplot(exprs(JIA),outline=FALSE)

#Get GSE45919 dataset
EBV <- getGEO("GSE45919", GSEMatrix =TRUE, AnnotGPL = TRUE)[[1]]

#Ensure syntax of names is valid for R and view phenotype data
fvarLabels(EBV) <- make.names(fvarLabels(EBV))

#Display table of characteristic data 
table(pData(EBV)$characteristics_ch1)

#characteristics_ch1 will need to be filtered
#characteristics_ch1 describes control and disease state groups

#Filter out "following resolution of infectious mononucleosis" from characteristics_ch1 column
EBV <- EBV[, pData(EBV)$characteristics_ch1 != "sample time: following resolution of infectious mononucelosis" ]

#Display table of filtered data
table(pData(EBV)$characteristics_ch1)

#Check normalization
summary(exprs(EBV))
boxplot(exprs(EBV),outline=FALSE)

#Log2 transformed values will be between 0 and 16, so EBV has not been log2 transformed
#Apply log2 transformation
exprs(EBV) <- log2(exprs(EBV))
summary(exprs(EBV))
boxplot(exprs(EBV),outline=FALSE)

#Get GSE34205 dataset
H1N1 <- getGEO("GSE34205", GSEMatrix =TRUE, AnnotGPL = TRUE)[[1]]

#Ensure syntax of names is valid for R and view phenotype data
fvarLabels(H1N1) <- make.names(fvarLabels(H1N1))

#Filter out "infection: rsv" from characteristics_ch1.3 column
H1N1 <- H1N1[, pData(H1N1)$characteristics_ch1.3 != "infection: rsv" ]

#Display table of characteristic data 
table(pData(H1N1)$characteristics_ch1.3)

#Check normalization
summary(exprs(H1N1))
boxplot(exprs(H1N1), outline=FALSE)

#Log2 transformed values will be between 0 and 16, so H1N1 has not been log2 transformed
#Apply log2 transformation
exprs(H1N1) <- log2(exprs(H1N1))
summary(exprs(H1N1))
boxplot(exprs(H1N1),outline=FALSE)

#Group GSE67596
gsms <- "11111111111111000000000000000"
sml <- strsplit(gsms, split="")[[1]]
gs <- factor(sml)
groups <- make.names(c("polyJIA","control"))
levels(gs) <- groups
JIA$group <- gs
design <- model.matrix(~group + 0, JIA)
colnames(design) <- levels(gs)

#Skip missing values
JIA <- JIA[complete.cases(exprs(JIA)), ]

#Fit linear model
JIAlm <- lmFit(JIA, design)

#Set up pairwise comparisons and refit model
cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
JIAlm2 <- contrasts.fit(JIAlm, cont.matrix)

#Empirical Bayes analysis with Benjamini & Hochberg adjustment "fdr"
JIAlm2 <- eBayes(JIAlm2, 0.01)
tT_JIA <- topTable(JIAlm2, adjust="fdr", sort.by="B", number=250, p.value=0.005)

#Filter and display table of results
tT_JIA <- subset(tT_JIA, select=c("ID","adj.P.Val","P.Value","Gene.symbol","Gene.title"))
tT_JIA <- tT_JIA[!is.na(tT_JIA$Gene.symbol) &
                   tT_JIA$Gene.symbol != "" &
                   trimws(tT_JIA$Gene.symbol) != "", ]
View(tT_JIA)

#Group GSE45919
gsms <- "101000100000"
sml <- strsplit(gsms, split="")[[1]]
gs <- factor(sml)
groups <- make.names(c("during","control"))
levels(gs) <- groups
EBV$group <- gs
design <- model.matrix(~group + 0, EBV)
colnames(design) <- levels(gs)

#Skip missing values
EBV <- EBV[complete.cases(exprs(EBV)), ]

#Fit linear model
EBVlm <- lmFit(EBV, design)

#Set up pairwise comparisons and refit model
cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
EBVlm2 <- contrasts.fit(EBVlm, cont.matrix)

#Empirical Bayes analysis with Benjamini & Hochberg adjustment "fdr"
EBVlm2 <- eBayes(EBVlm2, 0.01)
tT_EBV <- topTable(EBVlm2, adjust="fdr", sort.by="B", number=250, p.value=0.005)

#Filter and display table of results
tT_EBV <- subset(tT_EBV, select=c("ID","adj.P.Val","P.Value","Gene.symbol","Gene.title"))
View(tT_EBV)

#Group GSE34205
gsms <- "00000000011111111111100111000000000001111111111111"
sml <- strsplit(gsms, split="")[[1]]
gs <- factor(sml)
groups <- make.names(c("control","infection"))
levels(gs) <- groups
H1N1$group <- gs
design <- model.matrix(~group + 0, H1N1)
colnames(design) <- levels(gs)

#Skip missing values
H1N1 <- H1N1[complete.cases(exprs(H1N1)), ]

#Fit linear model
H1N1lm <- lmFit(H1N1, design)

#Set up pairwise comparisons and refit model
cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
H1N1lm2 <- contrasts.fit(H1N1lm, cont.matrix)

#Empirical Bayes analysis with Benjamini & Hochberg adjustment "fdr"
H1N1lm2 <- eBayes(H1N1lm2, 0.01)
tT_H1N1 <- topTable(H1N1lm2, adjust="fdr", sort.by="B", number=500, p.value=0.005)

#Filter and display table of results
tT_H1N1 <- subset(tT_H1N1, select=c("ID","adj.P.Val","P.Value","Gene.symbol","Gene.title"))
tT_H1N1 <- tT_H1N1[!is.na(tT_H1N1$Gene.symbol) &
                   tT_H1N1$Gene.symbol != "" &
                   trimws(tT_H1N1$Gene.symbol) != "", ]
View(tT_H1N1)
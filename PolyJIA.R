library(GEOquery)
library(Biobase)

#Get GSE67596 dataset
JIA <- getGEO("GSE67596", GSEMatrix =TRUE, AnnotGPL = TRUE)[[1]]

#Ensure syntax of names is valid for R and view phenotype data
fvarLabels(JIA) <- make.names(fvarLabels(JIA))
View(pData(JIA))

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
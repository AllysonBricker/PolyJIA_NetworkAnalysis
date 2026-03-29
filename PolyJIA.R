library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
library(WGCNA)
library(igraph)

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

#Collapse probes to genes
exprJIA <- exprs(JIA)
featureJIA <- pData(featureData(JIA))
geneSymbolsJIA <- featureJIA$Gene.symbol
#Filter out rows with missing gene symbols
keep <- !is.na(geneSymbolsJIA) & geneSymbolsJIA != "" & trimws(geneSymbolsJIA) != ""
exprJIA <- exprJIA[keep, ]
geneSymbolsJIA <- geneSymbolsJIA[keep]
collapsed <- collapseRows(
  exprJIA, 
  rowGroup = geneSymbolsJIA, 
  rowID = rownames(exprJIA), 
  method = 'maxRowVariance'
)
collapsed_JIA <- collapsed$datETcollapsed

# Subset feature data to match the collapsed genes
fd <- fData(JIA)
fd <- fd[match(rownames(collapsed_JIA), fd$Gene.symbol), ]
rownames(fd) <- rownames(collapsed_JIA)

# Rebuild ExpressionSet
JIA <- ExpressionSet(
  assayData = collapsed_JIA,
  phenoData = phenoData(JIA),
  featureData = AnnotatedDataFrame(fd)
)

#Group JIA by disease and control
gsms <- "11111111111111000000000000000"
sml <- strsplit(gsms, split="")[[1]]
gs <- factor(sml)
groups <- make.names(c("polyJIA","control"))
levels(gs) <- groups
JIA$group <- gs
design <- model.matrix(~group + 0, JIA)
colnames(design) <- levels(gs)

#Fit linear model
JIAlm <- lmFit(JIA, design)

#Set up pairwise comparisons and refit model
cont.matrix <- makeContrasts(polyJIA - control, levels=design)
JIAlm2 <- contrasts.fit(JIAlm, cont.matrix)

#Empirical Bayes analysis with Benjamini & Hochberg adjustment "fdr"
JIAlm2 <- eBayes(JIAlm2, 0.01)
tT_JIA <- topTable(JIAlm2, coef = 1, adjust="fdr", sort.by="B", number=1000, p.value=0.05)

#Filter table of results
tT_JIA <- subset(tT_JIA, select=c("ID","adj.P.Val","P.Value","Gene.symbol","Gene.title"))

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(JIAlm2, adjust.method="fdr", p.value=0.05, lfc=0)

#Create volcano plot
volcanoplot(JIAlm2, coef = 1, main = "Volcano Plot of Upregulated and Downregulated JIA Genes", pch = 20,
            highlight = length(which(dT[,1]!=0)), names = rep('+', nrow(JIAlm2)))

#Create a list of gene symbols to use to filter the original JIA gset
JIA_genes <- as.list(tT_JIA$Gene.symbol)
filtered_JIA <- JIA[fData(JIA)$Gene.symbol %in% JIA_genes, ]

#Create correlation matrix
corr_matrix <- cor(t(exprs(filtered_JIA)), method="pearson")
diag(corr_matrix) <- 0

#Create network graph object from adjacency matrix
g_JIA <- graph_from_adjacency_matrix(
  corr_matrix,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

#Assign gene names to vertices
V(g_JIA)$name <- fData(filtered_JIA)$Gene.symbol

#Remove multiple edges and self loops
g_JIA <- simplify(g_JIA, remove.multiple=TRUE, remove.loops=TRUE)

#Convert edge weights to absolute values
E(g_JIA)$weight <- abs(E(g_JIA)$weight)

# Remove edges below absolute Pearson correlation 0.8
g_JIA <- delete_edges(g_JIA, E(g_JIA)[E(g_JIA)$weight<0.8])

#Recompute edgeweights for plotting
edgeweights <- E(g_JIA)$weight

# Remove any vertices remaining that have no edges
g_JIA <- delete_vertices(g_JIA, degree(g_JIA)==0)

#Convert to minimum spanning tree with Prim's algorithm
mst_JIA <- mst(g_JIA, algorithm="prim")

#Plot the network graph
plot(
  mst_JIA,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="JIA Network Graph of All Nodes")

#Subset graph to nodes with 5 or more edges, plus their neighbor nodes
high_deg_nodes <- V(mst_JIA)[degree(mst_JIA) >= 5]
neighbor_nodes <- unlist(neighborhood(mst_JIA, order = 1, nodes = high_deg_nodes))
sub_nodes <- unique(neighbor_nodes)
subgraph_JIA <- induced_subgraph(mst_JIA, sub_nodes)

#Remove multiple edges and self loops
subgraph_JIA <- simplify(subgraph_JIA, remove.multiple=TRUE, remove.loops=TRUE)

#Convert edge weights to absolute values
E(subgraph_JIA)$weight <- abs(E(subgraph_JIA)$weight)

# Remove edges below absolute Pearson correlation 0.8
subgraph_JIA <- delete_edges(subgraph_JIA, E(subgraph_JIA)[E(subgraph_JIA)$weight<0.8])

#Recompute edgeweights for plotting
edgeweights <- E(subgraph_JIA)$weight

# Remove any vertices remaining that have no edges
subgraph_JIA <- delete_vertices(subgraph_JIA, degree(subgraph_JIA)==0)

#Convert to minimum spanning tree with Prim's algorithm
mst_JIA2 <- mst(subgraph_JIA, algorithm="prim")

#Set node color
V(mst_JIA2)$color <- "green"

#Plot the network graph
plot(
  mst_JIA2,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.label.dist=-0.3,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="JIA Network Graph of Nodes with at Least 5 Edges")

degrees <- degree(mst_JIA2)

degree_JIA <- data.frame(
  node = names(degrees),
  degree = as.numeric(degrees),
  stringsAsFactors = FALSE
)

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

#Collapse probes to genes
exprEBV <- exprs(EBV)
featureEBV <- pData(featureData(EBV))
geneSymbolsEBV <- featureEBV$Gene.symbol
#Filter out rows with missing gene symbols
keep <- !is.na(geneSymbolsEBV) & geneSymbolsEBV != "" & trimws(geneSymbolsEBV) != ""
exprEBV <- exprEBV[keep, ]
geneSymbolsEBV <- geneSymbolsEBV[keep]
collapsed <- collapseRows(
  exprEBV, 
  rowGroup = geneSymbolsEBV, 
  rowID = rownames(exprEBV), 
  method = 'maxRowVariance'
)
collapsed_EBV <- collapsed$datETcollapsed

# Subset feature data to match the collapsed genes
fd <- fData(EBV)
fd <- fd[match(rownames(collapsed_EBV), fd$Gene.symbol), ]
rownames(fd) <- rownames(collapsed_EBV)

# Rebuild ExpressionSet
EBV <- ExpressionSet(
  assayData = collapsed_EBV,
  phenoData = phenoData(EBV),
  featureData = AnnotatedDataFrame(fd)
)

#Group GSE45919
gsms <- "101000100000"
sml <- strsplit(gsms, split="")[[1]]
gs <- factor(sml)
groups <- make.names(c("during","control"))
levels(gs) <- groups
EBV$group <- gs
design <- model.matrix(~group + 0, EBV)
colnames(design) <- levels(gs)

#Fit linear model
EBVlm <- lmFit(EBV, design)

#Set up pairwise comparisons and refit model
cont.matrix <- makeContrasts(during - control, levels=design)
EBVlm2 <- contrasts.fit(EBVlm, cont.matrix)

#Empirical Bayes analysis with Benjamini & Hochberg adjustment "fdr"
EBVlm2 <- eBayes(EBVlm2, 0.01)
tT_EBV <- topTable(EBVlm2, coef = 1, adjust="fdr", sort.by="B", number=1000, p.value=0.05)

#Filter and display table of results
tT_EBV <- subset(tT_EBV, select=c("ID","adj.P.Val","P.Value","Gene.symbol","Gene.title"))

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(EBVlm2, adjust.method="fdr", p.value=0.05, lfc=0)

#Create volcano plot
volcanoplot(EBVlm2, coef = 1, main = "Volcano Plot of Upregulated and Downregulated EBV Genes", pch = 20,
            highlight = length(which(dT[,1]!=0)), names = rep('+', nrow(EBVlm2)))

#Create a list of gene symbols to use to filter the original EBV gset
EBV_genes <- as.list(tT_EBV$Gene.symbol)
filtered_EBV <- EBV[fData(EBV)$Gene.symbol %in% EBV_genes, ]

#Create correlation matrix
corr_matrix <- cor(t(exprs(filtered_EBV)), method="pearson")
diag(corr_matrix) <- 0

#Create network graph object from adjacency matrix
g_EBV <- graph_from_adjacency_matrix(
  corr_matrix,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

#Assign gene names to vertices
V(g_EBV)$name <- fData(filtered_EBV)$Gene.symbol

#Remove multiple edges and self loops
g_EBV <- simplify(g_EBV, remove.multiple=TRUE, remove.loops=TRUE)

#Convert edge weights to absolute values
E(g_EBV)$weight <- abs(E(g_EBV)$weight)

# Remove edges below absolute Pearson correlation 0.8
g_EBV <- delete_edges(g_EBV, E(g_EBV)[E(g_EBV)$weight<0.8])

#Recompute edgeweights for plotting
edgeweights <- E(g_EBV)$weight

# Remove any vertices remaining that have no edges
g_EBV <- delete_vertices(g_EBV, degree(g_EBV)==0)

#Convert to minimum spanning tree with Prim's algorithm
mst_EBV <- mst(g_EBV, algorithm="prim")

#Plot the network graph
plot(
  mst_EBV,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="EBV Network Graph of All Nodes")

#Subset graph to nodes with 5 or more edges, plus their neighbor nodes
high_deg_nodes <- V(mst_EBV)[degree(mst_EBV) >= 5]
neighbor_nodes <- unlist(neighborhood(mst_EBV, order = 1, nodes = high_deg_nodes))
sub_nodes <- unique(neighbor_nodes)
subgraph_EBV <- induced_subgraph(mst_EBV, sub_nodes)

#Remove multiple edges and self loops
subgraph_EBV <- simplify(subgraph_EBV, remove.multiple=TRUE, remove.loops=TRUE)

#Convert edge weights to absolute values
E(subgraph_EBV)$weight <- abs(E(subgraph_EBV)$weight)

# Remove edges below absolute Pearson correlation 0.8
subgraph_EBV <- delete_edges(subgraph_EBV, E(subgraph_EBV)[E(subgraph_EBV)$weight<0.8])

#Recompute edgeweights for plotting
edgeweights <- E(subgraph_EBV)$weight

# Remove any vertices remaining that have no edges
subgraph_EBV <- delete_vertices(subgraph_EBV, degree(subgraph_EBV)==0)

#Convert to minimum spanning tree with Prim's algorithm
mst_EBV2 <- mst(subgraph_EBV, algorithm="prim")

#Set node color
V(mst_EBV2)$color <- "cyan"

#Plot the network graph
plot(
  mst_EBV2,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.label.dist=-0.3,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="EBV Network Graph of Nodes with at Least 5 Edges")

degrees <- degree(mst_EBV2)

degree_EBV <- data.frame(
  node = names(degrees),
  degree = as.numeric(degrees),
  stringsAsFactors = FALSE
)

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

#Collapse probes to genes
exprH1N1 <- exprs(H1N1)
featureH1N1 <- pData(featureData(H1N1))
geneSymbolsH1N1 <- featureH1N1$Gene.symbol
#Filter out rows with missing gene symbols
keep <- !is.na(geneSymbolsH1N1) & geneSymbolsH1N1 != "" & trimws(geneSymbolsH1N1) != ""
exprH1N1 <- exprH1N1[keep, ]
geneSymbolsH1N1 <- geneSymbolsH1N1[keep]
collapsed <- collapseRows(
  exprH1N1, 
  rowGroup = geneSymbolsH1N1, 
  rowID = rownames(exprH1N1), 
  method = 'maxRowVariance'
)
collapsed_H1N1 <- collapsed$datETcollapsed

#Collapse probes to genes
exprH1N1 <- exprs(H1N1)
featureH1N1 <- pData(featureData(H1N1))
geneSymbolsH1N1 <- featureH1N1$Gene.symbol
collapsed <- collapseRows(
  exprH1N1, 
  rowGroup = geneSymbolsH1N1, 
  rowID = rownames(exprH1N1), 
  method = 'maxRowVariance'
)
collapsed_H1N1 <- collapsed$datETcollapsed

# Subset feature data to match the collapsed genes
fd <- fData(H1N1)
fd <- fd[match(rownames(collapsed_H1N1), fd$Gene.symbol), ]
rownames(fd) <- rownames(collapsed_H1N1)

# Rebuild ExpressionSet
H1N1 <- ExpressionSet(
  assayData = collapsed_H1N1,
  phenoData = phenoData(H1N1),
  featureData = AnnotatedDataFrame(fd)
)

#Group GSE34205
gsms <- "00000000011111111111100111000000000001111111111111"
sml <- strsplit(gsms, split="")[[1]]
gs <- factor(sml)
groups <- make.names(c("control","infection"))
levels(gs) <- groups
H1N1$group <- gs
design <- model.matrix(~group + 0, H1N1)
colnames(design) <- levels(gs)

#Fit linear model
H1N1lm <- lmFit(H1N1, design)

#Set up pairwise comparisons and refit model
cont.matrix <- makeContrasts(control - infection, levels=design)
H1N1lm2 <- contrasts.fit(H1N1lm, cont.matrix)

#Empirical Bayes analysis with Benjamini & Hochberg adjustment "fdr"
H1N1lm2 <- eBayes(H1N1lm2, 0.01)
tT_H1N1 <- topTable(H1N1lm2, coef = 1, adjust="fdr", sort.by="B", number=1500, p.value=0.05)

#Filter and display table of results
tT_H1N1 <- subset(tT_H1N1, select=c("ID","adj.P.Val","P.Value","Gene.symbol","Gene.title"))
View(tT_H1N1)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(H1N1lm2, adjust.method="fdr", p.value=0.05, lfc=0)

#Create volcano plot
volcanoplot(H1N1lm2, coef = 1, main = "Volcano Plot of Upregulated and Downregulated H1N1 Genes", pch = 20,
            highlight = length(which(dT[,1]!=0)), names = rep('+', nrow(H1N1lm2)))

#Create a list of gene symbols to use to filter the original H1N1 gset
H1N1_genes <- as.list(tT_H1N1$Gene.symbol)
filtered_H1N1 <- H1N1[fData(H1N1)$Gene.symbol %in% H1N1_genes, ]

#Create correlation matrix
corr_matrix <- cor(t(exprs(filtered_H1N1)), method="pearson")
diag(corr_matrix) <- 0

#Create network graph object from adjacency matrix
g_H1N1 <- graph_from_adjacency_matrix(
  corr_matrix,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

#Assign gene names to vertices
V(g_H1N1)$name <- fData(filtered_H1N1)$Gene.symbol

#Remove multiple edges and self loops
g_H1N1 <- simplify(g_H1N1, remove.multiple=TRUE, remove.loops=TRUE)

#Convert edge weights to absolute values
E(g_H1N1)$weight <- abs(E(g_H1N1)$weight)

# Remove edges below absolute Pearson correlation 0.8
g_H1N1 <- delete_edges(g_H1N1, E(g_H1N1)[E(g_H1N1)$weight<0.8])

#Recompute edgeweights for plotting
edgeweights <- E(g_H1N1)$weight

# Remove any vertices remaining that have no edges
g_H1N1 <- delete_vertices(g_H1N1, degree(g_H1N1)==0)

#Convert to minimum spanning tree with Prim's algorithm
mst_H1N1 <- mst(g_H1N1, algorithm="prim")

#Plot the network graph
plot(
  mst_H1N1,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="H1N1 Network Graph of All Nodes")

#Subset graph to nodes with 5 or more edges, plus their neighbor nodes
high_deg_nodes <- V(mst_H1N1)[degree(mst_H1N1) >= 5]
neighbor_nodes <- unlist(neighborhood(mst_H1N1, order = 1, nodes = high_deg_nodes))
sub_nodes <- unique(neighbor_nodes)
subgraph_H1N1 <- induced_subgraph(mst_H1N1, sub_nodes)

#Remove multiple edges and self loops
subgraph_H1N1 <- simplify(subgraph_H1N1, remove.multiple=TRUE, remove.loops=TRUE)

#Convert edge weights to absolute values
E(subgraph_H1N1)$weight <- abs(E(subgraph_H1N1)$weight)

# Remove edges below absolute Pearson correlation 0.8
subgraph_H1N1 <- delete_edges(subgraph_H1N1, E(subgraph_H1N1)[E(subgraph_H1N1)$weight<0.8])

#Recompute edgeweights for plotting
edgeweights <- E(subgraph_H1N1)$weight

# Remove any vertices remaining that have no edges
subgraph_JIA <- delete_vertices(subgraph_H1N1, degree(subgraph_H1N1)==0)

#Convert to minimum spanning tree with Prim's algorithm
mst_H1N12 <- mst(subgraph_H1N1, algorithm="prim")

#Set node color
V(mst_H1N12)$color <- "yellow"

#Plot the network graph
plot(
  mst_H1N12,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.label.dist=-0.3,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="H1N1 Network Graph of Nodes with at Least 5 Edges")

degrees <- degree(mst_H1N12)

degree_H1N1 <- data.frame(
  node = names(degrees),
  degree = as.numeric(degrees),
  stringsAsFactors = FALSE
)

#Compare all tables of all differentially expressed genes
df_EBV <- data.frame(tT_EBV)
df_JIA <- data.frame(tT_JIA)
df_H1N1 <- data.frame(tT_H1N1)
matching_JIA_EBV <- inner_join(df_EBV, df_JIA, by = "Gene.symbol")
matching_JIA_H1N1 <- inner_join(df_H1N1, df_JIA, by = "Gene.symbol")


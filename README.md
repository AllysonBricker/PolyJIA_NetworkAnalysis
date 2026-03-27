# PolyJIA_NetworkAnalysis
The following libraries were used in addition to the standard R libraries: GEOquery (Davis and Meltzer, 2007), Biobase (Huber et al., 2015),  BiocManager (Morgan & Ramos, 2025), dplyr (Wickham et al., 2026), WGCNA (Langfelder & Horvath, 2008) (Langfelder & Horvath, 2012), and igraph (Csárdi & Nepusz, 2006) (Antonov et al., 2023) (Csárdi et al., 2026). The polyarticular JIA, EBV, and H1N1 datasets all required filtering to remove disease state values that are not relevant to this study. The polyarticular JIA dataset also required filtering based on cell type. All three datasets had a log2 transformation applied. This is a recommended step to minimize variation before further analysis (Wilson et al., 2006). The code used to check the min, max, mean, median, 1st quartile, and 3rd quartile of each GEO accession number, log2 transform, and display the boxplot was based on code used in a tutorial by Mark Dunning (Dunning, 2020). The datasets were grouped by healthy and disease conditions and fit with initial linear models. Then pairwise comparisons were set up and the linear model was refit. Empirical Bayes analysis with Benjamini & Hochberg adjustment of p-values was performed to enable the selection of genes with statistically significant differential expression (Barrett et al., 2013). The code for loading the GSE datasets and labels, ensuring the syntax is correct for R, filtering the phenotype data, grouping, linear modeling, and empirical Bayes analysis with Benjamini & Hochberg adjustment was based on the code used by the GEO2R tool in the GEO database (Barrett et al., 2013). The genes with significant p-values (p < 0.05) were selected and correlation matrices were set up. The correlation matrices were then converted to adjacency matrices so the network graph objects could be created. The network graphs were subsetted to nodes with at least 5 edges, plus their neighboring nodes.

References:
Davis, S., & Meltzer, P. S. (2007). GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. Bioinformatics (Oxford, England), 23(14), 1846–1847. https://doi.org/10.1093/bioinformatics/btm254

Huber, W., Carey, V. J., Gentleman, R., Anders, S., Carlson, M., Carvalho, B. S., Bravo, H. C., Davis, S., Gatto, L., Girke, T., Gottardo, R., Hahne, F., Hansen, K. D., Irizarry, R. A., Lawrence, M., Love, M. I., MacDonald, J., Obenchain, V., Oleś, A. K., Pagès, H., … Morgan, M. (2015). Orchestrating high-throughput genomic analysis with Bioconductor. Nature methods, 12(2), 115–121. https://doi.org/10.1038/nmeth.3252

Morgan M, Ramos M (2025). _BiocManager: Access the Bioconductor Project Package Repository_. doi:10.32614/CRAN.package.BiocManager <https://doi.org/10.32614/CRAN.package.BiocManager>, R package version 1.30.27, https://CRAN.R-project.org/package=BiocManager.

Wickham H, François R, Henry L, Müller K, Vaughan D (2026). _dplyr: A Grammar of Data Manipulation_. doi:10.32614/CRAN.package.dplyr <https://doi.org/10.32614/CRAN.package.dplyr>, R package version 1.2.0, <https://CRAN.R-project.org/package=dplyr>.

Langfelder P, Horvath S (2008). “WGCNA: an R package for weighted correlation network analysis.” _BMC Bioinformatics_, 559. <https://link.springer.com/article/10.1186/1471-2105-9-559>.

Langfelder P, Horvath S (2012). “Fast R Functions for Robust Correlations and Hierarchical Clustering.” _Journal of Statistical Software_, *46*(11), 1-17. <https://www.jstatsoft.org/v46/i11/>.

Csárdi G, Nepusz T (2006). “The igraph software package for complex network research.”_InterJournal_, *Complex Systems*, 1695. <https://igraph.org>.
Antonov M, Csárdi G, Horvát S, Müller K, Nepusz T, Noom D, Salmon M, Traag V, Welles BF, Zanini F (2023). “igraph enables fast and robust network analysis across programming languages.” _arXiv preprint arXiv:2311.10260_. doi:10.48550/arXiv.2311.10260  <https://doi.org/10.48550/arXiv.2311.10260>.
Csárdi G, Nepusz T, Traag V, Horvát Sz, Zanini F, Noom D, Müller K, Schoch D, Salmon M (2026). _igraph: Network Analysis and Visualization in R_. doi:10.5281/zenodo.7682609  <https://doi.org/10.5281/zenodo.7682609>, R package version 2.2.2,  <https://CRAN.R-project.org/package=igraph>.

Wilson, C. H., Tsykin, A., Wilkinson, C. R., & Abbott, C. A. (2006). Experimental Design and Analysis of Microarray Data. Applied Mycology and Biotechnology, 6, 1–36. https://doi.org/10.1016/S1874-5334(06)80004-3

Dunning, M. (2020, June 30). Analysing data from GEO - Work in Progress. https://sbc.shef.ac.uk/geo_tutorial/tutorial.nb.html

Barrett, T., Wilhite, S.E., Ledoux, P., Evangelista, C., Kim, I.F., Tomashevsky, M., Marshall, K.A., Phillippy, K.H., Sherman, P.M., Holko, M., Yefanov, A., Lee, H., Zhang, N., Robertson, C.L., Serova, N., Davis, S., Soboleva, A. (2013)  NCBI GEO: archive for functional genomics data sets—update. Nucleic Acids Research, Volume 41, Issue D1, Pages D991–D995. https://doi.org/10.1093/nar/gks1193

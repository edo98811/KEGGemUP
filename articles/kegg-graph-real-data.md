# kegg-graph

If you dont have the package installed you can do so by running:

``` r
install.packages("remotes")     
remotes::install_github("edo98811/KEGGemUP")
```

The load it.

``` r
library("KEGGemUP")
```

KEGG pathways are used in many bioinformatics contexts, this package
package can make it very easy to implement these in the analysis of real
data. A very important part of bioinformatics analyses is both data
integration and visiazionation. KEGGemUP aims to facilitate these tasks
by providing functions to parse KEGG pathways and build graph objects,
the idea is then to use these graph to map on them thre results of
differential expression analyses.

These vignette assumes you have already performed a differential
expression analysis and have the results available as a `data.frame` or
a result object from `limma` or `DESeq2`. To reproduce this situation we
will now build a examplary differential expression results using the
macrophage dataset from the macrophage package.

Setting up the data

KEGG pathways are used in many bioinformatics contexts, this package
package can make it very easy to implement these in the analysis of real
data. A very important part of bioinformatics analyses is both data
integration and visiazionation. KEGGemUP aims to facilitate these tasks
by providing functions to parse KEGG pathways and build graph objects,
the idea is then to use these graph to map on them thre results of
differential expression analyses.

These vignette assumes you have already performed a differential
expression analysis and have the results available as a `data.frame` or
a result object from `limma` or `DESeq2`. To reproduce this situation we
will now build a examplary differential expression results using the
macrophage dataset from the macrophage package.

Setting up the data

``` r
message("--- Loading packages...")
#> --- Loading packages...

suppressPackageStartupMessages({
  library("macrophage")
  library("DESeq2")
  library("org.Hs.eg.db")
  library("AnnotationDbi")
  library("clusterProfiler")
  library("limma")
  library("edgeR")
})
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
message("- Done!")
#> - Done!

# Load the macrophage dataset ---------------------------------------------------
data(gse)
rownames(gse) <- substr(rownames(gse), 1, 15)  # truncate rownames at 15 characters

# DESeq2 analysis ---------------------------------------------------------------
dds <- DESeqDataSet(gse, design = ~ condition)
#> using counts and average transcript lengths from tximeta

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep, ]

dds <- DESeq(dds)
#> estimating size factors
#> using 'avgTxLength' from assays(dds), correcting for library size
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

# limma analysis ---------------------------------------------------------------
condition <- factor(colData(gse)[, "condition_name"])
design <- model.matrix(~0 + condition)

contrast.matrix <- makeContrasts(
  IFNg_vs_naive = conditionIFNg - conditionnaive,
  levels = design
)

dge <- DGEList(assay(gse))
dge <- calcNormFactors(dge)

cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
dge <- dge[-drop, ]

voom_mat <- voom(dge, design, plot = FALSE)
fit <- lmFit(voom_mat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)  # Empirical Bayes moderation

# Gene annotation ---------------------------------------------------------------

anns <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(dds),
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "ENSEMBL",
  multiVals = "first"
)
#> 'select()' returned 1:many mapping between keys and columns

# DESeq2 results ---------------------------------------------------------------
res_macrophage_IFNg_vs_naive_dds <- results(
  dds,
  contrast = c("condition", "IFNg", "naive"),
  lfcThreshold = 1,
  alpha = 0.05
)

summary(res_macrophage_IFNg_vs_naive_dds) 
#> 
#> out of 17806 with nonzero total read count
#> adjusted p-value < 0.05
#> LFC > 1.00 (up)    : 416, 2.3%
#> LFC < -1.00 (down) : 141, 0.79%
#> outliers [1]       : 95, 0.53%
#> low counts [2]     : 0, 0%
#> (mean count < 3)
#> [1] see 'cooksCutoff' argument of ?results
#> [2] see 'independentFiltering' argument of ?results

res_macrophage_IFNg_vs_naive_dds$SYMBOL <- rowData(dds)$SYMBOL
res_macrophage_IFNg_vs_naive_dds$ENTREZID <- anns$ENTREZID[match(rownames(res_macrophage_IFNg_vs_naive_dds), anns$ENSEMBL)]

# limma results ---------------------------------------------------------------
res_macrophage_IFNg_vs_naive_limma <- topTable(
  fit2,
  coef = "IFNg_vs_naive",
  adjust = "fdr",
  number = Inf,
  confint = TRUE
)

res_macrophage_IFNg_vs_naive_limma$ENTREZID <- anns$ENTREZID[match(rownames(res_macrophage_IFNg_vs_naive_limma), anns$ENSEMBL)]

# Enrichment analysis ----------------------------------------------------------
de_entrez_IFNg_vs_naive_genes <- anns $ENTREZID[
  (!is.na(res_macrophage_IFNg_vs_naive_dds$padj)) &
    (res_macrophage_IFNg_vs_naive_dds$padj <= 0.05)
]
```

## Functions to parse KGML files

KGML is the format that KEGG uses to save the pathway structure and it
is what this package interfaces itself with.

You can use these functions to parse KGML files directly. From these you
can build a graph object if you wish to do so and you have expertise
with graph analysis in R. The first thing you can do is download a KGML
file from KEGG. You can do this using the function
`download_kegg_kgml()`. The only parameter that needs to be provided is
the destination directory where to save the file. The function will
return the path to the downloaded file.

``` r
kgml_file <- get_and_cache_kgml("hsa04010", file_name = tempfile(fileext = "xml"))
#> Downloading KGML for hsa04010 ...
```

``` r
nodes_df <- parse_kgml_entries(kgml_file)
#> Parsed 134 nodes from KGML file.
edges_df <- parse_kgml_relations(kgml_file)
#> Parsed 193 edges from KGML file.
```

### The output data.frame frame for nodes

Here you can see the first 5 columns of the dataframe that you get by
parsing the nodes from a KGML file. It is a `data.frame` where each row
represents a node in the KEGG pathway graph. No operations are done on
the data, it is just a direct parsing of the KGML file.

``` r
knitr::kable(head(nodes_df))
```

| name | id  | kegg_name                                                                                   | type     | link                                                                                                                                 | reaction | graphics_name                                                            | label                                                                    | fgcolor  | bgcolor  | graphics_type | x   | y   | width | height | components |
|:-----|:----|:--------------------------------------------------------------------------------------------|:---------|:-------------------------------------------------------------------------------------------------------------------------------------|:---------|:-------------------------------------------------------------------------|:-------------------------------------------------------------------------|:---------|:---------|:--------------|:----|:----|:------|:-------|:-----------|
| 19   | 19  | cpd:C00338                                                                                  | compound | <https://www.kegg.jp/dbget-bin/www_bget?C00338>                                                                                      | NA       | C00338                                                                   | C00338                                                                   | \#000000 | \#FFFFFF | circle        | 138 | 743 | 8     | 8      | NA         |
| 20   | 20  | hsa:5923 hsa:5924                                                                           | gene     | <https://www.kegg.jp/dbget-bin/www_bget?hsa:5923+hsa:5924>                                                                           | NA       | RASGRF1, CDC25, CDC25L, GNRP, GRF1, GRF55, H-GRF55, PP13187, ras-GRF1…   | RASGRF1, CDC25, CDC25L, GNRP, GRF1, GRF55, H-GRF55, PP13187, ras-GRF1…   | \#000000 | \#BFFFBF | rectangle     | 392 | 236 | 46    | 17     | NA         |
| 21   | 21  | hsa:11221 hsa:1843 hsa:1844 hsa:1846 hsa:1847 hsa:1848 hsa:1849 hsa:1850 hsa:1852 hsa:80824 | gene     | <https://www.kegg.jp/dbget-bin/www_bget?hsa:11221+hsa:1843+hsa:1844+hsa:1846+hsa:1847+hsa:1848+hsa:1849+hsa:1850+hsa:1852+hsa:80824> | NA       | DUSP10, MKP-5, MKP5…                                                     | DUSP10, MKP-5, MKP5…                                                     | \#000000 | \#BFFFBF | rectangle     | 789 | 364 | 46    | 17     | NA         |
| 22   | 22  | hsa:1845 hsa:5778 hsa:5801 hsa:84867                                                        | gene     | <https://www.kegg.jp/dbget-bin/www_bget?hsa:1845+hsa:5778+hsa:5801+hsa:84867>                                                        | NA       | DUSP3, VHR…                                                              | DUSP3, VHR…                                                              | \#000000 | \#BFFFBF | rectangle     | 740 | 364 | 46    | 17     | NA         |
| 23   | 23  | hsa:5495                                                                                    | gene     | <https://www.kegg.jp/dbget-bin/www_bget?hsa:5495>                                                                                    | NA       | PPM1B, PP2C-beta, PP2C-beta-X, PP2CB, PP2CBETA, PPC2BETAX                | PPM1B, PP2C-beta, PP2C-beta-X, PP2CB, PP2CBETA, PPC2BETAX                | \#000000 | \#BFFFBF | rectangle     | 596 | 770 | 46    | 17     | NA         |
| 24   | 24  | hsa:356                                                                                     | gene     | <https://www.kegg.jp/dbget-bin/www_bget?hsa:356>                                                                                     | NA       | FASLG, ALPS1B, APT1LG1, APTL, CD178, CD95-L, CD95L, FASL, TNFSF6, TNLG1A | FASLG, ALPS1B, APT1LG1, APTL, CD178, CD95-L, CD95L, FASL, TNFSF6, TNLG1A | \#000000 | \#BFFFBF | rectangle     | 137 | 692 | 46    | 17     | NA         |

### The output data.frame frame for edges

Here you can see the first 5 columns of the dataframe that you get by
parsing the edges from a KGML file. It is a `data.frame` where each row
represents an edge in the KEGG pathway graph. No operations are done on
the data, it is just a direct parsing of the KGML file.

``` r
knitr::kable(head(edges_df))
```

| from | to  | type  | subtype    | rel_value |
|:-----|:----|:------|:-----------|:----------|
| 52   | 48  | PPrel | activation | –\>       |
| 47   | 49  | PPrel | activation | –\>       |
| 20   | 49  | PPrel | activation | –\>       |
| 50   | 49  | PPrel | activation | –\>       |
| 45   | 41  | PPrel | activation | –\>       |
| 111  | 108 | PPrel | activation | –\>       |

## Build a graph from a pathway ID and map results to nodes

To map the differential expression results to the nodes of a KEGG
pathway graph you can use
[`map_results_to_graph()`](https://edo98811.github.io/KEGGemUP/reference/map_results_to_graph.md).
The input of this function is a graph object built with
[`kegg_to_graph()`](https://edo98811.github.io/KEGGemUP/reference/kegg_to_graph.md)
and a list of differential expression results tables or a single
differential expression results table.

You can also pass as input to
[`map_results_to_graph()`](https://edo98811.github.io/KEGGemUP/reference/map_results_to_graph.md)
a single `data.frame` with the differential expression results. This
dataframe must contain at least two columns: one with the KEGG feature
IDs (without organism prefix) and another with the values to map to the
nodes (e.g., log2 fold changes). You can pass to the function the names
of these columns if they differ from the default ones with the
parameters `feature_column` and `value_column`.

Note that the KEGG IDs without organism prefix are the the ENTREZ IDs
for genes. For other feature types (e.g., compounds) you will need to
make sure that the IDs in your differential expression results table
match the KEGG IDs used in the graph.

The defualt parameters are:

- `feature_column`: “KEGG_ids”
- `value_column`: “log2FoldChange”

If you have multiple differential expression results tables (for example
if you have one metabolomics and one transcriptomics) to map to the
nodes you can pass a list of lists. Each sublist must contain the
following elements:

- `de_table`: a data.frame with the differential expression results.
- `value_column`: the name of the column in de_table containing the
  values to map to the nodes.
- `feature_column`: the name of the column in de_table containing the
  feature IDs (e.g., ENTREZ IDs) that correspond to the KEGG ids in the
  graph (without organism prefix).

``` r
de_results_list <-list(
  trans_limma = list(
    de_table = data.frame(res_macrophage_IFNg_vs_naive_limma),
    value_column = "logFC",
    feature_column = "ENTREZID"
  ),
  trans_deseq = list(
    de_table = data.frame(res_macrophage_IFNg_vs_naive_dds),
    value_column = "log2FoldChange",
    feature_column = "ENTREZID"
    )
)
```

### Example of usage with the list of DE results tables

We will now take a KEGG pathway from the enrichment results we built
earlier and map the differential expression results to its nodes. As you
can see you simply need to pass to the function the graph object and the
list of differential expression results tables that we defined before.

``` r
pathway <- "hsa00563"
graph <- kegg_to_graph(pathway, return_type = "igraph")
#> Downloading KGML for hsa00563 ...
#> Downloaded & cached: hsa00563
#> Parsed 126 nodes from KGML file.
#> Parsed 35 edges from KGML file.
#> Downloading KEGG compounds...
#> Downloading KEGG glycans...
graph_visnetwork <- map_results_to_graph(graph, de_results_list, return_type = "visNetwork")
#> Mapping differential expression results to nodes...
#> Warning in add_results_nodes(nodes_df, results_combined): Some nodes had
#> multiple matching KEGG IDs; only the first match was assigned a value.
graph_visnetwork
```

### Example of usage with a single DE results table

Let’s first build a filtered differential expression results table with
only the significant results.

``` r
de_results_limma <- res_macrophage_IFNg_vs_naive_limma
```

Then we will call the function with this table as input. The parameter
`feature_column` and `value_column` are set to match the column names in
our differential expression results table.

``` r
graph <- kegg_to_graph(pathway, return_type = "igraph")
#> Using cached KEGG KGML for hsa00563
#> Parsed 126 nodes from KGML file.
#> Parsed 35 edges from KGML file.
#> Loading KEGG compounds from cache...
#> Loading KEGG glycans from cache...
graph_visnetwork <- map_results_to_graph(graph, de_results_limma, feature_column = "ENTREZID", value_column = "logFC", return_type = "visNetwork")
#> Mapping differential expression results to nodes...
graph_visnetwork
```

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] edgeR_4.8.0                 limma_3.66.0               
#>  [3] clusterProfiler_4.18.2      org.Hs.eg.db_3.22.0        
#>  [5] AnnotationDbi_1.72.0        DESeq2_1.50.2              
#>  [7] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#>  [9] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#> [11] GenomicRanges_1.62.0        Seqinfo_1.0.0              
#> [13] IRanges_2.44.0              S4Vectors_0.48.0           
#> [15] BiocGenerics_0.56.0         generics_0.1.4             
#> [17] macrophage_1.26.0           KEGGemUP_0.1.0             
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3      jsonlite_2.0.0          tidydr_0.0.6           
#>   [4] magrittr_2.0.4          ggtangle_0.0.9          farver_2.1.2           
#>   [7] rmarkdown_2.30          fs_1.6.6                ragg_1.5.0             
#>  [10] vctrs_0.6.5             memoise_2.0.1           ggtree_4.0.1           
#>  [13] htmltools_0.5.9         S4Arrays_1.10.1         curl_7.0.0             
#>  [16] SparseArray_1.10.6      gridGraphics_0.5-1      sass_0.4.10            
#>  [19] bslib_0.9.0             htmlwidgets_1.6.4       desc_1.4.3             
#>  [22] plyr_1.8.9              httr2_1.2.2             cachem_1.1.0           
#>  [25] igraph_2.2.1            lifecycle_1.0.4         pkgconfig_2.0.3        
#>  [28] gson_0.1.0              Matrix_1.7-4            R6_2.6.1               
#>  [31] fastmap_1.2.0           digest_0.6.39           aplot_0.2.9            
#>  [34] enrichplot_1.30.4       ggnewscale_0.5.2        patchwork_1.3.2        
#>  [37] textshaping_1.0.4       RSQLite_2.4.5           filelock_1.0.3         
#>  [40] polyclip_1.10-7         httr_1.4.7              abind_1.4-8            
#>  [43] compiler_4.5.2          withr_3.0.2             bit64_4.6.0-1          
#>  [46] fontquiver_0.2.1        S7_0.2.1                BiocParallel_1.44.0    
#>  [49] DBI_1.2.3               ggforce_0.5.0           R.utils_2.13.0         
#>  [52] MASS_7.3-65             rappdirs_0.3.3          DelayedArray_0.36.0    
#>  [55] tools_4.5.2             otel_0.2.0              scatterpie_0.2.6       
#>  [58] ape_5.8-1               R.oo_1.27.1             glue_1.8.0             
#>  [61] nlme_3.1-168            GOSemSim_2.36.0         grid_4.5.2             
#>  [64] cluster_2.1.8.1         reshape2_1.4.5          fgsea_1.36.0           
#>  [67] gtable_0.3.6            R.methodsS3_1.8.2       tidyr_1.3.1            
#>  [70] data.table_1.17.8       xml2_1.5.1              XVector_0.50.0         
#>  [73] ggrepel_0.9.6           pillar_1.11.1           stringr_1.6.0          
#>  [76] yulab.utils_0.2.2       splines_4.5.2           tweenr_2.0.3           
#>  [79] dplyr_1.1.4             treeio_1.34.0           BiocFileCache_3.0.0    
#>  [82] lattice_0.22-7          bit_4.6.0               tidyselect_1.2.1       
#>  [85] fontLiberation_0.1.0    GO.db_3.22.0            locfit_1.5-9.12        
#>  [88] Biostrings_2.78.0       knitr_1.50              fontBitstreamVera_0.1.1
#>  [91] xfun_0.54               statmod_1.5.1           visNetwork_2.1.4       
#>  [94] stringi_1.8.7           lazyeval_0.2.2          ggfun_0.2.0            
#>  [97] yaml_2.3.12             evaluate_1.0.5          codetools_0.2-20       
#> [100] gdtools_0.4.4           tibble_3.3.0            qvalue_2.42.0          
#> [103] ggplotify_0.1.3         cli_3.6.5               systemfonts_1.3.1      
#> [106] jquerylib_0.1.4         Rcpp_1.1.0              dbplyr_2.5.1           
#> [109] png_0.1-8               parallel_4.5.2          pkgdown_2.2.0          
#> [112] ggplot2_4.0.1           blob_1.2.4              DOSE_4.4.0             
#> [115] tidytree_0.4.6          ggiraph_0.9.2           scales_1.4.0           
#> [118] purrr_1.2.0             crayon_1.5.3            rlang_1.1.6            
#> [121] cowplot_1.2.0           fastmatch_1.1-6         KEGGREST_1.50.0
```

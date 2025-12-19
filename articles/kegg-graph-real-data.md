# The \`KEGGemUP\` package: building KEGG pathway graphs and mapping real data results

**Compiled date**: 2025-12-19

**Last edited**: 15-12-2025

**License**: MIT + file LICENSE

## Installation

The first step is installing the package. If you dont have the package
installed you can do so by running:

``` r
install.packages("remotes")     
remotes::install_github("edo98811/KEGGemUP")
```

Then load it.

``` r
library("KEGGemUP")
```

KEGG pathways are used in many bioinformatics contexts, this package
package can make it very easy to implement these in the analysis of real
data.

A very important part of bioinformatics analyses is both data
integration and visiazionation. KEGGemUP aims to facilitate these tasks
by providing functions to parse KEGG pathways and build graph objects,
the idea is then to use these graph to map on them thre results of
differential expression analyses.

KEGG database uses pathway ids to identify each pathway, for example
“hsa00563” is the KEGG pathway ID for “Glycosylphosphatidylinositol
(GPI)-anchor biosynthesis” in humans. You can use these pathway IDs to
download the KGML files for each pathway and then parse them to build
graph objects. The id must have the organism prefix, for example “hsa”
for human pathways, “mmu” for mouse pathways, etc.

These vignette assumes you have already performed a differential
expression analysis and have the results available as a `data.frame`
(which you can easily obtain from edgeR or limma by running
[`as.data.frame()`](https://rdrr.io/pkg/BiocGenerics/man/as.data.frame.html)
on the table of the results for each object). To reproduce this
situation we will now build an example of a differential expression
result using the macrophage dataset from the macrophage package:
[bioconductor
page](https://www.bioconductor.org/packages/release/data/experiment/html/macrophage.html)

To learn more about graph manipulation you can refer to the igraph and
visNetwork documentation:

- [igraph](https://r.igraph.org) Is a package for creating and
  manipulating graphs in R.
- [visNetwork](https://datastorm-open.github.io/visNetwork/) Is a
  package for interactive visualization of graphs in R, an R interface
  to the javascript library vis.js.

``` r
message("--- Loading packages...")
#> --- Loading packages...

suppressPackageStartupMessages({
  library("macrophage")
  library("org.Hs.eg.db")
  library("SummarizedExperiment")
  library("AnnotationDbi")
  library("clusterProfiler")
  library("limma")
  library("edgeR")
})
message("- Done!")
#> - Done!

# Load the macrophage dataset ---------------------------------------------------
data(gse)
rownames(gse) <- substr(rownames(gse), 1, 15)  # truncate rownames at 15 characters

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
  keys = rownames(gse),
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "ENSEMBL",
  multiVals = "first"
)
#> 'select()' returned 1:many mapping between keys and columns

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
  (!is.na(res_macrophage_IFNg_vs_naive_limma$adj.P.Val)) &
    (res_macrophage_IFNg_vs_naive_limma$adj.P.Valj <= 0.05)
]
```

### Functions to parse KGML files

KGML is the format that KEGG uses to save the pathway structure and it
is what this package interfaces itself with. There is an exported
function to download KGML files from KEGG,
[`download_kgml()`](https://edo98811.github.io/KEGGemUP/reference/download_kgml.md),
which you can use to get the KGML file for a given pathway ID. You can
pass to the function the KEGG pathway ID and the directory where to save
the file.

As an example we will download the KGML file for the KEGG pathway
“Glycosylphosphatidylinositol (GPI)-anchor biosynthesis” in humans,
which has the KEGG pathway ID “hsa00563”.

``` r
kgml_file <- download_kgml("hsa00563", directory = tempdir())  # KEGG pathway ID for "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis"
#> Downloaded & saved in: /tmp/Rtmpy2WsO4/hsa00563.xml
kgml_file
#> [1] "/tmp/Rtmpy2WsO4/hsa00563.xml"
```

You can use these functions to parse KGML files directly. From these you
can build a graph object if you wish to do so and you have expertese
with graph analysis in R. The first function is
[`parse_kgml_entries()`](https://edo98811.github.io/KEGGemUP/reference/parse_kgml_entries.md)
which parses the nodes of the pathway graph from a KGML file. The second
function is
[`parse_kgml_relations()`](https://edo98811.github.io/KEGGemUP/reference/parse_kgml_relations.md)
which parses the edges of the pathway graph from a KGML file.

``` r
nodes_df <- parse_kgml_entries(kgml_file)
#> Parsed 126 nodes from KGML file.
edges_df <- parse_kgml_edges(kgml_file)
#> Parsed 35 relationship edges from KGML file.
#> Parsed 110 reaction edges from KGML file.
#> Total edges parsed from KGML file: 145
```

#### The output data.frame frame for nodes

Here you can see the first 5 columns of the dataframe that you get by
parsing the nodes from a KGML file. It is a `data.frame` where each row
represents a node in the KEGG pathway graph. No operations are done on
the data, it is just a direct parsing of the KGML file.

``` r
knitr::kable(head(nodes_df))
```

| name | id  | kegg_name | type | link                                               | reaction  | graphics_name                                         | label                                                 | fgcolor  | bgcolor  | graphics_type | x   | y   | width | height | components | plot_value | source | color    | text | group | fixed | widthConstraint | heightConstraint | size | shape | KEGG  |
|:-----|:----|:----------|:-----|:---------------------------------------------------|:----------|:------------------------------------------------------|:------------------------------------------------------|:---------|:---------|:--------------|:----|:----|:------|:-------|:-----------|-----------:|:-------|:---------|:-----|:------|:------|----------------:|-----------------:|-----:|:------|:------|
| 13   | 13  | hsa:84992 | gene | <https://www.kegg.jp/dbget-bin/www_bget?hsa:84992> | rn:R05916 | PIGY, HPMRS6, PIG-Y                                   | PIGY, HPMRS6, PIG-Y                                   | \#000000 | \#BFFFBF | rectangle     | 233 | 206 | 46    | 17     | NA         |         NA | NA     | \#BFFFBF |      | NA    | FALSE |              NA |               NA |   NA | NA    | 84992 |
| 14   | 14  | hsa:8818  | gene | <https://www.kegg.jp/dbget-bin/www_bget?hsa:8818>  | rn:R05916 | DPM2, CDG1U                                           | DPM2, CDG1U                                           | \#000000 | \#BFFFBF | rectangle     | 211 | 223 | 46    | 17     | NA         |         NA | NA     | \#BFFFBF |      | NA    | FALSE |              NA |               NA |   NA | NA    | 8818  |
| 15   | 15  | hsa:9091  | gene | <https://www.kegg.jp/dbget-bin/www_bget?hsa:9091>  | rn:R05916 | PIGQ, DEE77, EIEE77, GPI1, GPIBD19, MCAHS4, c407A10.1 | PIGQ, DEE77, EIEE77, GPI1, GPIBD19, MCAHS4, c407A10.1 | \#000000 | \#BFFFBF | rectangle     | 187 | 206 | 46    | 17     | NA         |         NA | NA     | \#BFFFBF |      | NA    | FALSE |              NA |               NA |   NA | NA    | 9091  |
| 16   | 16  | hsa:51227 | gene | <https://www.kegg.jp/dbget-bin/www_bget?hsa:51227> | rn:R05916 | PIGP, DCRC, DCRC-S, DEE55, DSCR5, DSRC, EIEE55, PIG-P | PIGP, DCRC, DCRC-S, DEE55, DSCR5, DSRC, EIEE55, PIG-P | \#000000 | \#BFFFBF | rectangle     | 233 | 189 | 46    | 17     | NA         |         NA | NA     | \#BFFFBF |      | NA    | FALSE |              NA |               NA |   NA | NA    | 51227 |
| 17   | 17  | hsa:5283  | gene | <https://www.kegg.jp/dbget-bin/www_bget?hsa:5283>  | rn:R05916 | PIGH, GPI-H                                           | PIGH, GPI-H                                           | \#000000 | \#BFFFBF | rectangle     | 187 | 189 | 46    | 17     | NA         |         NA | NA     | \#BFFFBF |      | NA    | FALSE |              NA |               NA |   NA | NA    | 5283  |
| 18   | 18  | hsa:5279  | gene | <https://www.kegg.jp/dbget-bin/www_bget?hsa:5279>  | rn:R05916 | PIGC, GPI2, GPIBD16, MRT62                            | PIGC, GPI2, GPIBD16, MRT62                            | \#000000 | \#BFFFBF | rectangle     | 233 | 172 | 46    | 17     | NA         |         NA | NA     | \#BFFFBF |      | NA    | FALSE |              NA |               NA |   NA | NA    | 5279  |

#### The output data.frame frame for edges

Here you can see the first 5 columns of the dataframe that you get by
parsing the edges from a KGML file. It is a `data.frame` where each row
represents an edge in the KEGG pathway graph. No operations are done on
the data, it is just parsing of the KGML file while keeping all the
information from it.

``` r
knitr::kable(head(edges_df))
```

| from | to  | type  | relation_subtype | relation_value | title | width | color | arrows | dashes | label | reaction_id | reaction_name | reaction_type | from_name | to_name |
|:-----|:----|:------|:-----------------|:---------------|:------|------:|:------|:-------|:-------|:------|:------------|:--------------|:--------------|:----------|:--------|
| 37   | 44  | ECrel | compound         | 28             | NA    |     1 | gray  | to     | FALSE  |       | NA          | NA            | NA            | NA        | NA      |
| 37   | 46  | ECrel | compound         | 28             | NA    |     1 | gray  | to     | FALSE  |       | NA          | NA            | NA            | NA        | NA      |
| 44   | 46  | ECrel | compound         | 28             | NA    |     1 | gray  | to     | FALSE  |       | NA          | NA            | NA            | NA        | NA      |
| 191  | 46  | ECrel | compound         | 26             | NA    |     1 | gray  | to     | FALSE  |       | NA          | NA            | NA            | NA        | NA      |
| 42   | 47  | ECrel | compound         | 31             | NA    |     1 | gray  | to     | FALSE  |       | NA          | NA            | NA            | NA        | NA      |
| 39   | 192 | ECrel | compound         | 34             | NA    |     1 | gray  | to     | FALSE  |       | NA          | NA            | NA            | NA        | NA      |

### Build a graph from a pathway ID and map results to nodes

To map the differential expression results to the nodes of a KEGG
pathway graph you can use `map_results_to_nodes()`. The input of this
function is a graph object built with
[`kegg_to_graph()`](https://edo98811.github.io/KEGGemUP/reference/kegg_to_graph.md)
and a list of differential expression results tables or a single
differential expression results table.

You can also pass as input to `map_results_to_nodes()` a single
`data.frame` with the differential expression results. This dataframe
must contain at least two columns: one with the KEGG feature IDs
(without organism prefix) and another with the values to map to the
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

This is an example of how to build such a list of differential
expression results tables. In this case there is only one element, but
ideally you would have one for each omics layer or source you want to
map to the graph nodes.

``` r
de_results_list <-list(
  trans_limma = list(
    de_table = data.frame(res_macrophage_IFNg_vs_naive_limma),
    value_column = "logFC",
    feature_column = "ENTREZID"
  )
)
```

#### Example of usage with the list of DE results tables

We will now take a KEGG pathway from the enrichment results we built
earlier and map the differential expression results to its nodes. As you
can see you simply need to pass to the function the graph object and the
list of differential expression results tables that we defined before.

The output of the first fucntion should be an igraph object, which is
what we specify with the parameter `return_type = "igraph"`. It can also
be a visNetwork object if you set `return_type = "visNetwork"`. The
visNetwork object is useful for interactive visualization of the graph,
but to map the results to the nodes we need to start from an igraph
object. The output of the second function will be a visNetwork object
with the differential expression results mapped to the nodes. In this
case we set `return_type = "visNetwork"` to get a visNetwork object as
output, which then we can visualize directly.

``` r
pathway <- "hsa00563"  
graph <- kegg_to_graph(pathway, return_type = "igraph")
#> adding rname 'https://rest.kegg.jp/get/hsa00563/kgml'
#> Downloaded & cached: hsa00563
#> Parsed 126 nodes from KGML file.
#> Parsed 35 relationship edges from KGML file.
#> Parsed 110 reaction edges from KGML file.
#> Total edges parsed from KGML file: 145
#> adding rname 'https://rest.kegg.jp/list/compound'
#> adding rname 'https://rest.kegg.jp/list/glycan'
graph_visnetwork <- map_results_to_graph(graph, de_results_list, return_type = "visNetwork")
#> Mapping differential expression results to nodes...
graph_visnetwork
```

We could also have the first function return a visNetwork object
directly by setting `return_type = "visNetwork"` in
[`kegg_to_graph()`](https://edo98811.github.io/KEGGemUP/reference/kegg_to_graph.md).
In this case we can immediately visualize the graph but to perfom the
mapping of the differential expression results to the nodes we need an
igraph object as input to
[`map_results_to_graph()`](https://edo98811.github.io/KEGGemUP/reference/map_results_to_graph.md).

#### Example of usage with a single DE results table

Let’s first build a filtered differential expression results table with
only the significant results.

``` r
de_results_limma <-  data.frame(res_macrophage_IFNg_vs_naive_limma)
```

Then we will call the function with this table as input. The parameter
`feature_column` and `value_column` are set to match the column names in
our differential expression results table.

``` r
graph <- kegg_to_graph(pathway, return_type = "igraph")
#> Downloaded & cached: hsa00563
#> Parsed 126 nodes from KGML file.
#> Parsed 35 relationship edges from KGML file.
#> Parsed 110 reaction edges from KGML file.
#> Total edges parsed from KGML file: 145
graph_visnetwork <- map_results_to_graph(graph, de_results_limma, feature_column = "ENTREZID", value_column = "logFC", return_type = "visNetwork")
#> Mapping differential expression results to nodes...
graph_visnetwork
```

You can also control the palette that is used to map the values to
colors on the nodes with the parameter `palette`. The default is “RdBu”
from RColorBrewer, but you can use any palette supported by
RColorBrewer. To see the available palettes you can run
[`RColorBrewer::display.brewer.all()`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html).
You can also visit this page: [RColorBrewer
palettes](https://r-graph-gallery.com/38-rcolorbrewers-palettes.html).

``` r
graph <- kegg_to_graph(pathway, return_type = "igraph")
#> Downloaded & cached: hsa00563
#> Parsed 126 nodes from KGML file.
#> Parsed 35 relationship edges from KGML file.
#> Parsed 110 reaction edges from KGML file.
#> Total edges parsed from KGML file: 145
graph_visnetwork <- map_results_to_graph(graph, de_results_limma, feature_column = "ENTREZID", value_column = "logFC", return_type = "visNetwork", palette = "PiYG")
#> Mapping differential expression results to nodes...
graph_visnetwork
```

### Session info

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
#>  [1] edgeR_4.8.1                 limma_3.66.0               
#>  [3] clusterProfiler_4.18.4      SummarizedExperiment_1.40.0
#>  [5] GenomicRanges_1.62.1        Seqinfo_1.0.0              
#>  [7] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#>  [9] org.Hs.eg.db_3.22.0         AnnotationDbi_1.72.0       
#> [11] IRanges_2.44.0              S4Vectors_0.48.0           
#> [13] Biobase_2.70.0              BiocGenerics_0.56.0        
#> [15] generics_0.1.4              macrophage_1.26.0          
#> [17] KEGGemUP_0.1.0             
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3      jsonlite_2.0.0          tidydr_0.0.6           
#>   [4] magrittr_2.0.4          ggtangle_0.0.9          farver_2.1.2           
#>   [7] rmarkdown_2.30          fs_1.6.6                ragg_1.5.0             
#>  [10] vctrs_0.6.5             memoise_2.0.1           ggtree_4.0.1           
#>  [13] htmltools_0.5.9         S4Arrays_1.10.1         curl_7.0.0             
#>  [16] SparseArray_1.10.7      gridGraphics_0.5-1      sass_0.4.10            
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
#>  [76] yulab.utils_0.2.3       splines_4.5.2           tweenr_2.0.3           
#>  [79] dplyr_1.1.4             treeio_1.34.0           BiocFileCache_3.0.0    
#>  [82] lattice_0.22-7          bit_4.6.0               tidyselect_1.2.1       
#>  [85] locfit_1.5-9.12         fontLiberation_0.1.0    GO.db_3.22.0           
#>  [88] Biostrings_2.78.0       knitr_1.50              fontBitstreamVera_0.1.1
#>  [91] xfun_0.55               statmod_1.5.1           visNetwork_2.1.4       
#>  [94] stringi_1.8.7           lazyeval_0.2.2          ggfun_0.2.0            
#>  [97] yaml_2.3.12             evaluate_1.0.5          codetools_0.2-20       
#> [100] gdtools_0.4.4           tibble_3.3.0            qvalue_2.42.0          
#> [103] BiocManager_1.30.27     ggplotify_0.1.3         cli_3.6.5              
#> [106] systemfonts_1.3.1       jquerylib_0.1.4         Rcpp_1.1.0             
#> [109] dbplyr_2.5.1            png_0.1-8               parallel_4.5.2         
#> [112] pkgdown_2.2.0           ggplot2_4.0.1           blob_1.2.4             
#> [115] DOSE_4.4.0              tidytree_0.4.6          ggiraph_0.9.2          
#> [118] scales_1.4.0            purrr_1.2.0             crayon_1.5.3           
#> [121] BiocStyle_2.38.0        rlang_1.1.6             cowplot_1.2.0          
#> [124] fastmatch_1.1-6         KEGGREST_1.50.0
```

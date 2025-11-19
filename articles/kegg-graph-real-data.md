# kegg-graph

``` r
library("KEGGemUP")
```

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

res_enrich_IFNg_vs_naive_dds <- enrichKEGG(
  gene = de_entrez_IFNg_vs_naive_genes,
  organism = "hsa",
  pvalueCutoff = 0.05
)@result
#> Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...
#> Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...
```

## Functions to parse KGML files

You can use these functions to parse KGML files directly. From these you
can build a graph object if you wish to do so.

``` r
nodes_df <- parse_kgml_entries(system.file("extdata", "hsa04010.xml", package = "KEGGemUP"))
edges_df <- parse_kgml_relations(system.file("extdata", "hsa04010.xml", package = "KEGGemUP"))
```

### The output tibble frame for nodes

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

### The output tibble frame for edges

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

You can build a graph object from a KEGG pathway ID and map your
differential expression results to the nodes of the graph. To do so, you
can provide a list of differentially expressed results tables, each with
the following structure:

- de_table: a data.frame with the differential expression results
- value_column: the name of the column in de_table containing the values
  to map to the nodes
- feature_column: the name of the column in de_table containing the
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

You can also provide a single differential expression results table as a
data.frame. For that the column names need to respect the default
values: - “KEGG_ids”: the column containing the KEGG feature IDs
(without organism prefix) - “log2FoldChange”: the column containing the
values to map to the nodes

``` r
de_results_limma <-  data.frame(res_macrophage_IFNg_vs_naive_limma)[res_macrophage_IFNg_vs_naive_limma$adj.P.Val < 0.05, ]
de_results_limma$log2FoldChange <- de_results_limma$logFC
de_results_limma$KEGG_ids <- de_results_limma$ENTREZID
```

### Example of usage

``` r
pathway <- rownames(res_enrich_IFNg_vs_naive_dds)[6]
graph <- kegg_to_graph(pathway)
#> Downloading KGML for hsa00563...
#> Downloaded & cached: hsa00563
#> Downloading KEGG compounds...
#> Downloading KEGG glycans...
graph <- map_results_to_nodes(graph, de_results_list)
#> Mapping differential expression results to nodes...
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:8705. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:8705. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:8818. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:8733. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:8733. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:2822. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:65258. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:65258. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:65258. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:65258. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:65258. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:80055. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:80055. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:27315. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:93210. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:93210. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:84302. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:84302. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:5277. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:5277. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:5277. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:9488. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:9488. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:9488. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:5279. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:5281. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:5281. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:5281. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:5281. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:54872. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:54872. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:54872. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:54872. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:5283. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:10026. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:10026. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:10026. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:10026. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:9487. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:93183. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:93183. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:93183. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:23556. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:23556. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:23556. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:84720. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:84720. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:84720. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:84720. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:51227. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:9091. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:94005. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:94005. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:51604. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:51604. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:128869. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:128869. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:55650. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:55650. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:55650. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:284098. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:54965. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:80235. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:80235. Keeping the first occurrence to color the node.
#> Warning in add_results_nodes(nodes_df, results_combined): Multiple results
#> mapped to node hsa:80235. Keeping the first occurrence to color the node.
graph
```

``` r
graph <- kegg_to_graph(pathway)
#> Using cached KEGG KGML for hsa00563
#> Loading KEGG compounds from cache...
#> Loading KEGG glycans from cache...
graph <- map_results_to_nodes(graph, de_results_limma)
#> Mapping differential expression results to nodes...
#> Warning in map_results_to_nodes(graph, de_results_limma): Using defaults. For
#> personalisation use a named list of de results.
graph
```

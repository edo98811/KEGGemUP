<!-- badges: start -->
[![R-CMD-check](https://github.com/edo98811/KEGGemUP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/edo98811/KEGGemUP/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

  
# KEGGemUP

A package to map you differential expression results on KEGG pathways.

## Installation
To install the package you need an up to date R version (R >= 4.5.0) and devtools or remotes package installed. Then you can install KEGGemUP from GitHub:
```R
install.packages("remotes")     
remotes::install_github("edo98811/KEGGemUP")
```
## Usage

To use the Kegg pathway visualization function you only need a kegg id of the pathway you want to visualize. You can run the function kegg_to_graph to get the igraph representation of the pathway. If you wish to plot it directly to visnetwork then it is enough to add the parameter `return_type = "visnetwork"`.

```R
library(KEGGemUP)

pathway <- "hsa04110"  # Example pathway ID
graph <- kegg_to_graph(pathway)

kegg_to_graph(pathway, return_type = "visnetwork")
```

After that you can use the function map_results_to_nodes to map your differential expression results to the nodes of the graph. You can provide either a single data frame or a list of data frames containing your differential expression results. Each data frame should have a column for the feature IDs (e.g., ENTREZID) and a column for the values you want to map (e.g., logFC or log2FoldChange).

```R
library(KEGGemUP)

pathway <- "hsa04110"  # Example pathway ID

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

graph <- kegg_to_graph(pathway)
graph <- map_results_to_nodes(graph, de_results_list)
graph

```
If using a single data frame you can do it like this:

```R
library(KEGGemUP)

pathway <- "hsa04110"  # Example pathway ID

graph <- kegg_to_graph(pathway)
graph <- map_results_to_nodes(graph, res_macrophage_IFNg_vs_naive_limma)
graph

```
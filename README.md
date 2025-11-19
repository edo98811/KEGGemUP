# KEGGemUP

A package to map you differential expression results on KEGG pathways

## Installation

```R
install.packages("devtools")     
devtools::install_github("edo98811/KEGGemUP")
```
## Usage

To use the Kegg pathway visualization function:

To use it with a single differential expression results table: 
Make sure there are the two dfault columns: "KEGG_ids" for the kegg ids and "log2FoldChange" for the values to map.

```R
library(KEGGemUP)

pathway <- "hsa04110"  # Example pathway ID
graph <- kegg_to_graph(pathway)
graph <- map_results_to_nodes(graph, de_results_table)
graph
```

or 

```R
library(KEGGemUP)

pathway <- "hsa04110"  # Example pathway ID
de_results_list <-list(
  trans_limma = list(
    de_table = data.frame(res_macrophage_IFNg_vs_naive_limma[]),
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
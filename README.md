<!-- badges: start -->
[![R-CMD-check](https://github.com/edo98811/KEGGemUP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/edo98811/KEGGemUP/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

  
# KEGGemUP


## Installation

`KEGGemUP` can be installed from Bioconductor with the following code:

```
if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install("KEGGemUP")
```

You can also install the development version of `KEGGemUP` from GitHub with:

```
# install.packages("remotes")
remotes::install_github("edo98811/KEGGemUP")
```

Load the package after installation with

```
library("KEGGemUP")
```

## KEGGemUP at a glance

The `KEGGemUP` package allows you to:

* Retrieve and create a KEGG pathway graph, from the KGML files (with `create_kegg_graph()`)

* Map some continuous values onto that graph (e.g. the logFoldChange, with the `map_results_to_graph()`)

* Render that graph interactively (via `render_kegg_graph()`, based on `visNetwork`)

* Focus either on a subset or on a highlighted portion of that graph (thanks to `subset_kegg_graph()` and `highlight_kegg_graph()`)

## Quick start

```
## load an example dataset, and format the DE results
data(res_de_macro_IFNg_vs_naive, package = "KEGGemUP")

head(res_de_macro_IFNg_vs_naive)

de_results_list <- list(
  rnaseq_limma = list(
    de_table = data.frame(res_de_macro_IFNg_vs_naive),
    value_column = "logFC",
    feature_column = "ENTREZID"
  )
)

## retrieve and create the graph
kmu_cellcycle <- create_kegg_graph(pathway_id = "hsa04110")

## map the DE values to the graph
kmu_cellcycle_mapped <- map_results_to_graph(g = kmu_cellcycle, 
                                             de_results = de_results_list)

## render the graph interactively
kmu_rendered <- render_kegg_graph(g = kmu_cellcycle_mapped)

kmu_rendered
```

You can see more examples and a detailed usage in the package vignette




# Map differential expression results to nodes

Map differential expression results to nodes

## Usage

``` r
map_results_to_graph(
  g,
  de_results,
  return_type = "visNetwork",
  feature_column = NULL,
  value_column = NULL,
  palette = "RdBu"
)
```

## Arguments

- g:

  An igraph object representing the pathway.

- de_results:

  Named list of differential expression results.

- return_type:

  Output type: 'igraph' or 'visNetwork'.

- feature_column:

  Column name in de_table containing KEGG IDs (if de_results is a single
  data.frame).

- value_column:

  Column name in de_table containing values to map (if de_results is a
  single data.frame).

- palette:

  Color palette for node coloring (default: "RdBu").

## Value

An igraph or visNetwork object with mapped results.

## Details

This functionmaps differential expression results onto the nodes of a
KEGG pathway graph. The pathwhay given as input must be the output of
the function `kegg_to_graph`.

This function can be used to map the differential expression results to
the graph, the input of the graph must be the output of the function
`kegg_to_graph` in the igraph format. The results to be mapped can be
provided either as a list or as a single data.frame. If a single
data.frame is provided, the default column names for KEGG IDs and values
are 'KEGG_ids' and 'log2FoldChange', respectively, but these can be
changed using the `feature_column` and `value_column` parameters.

## Examples

``` r
pathway <- "hsa04110" # Example pathway ID
graph <- kegg_to_graph(pathway, return_type = "igraph")
#> Downloaded & cached: hsa04110
#> Parsed 134 nodes from KGML file.
#> Parsed 119 relationship edges from KGML file.
#> Parsed 0 reaction edges from KGML file.
#> Total edges parsed from KGML file: 119
# Example differential expression results
de_results <- data.frame(
  KEGG_ids = c("hsa:1234", "hsa:5678", "cpd:C00022"),
  log2FoldChange = c(1.5, -2.0, 0.5)
)
vis_graph <- map_results_to_graph(graph, de_results, return_type = "visNetwork")
#> Mapping differential expression results to nodes...
#> de_results provided as a single data.frame. Using default column names: value_column = 'log2FoldChange', feature_column = 'KEGG_ids'.
```

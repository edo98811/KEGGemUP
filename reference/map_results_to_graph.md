# Map continuous values to graph nodes

Map differential expression results to nodes

## Usage

``` r
map_results_to_graph(
  g,
  de_results,
  feature_column = NULL,
  value_column = NULL,
  verbose = FALSE,
  palette = NULL,
  palette_limit = NULL,
  palettes_list = list(NA_character_),
  palettes_limits_list = c(NA_real_)
)
```

## Arguments

- g:

  An igraph object representing the KEGG pathway graph.

- de_results:

  A data frame or a list of data frames containing differential
  expression results

- feature_column:

  Name of the column in de_results that contains KEGG IDs

- value_column:

  Name of the column in de_results that contains values to map

- verbose:

  Logical indicating whether to print verbose messages (default: FALSE)

- palette:

  Optional color palette for mapping values (default: NULL, will use a
  default palette)

- palette_limit:

  Optional numeric limit for the color palette (default: NULL, will be
  determined from data)

- palettes_list:

  Optional list of color palettes if de_results is a list (default:
  NULL)

- palettes_limits_list:

  Optional list of numeric limits for multiple palettes if de_results is
  a list (default: NULL)

## Value

An igraph object with differential expression results mapped to node
attributes

## Details

This function can be used to map the differential expression results to
the graph, the input of the graph must be the output of the function
`create_kegg_graph` in the igraph format. The results to be mapped can
be provided either as a list or as a single data.frame. If a single
data.frame is provided, the default column names that it will look for
are KEGG IDs and values are 'KEGG_ids' and 'log2FoldChange',
respectively, but these can be changed using the `feature_column` and
`value_column` parameters.

## Examples

``` r
pathway <- "hsa04110" # Example pathway ID
graph <- create_kegg_graph(pathway_id = pathway)
# Example differential expression results
de_results <- data.frame(
  KEGG_ids = c("hsa:1234", "hsa:5678", "cpd:C00022"),
  log2FoldChange = c(1.5, -2.0, 0.5)
)
vis_graph <- map_results_to_graph(
  graph,
  de_results,
  feature_column = "KEGG_ids",
  value_column = "log2FoldChange"
)
#> de_results provided as a single data.frame. Using provided value_column: 'log2FoldChange' and feature_column: 'KEGG_ids'.
```

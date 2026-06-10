# Create an interactive visualization of KEGG pathways with visNetwork

Plot KEGG pathway graph using visNetwork

## Usage

``` r
render_kegg_graph(
  g,
  graph_title = NULL,
  scaling_factor = 1.5,
  relationships = c("all", "reactions", "relations", "none"),
  visualization_type = c("standard", "positions", "node_name", "node_size")
)
```

## Arguments

- g:

  An igraph object representing the KEGG pathway graph.

- graph_title:

  Character string, used as a title for the rendered graph. Defaults to
  NULL, which would fall back to the name specified in the title
  attribute of the graph.

- scaling_factor:

  Numeric factor to scale node sizes (default: 1.5).

- relationships:

  Character specifying which relationships to include in edges ("all",
  "reactions", "relations", "none"; default: "all").

- visualization_type:

  Character specifying the type of visualization for nodes: "standard",
  "positions", "node_name", or "node_size" (default: "standard").

## Value

A visNetwork object representing the KEGG pathway graph with mapped
results.

## Details

This function takes an igraph object representing a KEGG pathway graph
and creates a visNetwork visualization. It maps node attributes to
visual properties.

## Examples

``` r
pathway <- "hsa04110" # Example pathway ID
graph <- create_kegg_graph(pathway_id = pathway)
# Example differential expression results
de_results <- data.frame(
  KEGG_ids = c("hsa:1234", "hsa:5678", "cpd:C00022"),
  log2FoldChange = c(1.5, -2.0, 0.5)
)
graph <- map_results_to_graph(
  graph,
  de_results,
  feature_column = "KEGG_ids",
  value_column = "log2FoldChange")
#> de_results provided as a single data.frame. Using provided value_column: 'log2FoldChange' and feature_column: 'KEGG_ids'.

vis_graph <- render_kegg_graph(graph, scaling_factor = 1.5,
relationships = "all", visualization_type = "standard")
```

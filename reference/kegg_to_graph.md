# Transform a ggkegg graph to igraph or visNetwork

Transform a ggkegg graph to igraph or visNetwork

## Usage

``` r
kegg_to_graph(pathway_id, return_type = "igraph", scaling_factor = 1.5)
```

## Arguments

- pathway_id:

  KEGG pathway ID (e.g., 'hsa:04110' or '04110').

- return_type:

  Output type: 'igraph' or 'visNetwork'.

- scaling_factor:

  Numeric factor to scale node sizes.

## Value

An igraph or visNetwork object representing the pathway.

## Details

This function downloads the KGML file for the specified KEGG pathway,
then parses it to generate a graph representation using either the
igraph or visNetwork package. It styles nodes and edges based on their
types the output can be used for visualization or further analysis. If
differential expression results are provided, they can be mapped to the
nodes using the function `map_results_to_graph`.

## Examples

``` r
pathway <- "hsa04110" # Example pathway ID
graph <- kegg_to_graph(pathway)
#> adding rname 'https://rest.kegg.jp/get/hsa04110/kgml'
#> 
#> Downloaded & cached: hsa04110
#> Parsed 134 nodes from KGML file.
#> Parsed 119 relationship edges from KGML file.
#> Warning: No entries found in kgml file.
#> Parsed 0 reaction edges from KGML file.
#> Total edges parsed from KGML file: 119
kegg_to_graph(pathway, return_type = "visNetwork")
#> Downloaded & cached: hsa04110
#> Parsed 134 nodes from KGML file.
#> Parsed 119 relationship edges from KGML file.
#> Warning: No entries found in kgml file.
#> Parsed 0 reaction edges from KGML file.
#> Total edges parsed from KGML file: 119
#> Error in visNetwork::visInteraction(v, dragNodes = TRUE, multiselect = TRUE,     selectable = TRUE) %>% visNetwork::visEvents(selectNode = "function(nodes) {\n        Shiny.setInputValue('graph_click', nodes.nodes, {priority: 'event'});\n      }",     deselectNode = "function(nodes) {\n        Shiny.setInputValue('graph_click', nodes.nodes, {priority: 'event'});\n      }"): could not find function "%>%"
```

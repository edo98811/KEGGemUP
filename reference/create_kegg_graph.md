# Create a graph for a KEGG pathway

Convert KEGG pathway to igraph object

## Usage

``` r
create_kegg_graph(pathway_id, kgml_file = NULL, verbose = FALSE)
```

## Arguments

- pathway_id:

  Character, KEGG pathway ID (e.g., 'hsa04110') or NULL if using
  kgml_file

- kgml_file:

  Path to a local KGML file (optional, if pathway_id is provided)

- verbose:

  Logical indicating whether to print verbose messages (default: FALSE)

## Value

An igraph object representing the KEGG pathway graph

## Examples

``` r
pathway <- "hsa04110" # Example pathway ID
graph <- create_kegg_graph(pathway_id = pathway, verbose = TRUE)
#> Downloading KGML file for pathway ID: hsa04110
#> Downloading KGML from: https://rest.kegg.jp/get/hsa04110/kgml
#> Cached: hsa04110
plot(graph) # Plot the graph using igraph's plotting functions
```

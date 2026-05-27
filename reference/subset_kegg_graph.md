# Subset a graph object for a KEGG pathway

Create igraph visualization with improved layout

## Usage

``` r
subset_kegg_graph(g, ids_to_include)
```

## Arguments

- g:

  An igraph object to visualize. Must have vertex attributes 'x' and 'y'
  for layout.

- ids_to_include:

  Character vector of KEGG IDs to include in the subset graph.

## Value

A plot of the igraph object with improved layout.

## Details

All the edges between the vertices are plotted automatically.

## Examples

``` r
pathway <- "mmu00230"
g <- create_kegg_graph(pathway)
KEGG_to_include <- c("C00262", "C00385", "C00366", "C00294", "C00387",
                 "C01762", "C05512", "C00301", "C01185", "C00455",
                 "22436", "14544", "18950", "11486", "80285", "59027")
subg <- subset_kegg_graph(g, KEGG_to_include)
# plot(g)
# plot(subg)
```

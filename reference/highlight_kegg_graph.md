# Highlight a subset of the KEGG graph

Highlight a subset of the graph based on KEGG IDs

## Usage

``` r
highlight_kegg_graph(g, ids_to_highlight)
```

## Arguments

- g:

  An igraph object to visualize. Must have vertex attributes 'x' and 'y'
  for layout.

- ids_to_highlight:

  Character vector of KEGG IDs to include in the highlighted subset.

## Value

An igraph object with highlighted nodes and faded non-highlighted nodes
and edges.

## Details

This function highlights the nodes corresponding to the provided KEGG
IDs and fades the rest of the graph. It modifies vertex attributes to
achieve this effect.

## Examples

``` r
pathway <- "mmu00230"
g <- create_kegg_graph(pathway)
#> adding rname 'https://rest.kegg.jp/get/mmu00230/kgml'
#> 
KEGG_to_include <- c("C00262", "C00385", "C00366", "C00294", "C00387",
                 "C01762", "C05512", "C00301", "C01185", "C00455",
                 "22436", "14544", "18950", "11486", "80285", "59027")
highlighted_subg <- highlight_kegg_graph(g, KEGG_to_include)
```

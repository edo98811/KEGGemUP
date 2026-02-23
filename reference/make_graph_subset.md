# Create igraph visualization with improved layout

Create igraph visualization with improved layout

## Usage

``` r
make_graph_subset(g, ids_to_include)
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

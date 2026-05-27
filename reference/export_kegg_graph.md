# Export a KEGG graph

Exporting a KEGG graph into its components, nodes and edges, as
tab-separated text files (having them represented as dataframes for max
portability)

## Usage

``` r
export_kegg_graph(g, basename)
```

## Arguments

- g:

  An igraph graph object, e.g. created with KEGGemUP

- basename:

  Character string, specifying the base name for the files to write the
  two individual data frames, for nodes and edges

## Value

NULL, invisibly

## Examples

``` r

g <- create_kegg_graph(pathway_id = "hsa04110")
export_kegg_graph(g, basename = tempfile())
#> Exported graph components in /tmp/Rtmp4cNBm3/file1e143d8e5151_nodes.tsv and /tmp/Rtmp4cNBm3/file1e143d8e5151_edges.tsv
```

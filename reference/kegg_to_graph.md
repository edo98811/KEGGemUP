# Convert KEGG pathway to graph

Convert KEGG pathway to graph

## Usage

``` r
kegg_to_graph(
  pathway_id,
  kgml_file = NULL,
  scaling_factor = 3,
  verbose = FALSE
)
```

## Arguments

- pathway_id:

  KEGG pathway ID (e.g., "hsa04110").

- kgml_file:

  Optional local path to a KGML file. If provided, the function will use
  this file instead of downloading it.

- scaling_factor:

  Numeric scaling factor for node sizes (default: 2).

- verbose:

  Logical, if TRUE, print additional messages.

## Value

An igraph object representing the KEGG pathway graph.

## Details

This function downloads the KGML file for a given KEGG pathway ID,
parses it and constructs an igraph object.

## Examples

``` r
pathway_graph <- kegg_to_graph(
  pathway_id = "hsa04110"
)
#> adding rname 'https://rest.kegg.jp/get/hsa04110/kgml'
#> 
#> adding rname 'https://rest.kegg.jp/list/compound'
#> 
#> adding rname 'https://rest.kegg.jp/list/glycan'
#> 
#> adding rname 'https://rest.kegg.jp/list/ko'
#> 
#> adding rname 'https://rest.kegg.jp/list/enzyme'
#> 
#> adding rname 'https://rest.kegg.jp/list/reaction'
#> 
```

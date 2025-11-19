# Map differential expression results to nodes

Map differential expression results to nodes

## Usage

``` r
map_results_to_nodes(g, de_results, return_type = "visNetwork")
```

## Arguments

- g:

  An igraph object representing the pathway.

- de_results:

  Named list of differential expression results.

- return_type:

  Output type: "igraph" or "visNetwork".

## Value

An igraph or visNetwork object with mapped results.

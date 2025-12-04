# Map differential expression results to nodes

Map differential expression results to nodes

## Usage

``` r
map_results_to_nodes(
  g,
  de_results,
  return_type = "visNetwork",
  feature_column = NULL,
  value_column = NULL
)
```

## Arguments

- g:

  An igraph object representing the pathway.

- de_results:

  Named list of differential expression results.

- return_type:

  Output type: "igraph" or "visNetwork".

- feature_column:

  Column name in de_table containing KEGG IDs (if de_results is a single
  data.frame).

- value_column:

  Column name in de_table containing values to map (if de_results is a
  single data.frame).

## Value

An igraph or visNetwork object with mapped results.

# Transform a ggkegg graph to igraph or visNetwork

Transform a ggkegg graph to igraph or visNetwork

## Usage

``` r
kegg_to_graph(
  path_id,
  de_results = NULL,
  return_type = "igraph",
  scaling_factor = 1.5
)
```

## Arguments

- path_id:

  KEGG pathway ID (e.g., "hsa:04110" or "04110").

- de_results:

  Named list of differential expression results. Each entry should be a
  list with elements: de_table (data.frame), value_column (character),
  feature_column (character), threshold (numeric).

- return_type:

  Output type: "igraph" or "visNetwork".

- scaling_factor:

  Numeric factor to scale node sizes.

## Value

An igraph or visNetwork object representing the pathway.

# Add results from combined results data frame to nodes data frame

Add results from combined results data frame to nodes data frame

## Usage

``` r
add_results_nodes(nodes_df, results_combined)
```

## Arguments

- nodes_df:

  Data frame of nodes with a column 'KEGG'

- results_combined:

  Data frame with combined results containing columns: KEGG, value,
  source

## Value

Updated nodes data frame with added columns: value, color, source, text

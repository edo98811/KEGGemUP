# Style nodes based on their type for visNetwork visualization.

Style nodes based on their type for visNetwork visualization.

## Usage

``` r
style_nodes(nodes_df, node_size_multiplier = 1.2)
```

## Arguments

- nodes_df:

  Data frame of nodes with a column 'type'.

- node_size_multiplier:

  Numeric factor to scale node sizes (default: 1.2).

## Value

nodes_df with added visual styling columns: shape, fixed,
widthConstraint, heightConstraint, size.

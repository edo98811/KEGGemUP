# Add compound names to compound nodes in the nodes data frame.

Add compound names to compound nodes in the nodes data frame.

## Usage

``` r
add_compound_names(nodes_df, bfc)
```

## Arguments

- nodes_df:

  Data frame of nodes with a column 'type' indicating node type.

- bfc:

  BiocFileCache object for caching KEGG compound mappings.

## Value

Updated nodes data frame with compound names added to compound nodes.

# Add tooltips to nodes for visNetwork visualization.

Add tooltips to nodes for visNetwork visualization.

## Usage

``` r
add_tooltip(nodes_df)
```

## Arguments

- nodes_df:

  Data frame of nodes with columns: KEGG, label, source, value.

## Value

nodes_df with added 'title' column for tooltips.

## Details

The tooltip includes a button to the specific KEGG entry page. If
multiple KEGG IDs are present, they are concatenated with '+' in the
URL. It also adds information about the node name, source of
differential expression data, and value.

# Handle multiple KEGG IDs in a single string separated by ";"

Handle multiple KEGG IDs in a single string separated by ";"

## Usage

``` r
expand_keggs(kegg_df)
```

## Arguments

- kegg_df:

  Data frame with at least two columns: 'name' and 'KEGG'

## Value

Data frame with a single column 'KEGG' containing individual KEGG IDs

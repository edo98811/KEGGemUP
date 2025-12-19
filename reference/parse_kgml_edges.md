# Parse KEGG KGML files to extract combined edges data frame.

Parse KEGG KGML files to extract combined edges data frame.

## Usage

``` r
parse_kgml_edges(file)
```

## Arguments

- file:

  Path to the KGML XML file.

## Value

A data.frame combining both relations and reactions edges.

## Details

This function combines the outputs of `parse_kgml_relations` and
`parse_kgml_reactions` to provide a comprehensive edges data frame
representing all interactions defined in the KGML file.

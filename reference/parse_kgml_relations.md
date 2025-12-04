# Parse KEGG KGML files to extract relations and edges data frames.

Parse KEGG KGML files to extract relations and edges data frames.

## Usage

``` r
parse_kgml_relations(file)
```

## Arguments

- file:

  Path to the KGML XML file.

## Value

A data.frame with the following columns:

- from:

  The ID of the source entry (node) in the pathway, corresponding to the
  `entry1` attribute in the KGML `<relation>` tag.

- to:

  The ID of the target entry (node) in the pathway, corresponding to the
  `entry2` attribute in the KGML `<relation>` tag.

- type:

  The general type of relationship between the two entries (e.g.,
  `"ECrel"`, `"PPrel"`, `"GErel"`, `"PCrel"`, `"maplink"`). These types
  describe the biological nature of the connection, such as
  enzyme-enzyme relation or protein-protein interaction.

- subtype:

  A more specific subtype of the relation, derived from the `<subtype>`
  child elements of the `<relation>` node (e.g., `"activation"`,
  `"inhibition"`, `"expression"`, `"compound"`). If no subtype is
  defined, this will be `NA`.

- rel_value:

  A categorical value associated with the subtype, which encodes the
  apearence of the arrow (from the `value` attribute of `<subtype>`). If
  not present, this will be `NA`.

## Details

The function reads a KEGG KGML (KEGG Markup Language) file, which
encodes pathway information as XML, and extracts all `<relation>`
elements that describe the interactions or relationships between
entities in the pathway. Each `<relation>` may contain one or more
`<subtype>` elements that provide additional details about the
interaction. The result is a tidy data.frame suitable for network
analysis or visualization, where each row represents one
relation–subtype pair.

## Examples

``` r
if (FALSE) { # \dontrun{
relations_df <- parse_kgml_relations("pathway.xml")
head(relations_df)
} # }
```

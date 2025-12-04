# Parse KEGG KGML files to extract nodes data frame.

Parse KEGG KGML files to extract nodes data frame.

## Usage

``` r
parse_kgml_entries(file)
```

## Arguments

- file:

  Path to the KGML XML file.

## Value

A data.frame with the following columns:

- id:

  Unique identifier of the entry within the KGML pathway (from the `id`
  attribute).

- kegg_name:

  The KEGG-specific name or identifier of the entity (from the `name`
  attribute). This may include one or more KEGG identifiers such as gene
  IDs, compound IDs, or enzyme EC numbers. Preceded by organism ID

- type:

  Type of the node (from the `type` attribute), indicating the
  biological entity class such as `"gene"`, `"enzyme"`, `"compound"`,
  `"map"`, `"ortholog"`, or `"group"`.

- link:

  URL linking to the KEGG resource for this entry, if available (from
  the `link` attribute).

- reaction:

  Associated reaction ID(s), if any (from the `reaction` attribute).
  Typically present for enzyme entries.

- graphics_name:

  Display name for the entry, taken from the `name` attribute of the
  `<graphics>` node.

- label:

  Text label for visualization purposes.

- fgcolor:

  Foreground color of the graphical element (from the `fgcolor`
  attribute of `<graphics>`).

- bgcolor:

  Background color of the graphical element (from the `bgcolor`
  attribute of `<graphics>`).

- graphics_type:

  Shape or representation type of the graphical element (from the `type`
  attribute of `<graphics>`), such as `"rectangle"`, `"circle"`, or
  `"line"`.

- x:

  X-coordinate of the node’s position in the pathway diagram (from the
  `x` attribute of `<graphics>`).

- y:

  Y-coordinate of the node’s position in the pathway diagram (from the
  `y` attribute of `<graphics>`).

- width:

  Width of the graphical element (from the `width` attribute of
  `<graphics>`).

- height:

  Height of the graphical element (from the `height` attribute of
  `<graphics>`).

## Details

The function parses a KEGG KGML (KEGG Markup Language) XML file and
extracts all `<entry>` elements, each representing a biological entity
in a KEGG pathway diagram. Each node may contain nested `<graphics>`
elements defining visual properties (such as position, size, and colors)
and `<component>` elements that define group membership for composite
entities.

The resulting data provides a tidy, one-row-per-entry representation
suitable for integration with relational data models or network
visualization frameworks (e.g., `igraph`or `visNetwork`).

## Examples

``` r
if (FALSE) { # \dontrun{
nodes_df <- parse_kgml_entries("pathway.xml")
head(nodes_df)
} # }
```

# The \`KEGGemUP\` User's Guide

**Compiled date**: 2026-07-14

**Last edited**: 2026-04-24

**License**: MIT + file LICENSE

## Introduction

Pathway-based analysis is a central strategy in biomedical research,
enabling the interpretation of datasets stemming from different omics
techniques, by linking quantitative high-dimensional profiles of
molecular features changes to higher-level biological mechanisms. The
Kyoto Encyclopedia of Genes and Genomes (KEGG, (Kanehisa 2000)) is one
of the most widely used pathway resources, also thanks to the support
offered for a large variety of species, and including different
biological entities (genes, proteins, metabolites). However, the
integration of KEGG pathways into omics workflows is currently either
relying on manual operations, or limited to generating static
visualizations which can be difficult to explore in further detail.

We addressed this problem by leveraging interactive visualization
capabilities of R, coupled to a streamlined process to retrieve the
pathway information and contextualize any omics data at hand. The
`KEGGemUP` R package that can parse the KGML files provided by the KEGG
API and builds an igraph representation of the pathway (Antonov et al.
2023), on which different data and results can be seamlessly mapped. The
resulting pathways can be then rendered and efficiently explorer in
interactive widgets that also support further drilldown operations
(e.g. by linking to existing databases), with tooltips and bindings
compatible within Shiny.

We believe that this can simplify and speed up the iterative process of
interpretation of omics data (provided as simple data formats widely
adopted within Bioconductor), either for single features or within the
frame of KEGG pathways. Leveraging a caching mechanism for efficient
retrieval of KGML files (with `BiocFileCache`), the `igraph` framework
to represent the information-rich extracted pathways, and a powerful
HTML-based network visualization system (via `visNetwork`) , `KEGGemUP`
delivers a modern, **interactive** way to **create**, **visualize**, and
**explore** KEGG pathways from the point of view of integrated omics
layers.

## Installation

`KEGGemUP` can be installed from Bioconductor with the following code:

``` r

if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install("KEGGemUP")
```

Load the package after installation with

``` r

library("KEGGemUP")
```

### Example dataset

We will demonstrate the features of `KEGGemUP` on a small dataset,
extracted from the differential expression analysis workflow, running
the limma-voom pipeline on the
*[macrophage](https://bioconductor.org/packages/3.23/macrophage)*
dataset. Specifically, we have compared here the Interferon gamma
treated samples vs the naive ones, while accounting for the cell line of
origin - we refer the reader to the script in the `inst/scripts` folder
to see how that object has been generated.

``` r

data(res_de_macro_IFNg_vs_naive, package = "KEGGemUP")

## Inspect this briefly
head(res_de_macro_IFNg_vs_naive)
#>                    logFC     CI.L     CI.R  AveExpr        t      P.Value
#> ENSG00000125347 5.616091 5.189756 6.042426 8.371198 27.61754 1.603356e-16
#> ENSG00000137496 3.897544 3.557251 4.237838 6.320924 24.01261 2.004938e-15
#> ENSG00000111181 5.055709 4.619754 5.491665 2.794973 24.31318 1.602755e-15
#> ENSG00000244731 5.626104 5.129536 6.122672 2.837925 23.75369 2.436608e-15
#> ENSG00000162645 6.734081 6.121447 7.346715 8.131362 23.04512 4.197967e-15
#> ENSG00000204257 4.071624 3.687518 4.455729 4.328773 22.22384 8.045146e-15
#>                    adj.P.Val        B ENTREZID gene_name
#> ENSG00000125347 2.649385e-12 27.95382     3659      IRF1
#> ENSG00000137496 1.006563e-11 25.54771    10068    IL18BP
#> ENSG00000111181 1.006563e-11 25.52206     6539   SLC6A12
#> ENSG00000244731 1.006563e-11 24.97471      720       C4A
#> ENSG00000162645 1.387344e-11 24.79998     2634      GBP2
#> ENSG00000204257 2.215633e-11 24.15399     3108   HLA-DMA
```

As you can see, this is simply a `data.frame` with the table of the
top-ranked genes from a linear model fit. If using the
*[DESea2](https://bioconductor.org/packages/3.23/DESea2)* framework, you
can obtain something similar by running
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on the
output of the \`results()\`\` function.

In particular, we will focus on the regulation occurring between these
two conditions, and we will use the logFC continuous value as a measure
of the effect size.

## KEGGemUP at a glance

The `KEGGemUP` package allows you to:

- Retrieve and create a KEGG pathway graph, from the KGML files (with
  [`create_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/create_kegg_graph.md))

- Map some continuous values onto that graph (e.g. the logFoldChange,
  with the
  [`map_results_to_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/map_results_to_graph.md))

- Render that graph interactively (via
  [`render_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/render_kegg_graph.md),
  based on `visNetwork`)

- Focus either on a subset or on a highlighted portion of that graph
  (thanks to
  [`subset_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/subset_kegg_graph.md)
  and
  [`highlight_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/highlight_kegg_graph.md))

The remainder of this vignette will cover how to perform these
operations to explore omics data in the context of KEGG pathway
diagrams.

## Building KEGG pathway graphs

KEGG database uses pathway identifiers to refer to each pathway - for
example, `hsa00563` is the KEGG pathway ID for
“Glycosylphosphatidylinositol (GPI)-anchor biosynthesis” in humans.

`KEGGemUP` retrieves for you (supported by a caching mechanism, by
default) the KGML files, parses them, and builds information-rich
`igraph` graph objects.

The identifier must contain the organism prefix, such as “hsa” for human
pathways, “mmu” for mouse pathways, etc. - Please refer to the [KEGG
website](https://www.kegg.jp/) to access the main components of the
database.

``` r

kmu_graph <- create_kegg_graph(pathway_id = "hsa00563")
#> adding rname 'https://rest.kegg.jp/get/hsa00563/kgml'

kmu_graph
#> IGRAPH eedddd7 DN-B 147 384 -- 
#> + attr: title (g/c), type (g/c), name (v/c), type (v/c), link (v/c),
#> | reaction (v/c), reaction_link (v/c), reaction_label (v/c),
#> | graphics_name (v/c), node_name (v/c), KEGG (v/c), components (v/c),
#> | line_id (v/c), coords (v/c), graphics_type (v/c), point_index (v/n),
#> | de_value (v/n), de_text (v/c), de_source (v/c), de_name (v/c),
#> | ids_for_mapping (v/c), label (v/c), x (v/n), y (v/n), width (v/n),
#> | height (v/n), source (v/c), fgcolor (v/c), bgcolor (v/c), group
#> | (v/c), vertex.color (v/c), size (v/n), fixed (v/l), shape (v/c), type
#> | (e/c), relation_type (e/c), relation_subtype_name (e/c),
#> | relation_subtype_value (e/c), reaction_alt_name_substrate (e/c),
#> | reaction_alt_name_product (e/c), reaction_id (e/c), reaction_name
#> | (e/c), reaction_type (e/c), reaction_from_name (e/c),
#> | reaction_to_name (e/c), point_index (e/n), directed (e/l), color
#> | (e/c), width (e/n), value (e/n), label (e/c), lty (e/n), arrow.mode
#> | (e/c), dashes (e/l), title (e/c), font.face (e/c), borderWidth (e/n),
#> | borderWidthSelected (e/n), from (e/c), to (e/c)
#> + edges from eedddd7 (vertex names):
```

While one can directly use the `plot` routine to visualize this graph,
most of its attributes are set to behave correctly in the interactive
viewer provided by
*[visNetwork](https://CRAN.R-project.org/package=visNetwork)*.
Nonetheless, using `igraph` objects as a container makes it easy to
convert/reuse by other graph representations.

The optional `kgml_file` parameter provides the option to specify the
path to a local KGML file, while `verbose` can print some messages
during the execution.

## Rendering KEGG pathway graphs

The advantage of having parsed all elements of a KGML file into an
`igraph` object directly tailored to be used within `visNetwork` is that
we can simply create an interactive view of that pathway with the
[`render_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/render_kegg_graph.md)
command.

``` r

render_kegg_graph(g = kmu_graph)
```

This is pretty much a vanilla representation of the KEGG pathway, but it
already contains some tooltips, displayed upon hovering the mouse on the
nodes, that contain buttons linking to the respective KEGG database
entries.

Other parameters such as `scaling_factor`, `relationships` and
`visualization_type` define the aspect of the rendered graph. In most
situations, these can be left as default values, which tends to generate
a truthful representation of the KEGG pathway graph.

In the Viewer pane of an IDE such as RStudio, it is possible to interact
with the graph with operations such as zoom, pan, select, hover -
instead of having a static view, this might become extremely useful to
obtain a deeper understanding.

Optionally, a pathway graph can also be rendered without the node
reporting the title itself as a vertex. This can be achieved combining
the
[`cleanup_title_node()`](https://imbeimainz.github.io/KEGGemUP/reference/cleanup_title_node.md)
function before rendering.

``` r

render_kegg_graph(g = cleanup_title_node(kmu_graph))
```

## Mapping values onto KEGG pathway graphs

To map continuous values (such as the logFoldChange obtained from
differential expression results) to the nodes of a KEGG pathway graph,
you can use
[`map_results_to_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/map_results_to_graph.md).

Its main parameter, `g`, is the `igraph` object returned by
[`create_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/create_kegg_graph.md).
Notably, there are two ways to specify the values to map onto it, using
the `de_results` parameter:

- A single `data.frame` with the differential expression results. This
  dataframe must contain at least two columns: one with the KEGG feature
  IDs (without organism prefix), and another with the values to map to
  the nodes. These two column names need to be specified by the
  parameters `feature_column` and `value_column` (otherwise defaulting
  to NULL).

- A nested list, which can also handle multiple differential expression
  results tables (e.g., one for RNA-seq and one for metabolomics). Each
  element of the list can have a name, and has to be structured itself
  as a list, containing the following elements:

  - `de_table`: a data.frame with the differential expression results
    (as in the single case)
  - `value_column`: the name of the column in `de_table` containing the
    values to map to the nodes
  - `feature_column`: the name of the column in de_table containing the
    feature IDs (e.g., ENTREZ IDs) corresponding to the KEGG ids in the
    graph, again without organism prefix)

Some typical values for `value_column` might be “logFC” or
“logFoldChange”.

Usually, the KEGG identifiers without organism prefix correspond to the
ENTREZ IDs for the genes. For other feature types (e.g., compounds) you
will need to make sure that the IDs in your differential expression
results table match the KEGG IDs used in the graph.  
Note: If your data frame does not contain a compatible identifier, it
might still be straightforward to add such a column with packages such
as AnnotationDbi and the orgDb packages provided by Bioconductor.

In the case of a single result (`res_de_macro_IFNg_vs_naive`):

``` r

head(res_de_macro_IFNg_vs_naive)
#>                    logFC     CI.L     CI.R  AveExpr        t      P.Value
#> ENSG00000125347 5.616091 5.189756 6.042426 8.371198 27.61754 1.603356e-16
#> ENSG00000137496 3.897544 3.557251 4.237838 6.320924 24.01261 2.004938e-15
#> ENSG00000111181 5.055709 4.619754 5.491665 2.794973 24.31318 1.602755e-15
#> ENSG00000244731 5.626104 5.129536 6.122672 2.837925 23.75369 2.436608e-15
#> ENSG00000162645 6.734081 6.121447 7.346715 8.131362 23.04512 4.197967e-15
#> ENSG00000204257 4.071624 3.687518 4.455729 4.328773 22.22384 8.045146e-15
#>                    adj.P.Val        B ENTREZID gene_name
#> ENSG00000125347 2.649385e-12 27.95382     3659      IRF1
#> ENSG00000137496 1.006563e-11 25.54771    10068    IL18BP
#> ENSG00000111181 1.006563e-11 25.52206     6539   SLC6A12
#> ENSG00000244731 1.006563e-11 24.97471      720       C4A
#> ENSG00000162645 1.387344e-11 24.79998     2634      GBP2
#> ENSG00000204257 2.215633e-11 24.15399     3108   HLA-DMA
```

… we would use “ENTREZID” as `feature_name` and “logFC” as `value_name`.

Otherwise, one can build the list-based object with this command:

``` r

de_results_list <- list(
  rnaseq_limma = list(
    de_table = data.frame(res_de_macro_IFNg_vs_naive),
    value_column = "logFC",
    feature_column = "ENTREZID"
  )
)
```

In this case, only one element is specified (named “rnaseq_limma”), but
additional entries can be specified afterwards extending the
`de_results_list` object.

Applying this step onto the previously generated `kmu_graph`, we would
simply need to run the following commands:

``` r

kmu_graph_mapped <- map_results_to_graph(g = kmu_graph, 
                                         de_results = de_results_list)
kmu_rendered <- render_kegg_graph(g = kmu_graph_mapped, 
                                  scaling_factor = 1.3)
kmu_rendered
```

As you can see,
[`map_results_to_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/map_results_to_graph.md)
utilizes by default a red-to-blue palette for the mapping of the values
(`RdBu` from the RColorBrewer package), detecting the range directly
from the data. This behavior can be controlled via the `palette` and
`palette_limit` parameters (or `palettes_list` and
`palettes_limits_list` if using multiple palettes for multiple entries
of `de_results`).

For example, using the single `res_de_macro_IFNg_vs_naive` DE result,
one could call:

``` r

kmu_graph_mapped_single <- map_results_to_graph(g = kmu_graph, 
                                                de_results = res_de_macro_IFNg_vs_naive,
                                                feature_column = "ENTREZID",
                                                value_column = "logFC", 
                                                palette = "Spectral")
#> de_results provided as a single data.frame. Using provided value_column: 'logFC' and feature_column: 'ENTREZID'.
render_kegg_graph(g = kmu_graph_mapped_single, scaling_factor = 1.3)
```

… to obtain an equivalent view, using the `Spectral` palette. Any
palette supported by RColorBrewer can be used, see the output of
[`RColorBrewer::display.brewer.all()`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)
for a complete reference.

## Focusing on subsets of a pathway graph

Sometimes one might be interested in focusing on a specific subset of
the graph, as pathway graphs can be almost too comprehensive (especially
in the case of larger pathways).

`KEGGemUP` provides two approaches to do so, with the
[`subset_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/subset_kegg_graph.md)
and the
[`highlight_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/highlight_kegg_graph.md)
functions. Their interface is similar, feeding on the `g` graph object,
created in the first steps via
[`create_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/create_kegg_graph.md),
and on a vector of nodes to keep or highlight. Both these functions
return a modified `igraph` object, that again can conveniently be
rendered by
[`render_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/render_kegg_graph.md) -
the chunks below illustrate these use cases.

Let’s work with the Cell cycle pathway, for human - “hsa04110” in the
KEGG database:

``` r

kmu_cellcycle <- create_kegg_graph(pathway_id = "hsa04110")
kmu_cellcycle_mapped <- map_results_to_graph(g = kmu_cellcycle, 
                                             de_results = de_results_list)
#> Warning in add_results_nodes(vertices_df, results_combined, verbose = verbose):
#> Some nodes had multiple matching IDs; only the first match was used for
#> de_value/de_source.

## Selecting a subset of the nodes to keep
KEGGids_to_include <- igraph::V(kmu_cellcycle_mapped)$ids_for_mapping |>
  tail(35) |> 
  strsplit(split = ";") |>
  unlist() |>
  unique()
head(KEGGids_to_include)
#> [1] "5925"   "5933"   "5934"   "2810"   "151648" "6500"
length(KEGGids_to_include)
#> [1] 38

kmu_cellcycle_subset <- subset_kegg_graph(g = kmu_cellcycle_mapped, 
                                          ids_to_include = KEGGids_to_include)

## Specifying a custom title with the `graph_title` parameter...
render_kegg_graph(kmu_cellcycle_subset, graph_title = "The subset of the original graph")
```

By keeping the nodes, also the edges that exist and connect them are
retained.

As you have seen, the names to pick the node from are coming from the
pool included in the `ids_for_mapping` attribute of the graph nodes, and
in our example we are simply keeping the last 35 (corresponding to 38
identifiers).

Some intermediate operations might be needed to convert from the more
classical gene symbols to the IDs used within KEGG (that often do
correspond to ENTREZ IDs).

Another step could be to select a subset of nodes to highlight - or in
other words, “grey out” the nodes that are not of interest (and with
them, also the edges that are not involving the selected nodes at both
ends). Reusing the object created above, you can run the following chunk
to create the “highlighted version”:

``` r

kmu_cellcycle_highlighted <- highlight_kegg_graph(g = kmu_cellcycle_mapped, 
                                                  ids_to_highlight = KEGGids_to_include)
render_kegg_graph(kmu_cellcycle_highlighted)
```

## Working further with KEGG graphs

Of course, the graphs created with `KEGGemUP` can also be processed
further with alternative packages (if these are compatible with `igraph`
objects). To ensure interoperability, the
[`export_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/export_kegg_graph.md)
function provides a compact wrapper to create two tab-separated text
files, where all the info on nodes and edges are written.

``` r

file_prefix <- tempfile()
export_kegg_graph(g = kmu_cellcycle_highlighted, basename = file_prefix)
#> Exported graph components in /tmp/Rtmpclu4Bi/file1f0b12dbdbe0_nodes.tsv and /tmp/Rtmpclu4Bi/file1f0b12dbdbe0_edges.tsv
```

Alternatively, the `igraph` objects can also be directly sent to
external software such as Cytoscape via the RESTful API implemented in
`RCy3` - this can also allow some settings of the nodes to be
conveniently set in an automated manner. We refer to the `RCy3` package
vignettes for further reading - possibly, this can be the best place to
start: [02. Cytoscape and
igraph](https://bioconductor.org/packages/3.23/RCy3/vignettes/Cytoscape-and-iGraph.html).

## Caching information with KEGGemUP

You might have noticed that when creating the KEGG graph, `KEGGemUP`
does not ask you where to store the retrieved KGML files, and neither
does it work in the current working directory. This is because
`KEGGemUP` makes extensive use of the functionality provided by
*[BiocFileCache](https://bioconductor.org/packages/3.23/BiocFileCache)*
package - and chooses this by default, instead of downloading to a local
folder. This ideally avoids the hassle of having to perform repeated
queries for files that have been already retrieved in other parallel
projects.

The functions
[`retrieve_kgml()`](https://imbeimainz.github.io/KEGGemUP/reference/retrieve_kgml.md)
and
[`get_kegg_db()`](https://imbeimainz.github.io/KEGGemUP/reference/get_kegg_db.md)
all make use of the `bfc` parameter, together with the `path`
specification, to control this aspect. The convenient wrapper
[`retrieve_all_pathways()`](https://imbeimainz.github.io/KEGGemUP/reference/retrieve_all_pathways.md)
retrieves in a single call all pathways for a given organism (to be
specified in the `org` parameter with a KEGG organism code, such as
“hsa”, “mmu”, “dme”, etc.).

``` r

local_kgml_file <- retrieve_kgml(pathway_id = "hsa04110", path = tempdir())
local_kgml_file
#> [1] "/tmp/Rtmpclu4Bi/hsa04110.xml"

graph_parsed <- create_kegg_graph(pathway_id = "local_kgml", kgml_file = local_kgml_file)
#> Warning in value[[3L]](cond): Could not retrieve pathway name for ID:
#> local_kgml
graph_parsed
#> IGRAPH de0b9e2 DN-B 134 119 -- 
#> + attr: title (g/c), type (g/c), name (v/c), type (v/c), link (v/c),
#> | reaction (v/c), reaction_link (v/c), reaction_label (v/c),
#> | graphics_name (v/c), node_name (v/c), KEGG (v/c), components (v/c),
#> | line_id (v/c), coords (v/c), graphics_type (v/c), point_index (v/n),
#> | de_value (v/n), de_text (v/c), de_source (v/c), de_name (v/c),
#> | ids_for_mapping (v/c), label (v/c), x (v/n), y (v/n), width (v/n),
#> | height (v/n), source (v/c), fgcolor (v/c), bgcolor (v/c), group
#> | (v/c), vertex.color (v/c), size (v/n), fixed (v/l), shape (v/c), type
#> | (e/c), relation_type (e/c), relation_subtype_name (e/c),
#> | relation_subtype_value (e/c), reaction_alt_name_substrate (e/c),
#> | reaction_alt_name_product (e/c), reaction_id (e/c), reaction_name
#> | (e/c), reaction_type (e/c), reaction_from_name (e/c),
#> | reaction_to_name (e/c), point_index (e/n), directed (e/l), color
#> | (e/c), width (e/n), value (e/n), label (e/c), lty (e/n), arrow.mode
#> | (e/c), dashes (e/l), title (e/c), font.face (e/c), borderWidth (e/n),
#> | borderWidthSelected (e/n), from (e/c), to (e/c)
#> + edges from de0b9e2 (vertex names):

cpd_db <- get_kegg_db(db_name = "compound", path = tempdir())
head(cpd_db)
#>   kegg_id
#> 1  C00001
#> 2  C00002
#> 3  C00003
#> 4  C00004
#> 5  C00005
#> 6  C00006
#>                                                                                                                                                  description
#> 1                                                                                                                                                 H2O; Water
#> 2                                                                                                                             ATP; Adenosine 5'-triphosphate
#> 3                                                         NAD+; NAD; Nicotinamide adenine dinucleotide; DPN; Diphosphopyridine nucleotide; Nadide; beta-NAD+
#> 4                                                                                                      NADH; DPNH; Reduced nicotinamide adenine dinucleotide
#> 5                                                                                           NADPH; TPNH; Reduced nicotinamide adenine dinucleotide phosphate
#> 6 NADP+; NADP; Nicotinamide adenine dinucleotide phosphate; beta-Nicotinamide adenine dinucleotide phosphate; TPN; Triphosphopyridine nucleotide; beta-NADP+
```

The functions
[`display_cache_KEGGemUP()`](https://imbeimainz.github.io/KEGGemUP/reference/display_cache_KEGGemUP.md)
and
[`reset_cache_KEGGemUP()`](https://imbeimainz.github.io/KEGGemUP/reference/reset_cache_KEGGemUP.md)
show and delete, respectively, all entries for KEGGemUP in the standard
BiocFileCache location.

``` r

kmu_cached <- display_cache_KEGGemUP()
head(kmu_cached)
#> $kegg
#> # A tibble: 3 × 10
#>   rid   rname create_time access_time rpath rtype fpath last_modified_time etag 
#>   <chr> <chr> <chr>       <chr>       <chr> <chr> <chr>              <dbl> <chr>
#> 1 BFC1  http… 2026-07-14… 2026-07-14… /hom… web   http…                 NA NA   
#> 2 BFC2  http… 2026-07-14… 2026-07-14… /hom… web   http…                 NA NA   
#> 3 BFC3  http… 2026-07-14… 2026-07-14… /hom… web   http…                 NA NA   
#> # ℹ 1 more variable: expires <dbl>
#> 
#> $mappings
#> # A tibble: 5 × 10
#>   rid   rname create_time access_time rpath rtype fpath last_modified_time etag 
#>   <chr> <chr> <chr>       <chr>       <chr> <chr> <chr>              <dbl> <chr>
#> 1 BFC1  http… 2026-07-14… 2026-07-14… /hom… web   http…                 NA NA   
#> 2 BFC2  http… 2026-07-14… 2026-07-14… /hom… web   http…                 NA NA   
#> 3 BFC3  http… 2026-07-14… 2026-07-14… /hom… web   http…                 NA NA   
#> 4 BFC4  http… 2026-07-14… 2026-07-14… /hom… web   http…                 NA NA   
#> 5 BFC5  http… 2026-07-14… 2026-07-14… /hom… web   http…                 NA NA   
#> # ℹ 1 more variable: expires <dbl>
```

``` r

## This will delete all previously retrieved KGML files, use this judiciously
reset_cache_KEGGemUP()
```

## Related methods & packages

KEGG pathways are used in many bioinformatics contexts and have been
around for more than two decades. Given their importance in data
integration and visualization, it is not surprising that a number of
other software packages have been developed to interface to this
resource. Our aim with `KEGGemUP` was to stress the usefulness of KEGG
pathway maps within interactive components of analysis workflows.

We are reporting them in this section for completeness, as many of them
can beautifully interoperate with `KEGGemUP`, possibly a few conversion
steps away.

- [pathview](https://bioconductor.org/packages/pathview), which is an
  R/Biocondutor package for pathway-based data integration and
  visualization - limited to static visualization onto png/pdf plots
  (Luo and Brouwer 2013).
- [KEGGgraph](https://bioconductor.org/packages/KEGGgraph), also on
  Bioconductor, offering a graph approach to KEGG pathways (based on
  `graph` and `Rgraphviz`) (Zhang and Wiemann 2009).
- [graphite](https://bioconductor.org/packages/graphite), a Bioconductor
  package to handle graph interactions from the pathway topological
  environments (Sales et al. 2012).
- [ggkegg](https://bioconductor.org/packages/ggkegg), a more recent
  Bioconductor package bringing the grammar of graphics principles (like
  for `ggplot2`) to the task of analyzing and visualizing KEGG
  information (Sato et al. 2023).
- [KEGGREST](https://bioconductor.org/packages/KEGGREST), providing
  low-level client-side REST access to the KEGG database.
- [PinPath](https://github.com/SyNUM-lab/PinPath), a package for
  visualizing omics data onto pathway diagrams (from KEGG or
  WikiPathways), and pinpoint where in the pathway the relevant changes
  occur. PinPath offers the option to save svg vectorized images of the
  plots created.
- [Cytoscape](https://cytoscape.org/), a general purpose open source
  software platform (Shannon et al. 2003), used among others for
  visualizing complex networks and integrating these with any type of
  attribute data - for which, as mentioned above, the RCy3 package
  exists to access and control directly within R.

**Disclaimer**  
`KEGGemUP` uses KEGG API, and therefore is provided for academic use by
academic users belonging to academic institutions. See
<https://www.kegg.jp/kegg/rest/> for more information.

## FAQs

**Does `KEGGemUP` play along nicely with other packages?**

Yes, since it leverages the class from the widely adopted `igraph`
package, for which implementations exist in a variety of languages.

**Why is the information obtained from parsing KGML files not 100%
coinciding with the png files provided by KEGG?**

Some inconsistencies have been identified while parsing KGML files, and
not all elements of the image files are fully included in the KGML
format. For this reason, you might notice some lines/graphical elements
missing. While this is overcome by other packages that provide their
representations on top of the original png, they are all static views,
that do not play well in an interactive HTML widget.

**I do like KEGG pathways layouts, but can I simply use other layouts?**

Yes, since we are using `igraph` to represent these objects, it is
fairly easy to switch to alternative views on the same graphs,
overriding the x and y coordinates we have retrieved from the KGML
files. We opted to keep that view as default as it is the one which
matches the “expected” aspect of the KEGG pathway map of interest.

**Is there support to create reports on KEGG enrichment analyses?**

Yes, you can simply loop over the content of a KEGG enrichment table,
retrieve their graphs, map the DE values onto them, and render them (or
some subsets of interest), all within the same analysis notebook - HTML
widgets are seamlessly working in RMarkdown/Quarto documents!

**Can I use the rendered graphs within Shiny and interact with them?**

Yes! Since the rendered graphs are created within `visNetwork`, there
are many options to do so, that fully leverage the capability of having
shiny bindings that can be used in other elements of the apps. If you
want to see some examples for this, a very good starting point is [the
visNetwork section on
Shiny](https://datastorm-open.github.io/visNetwork/shiny.html).

**Can I do customize \[this/that\] within the `visNetwork` framework?**

Possibly - a number of options are actually exposed to the end user in a
convenient manner via parameters, but of course one can achieve the
highest level of flexibility and customization interfacing directly to
the vis.js Javascript library. This, still, might not be desired, as
many of the information one wants to access and view are coming from R
and Bioconductor-centric workflows.

## Session info

This vignette was compiled on the following system:

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] KEGGemUP_0.99.1  BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] KEGGREST_1.52.2     gtable_0.3.6        xfun_0.60          
#>  [4] bslib_0.11.0        ggplot2_4.0.3       httr2_1.2.3        
#>  [7] visNetwork_2.1.4    htmlwidgets_1.6.4   vctrs_0.7.3        
#> [10] tools_4.6.1         generics_0.1.4      stats4_4.6.1       
#> [13] curl_7.1.0          tibble_3.3.1        RSQLite_3.53.3     
#> [16] blob_1.3.0          pkgconfig_2.0.3     dbplyr_2.6.0       
#> [19] RColorBrewer_1.1-3  S7_0.2.2            desc_1.4.3         
#> [22] S4Vectors_0.50.1    lifecycle_1.0.5     compiler_4.6.1     
#> [25] farver_2.1.2        textshaping_1.0.5   Biostrings_2.80.1  
#> [28] Seqinfo_1.2.0       htmltools_0.5.9     sass_0.4.10        
#> [31] yaml_2.3.12         pillar_1.11.1       pkgdown_2.2.1      
#> [34] crayon_1.5.3        jquerylib_0.1.4     cachem_1.1.0       
#> [37] tidyselect_1.2.1    digest_0.6.39       purrr_1.2.2        
#> [40] dplyr_1.2.1         bookdown_0.47       labeling_0.4.3     
#> [43] fastmap_1.2.0       grid_4.6.1          cli_3.6.6          
#> [46] magrittr_2.0.5      utf8_1.2.6          withr_3.0.3        
#> [49] filelock_1.0.3      scales_1.4.0        rappdirs_0.3.4     
#> [52] bit64_4.8.2         rmarkdown_2.31      XVector_0.52.0     
#> [55] httr_1.4.8          igraph_2.3.3        bit_4.6.0          
#> [58] otel_0.2.0          ragg_1.5.2          png_0.1-9          
#> [61] memoise_2.0.1       evaluate_1.0.5      knitr_1.51         
#> [64] IRanges_2.46.0      BiocFileCache_3.2.0 rlang_1.3.0        
#> [67] glue_1.8.1          DBI_1.3.0           xml2_1.6.0         
#> [70] BiocManager_1.30.27 BiocGenerics_0.58.1 jsonlite_2.0.0     
#> [73] R6_2.6.1            systemfonts_1.3.2   fs_2.1.0
```

## References

Antonov, Michael, Gábor Csárdi, Szabolcs Horvát, et al. 2023. “Igraph
Enables Fast and Robust Network Analysis Across Programming Languages.”
*arXiv Preprint arXiv:2311.10260*, ahead of print.
<https://doi.org/10.48550/arXiv.2311.10260>.

Kanehisa, M. 2000. “KEGG: Kyoto Encyclopedia of Genes and Genomes.”
*Nucleic Acids Research* 28 (1): 27–30.
<https://doi.org/10.1093/nar/28.1.27>.

Luo, Weijun, and Cory Brouwer. 2013. “Pathview: An r/Bioconductor
Package for Pathway-Based Data Integration and Visualization.”
*Bioinformatics* 29 (14): 1830–31.
<https://doi.org/10.1093/bioinformatics/btt285>.

Sales, Gabriele, Enrica Calura, Duccio Cavalieri, and Chiara Romualdi.
2012. “Graphite - a Bioconductor Package to Convert Pathway Topology to
Gene Network.” *BMC Bioinformatics* 13 (1).
<https://doi.org/10.1186/1471-2105-13-20>.

Sato, Noriaki, Miho Uematsu, Kosuke Fujimoto, Satoshi Uematsu, and Seiya
Imoto. 2023. “\<I\>ggkegg\</i\>: Analysis and Visualization of KEGG Data
Utilizing the Grammar of Graphics.” *Bioinformatics* 39 (10).
<https://doi.org/10.1093/bioinformatics/btad622>.

Shannon, Paul, Andrew Markiel, Owen Ozier, et al. 2003. “Cytoscape: A
Software Environment for Integrated Models of Biomolecular Interaction
Networks.” *Genome Research* 13 (11): 2498–504.
<https://doi.org/10.1101/gr.1239303>.

Zhang, Jitao David, and Stefan Wiemann. 2009. “KEGGgraph: A Graph
Approach to KEGG PATHWAY in r and Bioconductor.” *Bioinformatics* 25
(11): 1470–71. <https://doi.org/10.1093/bioinformatics/btp167>.

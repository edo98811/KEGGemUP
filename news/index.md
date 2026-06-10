# Changelog

## KEGGemUP 0.99.1

- Finalized some changes upon the Bioconductor package reviewing process
  - [`render_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/render_kegg_graph.md)
    gains a `graph_title` parameter to override the displayed title
  - fix in the vignette for correct rendering of the link to the RCy3
    vignette
  - added information on the `fnd` in the Authors section + formally
    acknowledged the contributions of the wet-lab collaborators (as
    `ctb` authors)

## KEGGemUP 0.99.0

- Ready for the submission to Bioconductor!

## KEGGemUP 0.9.1

- Better tooltip implementation
- Highlighting of the nodes now takes care of graying out a bit more of
  elements that should indeed be not too much into the focus
- Added
  [`export_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/export_kegg_graph.md)
  as a simple function to export the graphs to text files, for maximum
  compatibility with any framework (possibly even within Cytoscape)
- Providing an example dataset already precomputed within the package to
  avoid additional dependency load; as a companion to this, added
  detailed information on how to create this in the `inst/scripts`
  folder
- Finalized a full version of the vignette
- Completed the configuration of the pkgdown documentation
- Exporting also the
  [`cleanup_title_node()`](https://imbeimainz.github.io/KEGGemUP/reference/cleanup_title_node.md)
  function to remove the title node in pathway graphs

## KEGGemUP 0.9.0

- Restructuring the API to use declarative function names, that possibly
  better convey the piece of functionality
  ([`create_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/create_kegg_graph.md)
  and
  [`render_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/render_kegg_graph.md)
  as the main workhorse)
- The functions to retrieve the kgml files are now a bit more verbose if
  needed
  ([`retrieve_kgml()`](https://imbeimainz.github.io/KEGGemUP/reference/retrieve_kgml.md)
  and
  [`retrieve_all_pathways()`](https://imbeimainz.github.io/KEGGemUP/reference/retrieve_all_pathways.md))
- Subsetting and highlighting KEGG graphs is now possible with
  [`subset_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/subset_kegg_graph.md)
  and
  [`highlight_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/highlight_kegg_graph.md),
  subsequently to be passed to
  [`render_kegg_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/render_kegg_graph.md)
- The cache & download is handled in a more homogeneous manner, with
  [`display_cache_KEGGemUP()`](https://imbeimainz.github.io/KEGGemUP/reference/display_cache_KEGGemUP.md)
  and
  [`reset_cache_KEGGemUP()`](https://imbeimainz.github.io/KEGGemUP/reference/reset_cache_KEGGemUP.md)
  to check and reset the info retrieved
- Mapping continuous values from different DE-results like containers is
  handled by
  [`map_results_to_graph()`](https://imbeimainz.github.io/KEGGemUP/reference/map_results_to_graph.md)

## KEGGemUP 0.2.0

- Essential functionality implemented, from kgml files all the way down
  to rendering interactively

## KEGGemUP 0.1.0

- Getting the package ready with the full set of original features!

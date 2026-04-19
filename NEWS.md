# KEGGemUP 0.9.0

* Restructuring the API to use declarative function names, that possibly better convey the piece of functionality (`create_kegg_graph()` and `render_kegg_graph()` as the main workhorse)
* The functions to retrieve the kgml files are now a bit more verbose if needed (`retrieve_kgml()` and `retrieve_all_pathways()`)
* Subsetting and highlighting KEGG graphs is now possible with `subset_kegg_graph()` and `highlight_kegg_graph()`, subsequently to be passed to `render_kegg_graph()`
* The cache & download is handled in a more homogeneous manner, with `display_cache_KEGGemUP()` and `reset_cache_KEGGemUP()` to check and reset the info retrieved
* Mapping continuous values from different DE-results like containers is handled by `map_results_to_graph()`

# KEGGemUP 0.2.0

* Essential functionality implemented, from kgml files all the way down to rendering interactively

# KEGGemUP 0.1.0

* Getting the package ready with the full set of original features!

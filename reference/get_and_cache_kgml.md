# Download and cache KEGG KGML files.

Download and cache KEGG KGML files.

## Usage

``` r
get_and_cache_kgml(pathway_id, bfc = NULL, file_name = NULL)
```

## Arguments

- pathway_id:

  KEGG pathway ID (e.g., "hsa04110").

- bfc:

  BiocFileCache object for caching KEGG KGML files.

- file_name:

  Optional file name to save the KGML file directly.

## Value

Path to the cached KGML file.

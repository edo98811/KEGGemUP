# Download a KGML file

Download and cache KEGG KGML files.

## Usage

``` r
retrieve_kgml(pathway_id, bfc = NULL, path = NULL, verbose = FALSE)
```

## Arguments

- pathway_id:

  Character, KEGG pathway ID (e.g., 'hsa04110').

- bfc:

  BiocFileCache object for caching KEGG KGML files. Defaults to NULL

- path:

  Optional path to save the KGML file if not using cache. Defaults to
  NULL

- verbose:

  Logical, if TRUE, print additional messages.

## Value

Path to the cached KGML file.

## Details

This function can save an individual KGML to the cache or to a path. If
specified, the BiocFileCache cache is prioritized. If a `path` if
specified and it is a directory, it will save the file there. If the
path specified a file path, this will be the location to store the file.
If it is left as NULL, it will save the KMGL file in the current working
directory

## Examples

``` r
data_dir <- tempdir()
kgml_path <- retrieve_kgml("hsa04110", path = data_dir, verbose = TRUE)
#> Downloading KGML from: https://rest.kegg.jp/get/hsa04110/kgml
#> Downloaded & saved in: /tmp/RtmpHtNn53/hsa04110.xml

cache_dir <- tempdir()
kgml_path_cached <- retrieve_kgml("hsa04110",
  bfc = BiocFileCache::BiocFileCache(cache_dir),
  verbose = TRUE
)
#> Downloading KGML from: https://rest.kegg.jp/get/hsa04110/kgml
#> adding rname 'https://rest.kegg.jp/get/hsa04110/kgml'
#> 
#> Cached: hsa04110
```

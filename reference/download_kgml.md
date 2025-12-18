# Download and cache KEGG KGML files.

Download and cache KEGG KGML files.

## Usage

``` r
download_kgml(pathway_id, bfc = NULL, directory = NULL)
```

## Arguments

- pathway_id:

  KEGG pathway ID (e.g., 'hsa04110').

- bfc:

  BiocFileCache object for caching KEGG KGML files.

- directory:

  Optional directory to save the KGML file if not using cache.

## Value

Path to the cached KGML file.

## Examples

``` r
data_dir <- tempdir()
kgml_path <- download_kgml("hsa04110", directory = data_dir)
#> Downloaded & saved in: /tmp/RtmpwO4Si1/hsa04110.xml
```

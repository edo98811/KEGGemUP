# Download and cache KEGG KGML files.

Download and cache KEGG KGML files.

## Usage

``` r
download_kgml(pathway_id, bfc = NULL, directory = NULL, verbose = FALSE)
```

## Arguments

- pathway_id:

  KEGG pathway ID (e.g., 'hsa04110').

- bfc:

  BiocFileCache object for caching KEGG KGML files.

- directory:

  Optional directory to save the KGML file if not using cache.

- verbose:

  Logical, if TRUE, print additional messages.

## Value

Path to the cached KGML file.

## Examples

``` r
data_dir <- tempdir()
kgml_path <- download_kgml("hsa04110", directory = data_dir, verbose = TRUE)
#> Downloading KGML from: https://rest.kegg.jp/get/hsa04110/kgml
#> Downloaded & saved in: /tmp/Rtmpf7h8iD/hsa04110.xml
```

# Return all cached KEGG and mapping files from BiocFileCache

Return all cached KEGG and mapping files from BiocFileCache

## Usage

``` r
return_all_cached()
```

## Value

A list containing data frames of cached KEGG and mapping files

## Details

This function retrieves information about all cached KEGG pathway files
and mapping files stored using BiocFileCache.

## Examples

``` r
cache_info <- return_all_cached()
print(cache_info$kegg)      # View cached KEGG pathway files
#> # A tibble: 1 × 10
#>   rid   rname create_time access_time rpath rtype fpath last_modified_time etag 
#>   <chr> <chr> <chr>       <chr>       <chr> <chr> <chr>              <dbl> <chr>
#> 1 BFC1  http… 2026-02-23… 2026-02-23… /hom… web   http…                 NA NA   
#> # ℹ 1 more variable: expires <dbl>
print(cache_info$mappings)  # View cached mapping files
#> # A tibble: 5 × 10
#>   rid   rname create_time access_time rpath rtype fpath last_modified_time etag 
#>   <chr> <chr> <chr>       <chr>       <chr> <chr> <chr>              <dbl> <chr>
#> 1 BFC1  http… 2026-02-23… 2026-02-23… /hom… web   http…                 NA NA   
#> 2 BFC2  http… 2026-02-23… 2026-02-23… /hom… web   http…                 NA NA   
#> 3 BFC3  http… 2026-02-23… 2026-02-23… /hom… web   http…                 NA NA   
#> 4 BFC4  http… 2026-02-23… 2026-02-23… /hom… web   http…                 NA NA   
#> 5 BFC5  http… 2026-02-23… 2026-02-23… /hom… web   http…                 NA NA   
#> # ℹ 1 more variable: expires <dbl>
```

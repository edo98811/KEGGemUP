# Reset all KEGGemUP caches

Reset the KEGG and mapping caches by deleting all cached files

## Usage

``` r
reset_cache_KEGGemUP()
```

## Value

(invisible) NULL

## Details

This function deletes all cached KEGG pathway files and mapping files
stored using BiocFileCache. It prompts the user for confirmation before
proceeding with the deletion.

## Examples

``` r
if (FALSE) { # \dontrun{
  reset_cache_KEGGemUP()
} # }
```

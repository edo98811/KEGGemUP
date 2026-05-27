# Download all pathways

Download all KEGG pathways for a given organism

## Usage

``` r
retrieve_all_pathways(org, wait = 0.5, verbose = FALSE)
```

## Arguments

- org:

  KEGG organism code (e.g., 'hsa' for human).

- wait:

  Numeric value, needs to be strictly positive. Indicates the amount in
  seconds to wait in between requests, being polite and respectful of
  the limit rates imposed by KEGG. Defaults to 0.5, which is safely a
  bit above the rate of 3max/sec.

- verbose:

  Logical, if TRUE, print additional messages.

## Value

The BiocFileCache object is returned invisibly

## Details

This function can take a while to run completely, but is an efficient
way to retrieve all the kgml files encoding for the KEGG pathways. Some
failures could be triggered by a rate limitation imposed by the KEGG
website. Since this caches the files locally, it is safe to rerun to
complete the operation without extra unneeded requests to the API

## Examples

``` r
# Download all pathways for human
if (FALSE) { # \dontrun{
  retrieve_all_pathways("hsa", verbose = TRUE)
} # }
```

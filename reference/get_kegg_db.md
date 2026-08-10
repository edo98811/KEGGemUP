# Get a KEGG database file

Get KEGG database, with caching.

## Usage

``` r
get_kegg_db(db_name = "compound", path = NULL, bfc = NULL, verbose = FALSE)
```

## Arguments

- db_name:

  Name of the KEGG database to retrieve (default: "compound").

- path:

  Optional path to save the KEGG database file if not using cache.
  Defaults to NULL.

- bfc:

  BiocFileCache object for caching KEGG database files. Defaults to NULL

- verbose:

  Logical, if TRUE, print additional messages.

## Value

A data frame with KEGG IDs and names.

## Details

The valid KEGG database names are: kegg \| pathway \| brite \| module \|
ko \| genes \| \| vg \| vp \| ag \| genome \| ligand \| compound \|
glycan \| reaction \| rclass \| enzyme \| network \| variant \| disease
\| drug \| dgroup

## Examples

``` r
# Saving in path
data_dir <- tempdir()
kegg_compounds <- get_kegg_db(
  db_name = "compound",
  path = data_dir, verbose = TRUE
)
#> Retrieving KEGG database: compound
#> Downloaded & saved KEGG database in: /tmp/RtmpkwIYXN/kegg_compound.tsv

# Just returning without saving
kegg_compounds_onthefly <- get_kegg_db(
  db_name = "compound",
  verbose = TRUE
)
#> Retrieving KEGG database: compound
#> No 'bfc' or 'path' provided.
head(kegg_compounds_onthefly)
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

# saving to cache (in a temp dir)
kegg_compounds_cached <- get_kegg_db(
  db_name = "compound",
  bfc = BiocFileCache::BiocFileCache(tempdir()),
  verbose = TRUE
)
#> Retrieving KEGG database: compound
#> adding rname 'https://rest.kegg.jp/list/compound'
#> 
#> Cached KEGG database: compound
```

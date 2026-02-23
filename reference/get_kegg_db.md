# Get KEGG db with caching.

Get KEGG db with caching.

## Usage

``` r
get_kegg_db(
  db_name = "compound",
  directory = NULL,
  bfc = NULL,
  verbose = FALSE
)
```

## Arguments

- db_name:

  Name of the KEGG database to retrieve (default: "compound").

- directory:

  Optional directory to save the KEGG database file if not using cache.

- bfc:

  BiocFileCache object for caching KEGG database files.

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
# Saving in directory
data_dir <- tempdir()
kegg_compounds <- get_kegg_db("compound", directory = data_dir, verbose = TRUE)
#> Retrieving KEGG database: compound
#> Downloaded & saved KEGG database in: /tmp/Rtmpf7h8iD/kegg_compound.tsv
# Just returning without saving
kegg_genes <- get_kegg_db("compound", verbose = TRUE)
#> Retrieving KEGG database: compound
#> No 'bfc' or 'directory' provided. Not saving KEGG database only downloading and returning.
```

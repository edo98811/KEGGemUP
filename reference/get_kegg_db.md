# Get KEGG db with caching.

Get KEGG db with caching.

## Usage

``` r
get_kegg_db(db_name = "compound", directory = NULL, bfc = NULL)
```

## Arguments

- db_name:

  Name of the KEGG database to retrieve (default: "compound").

- directory:

  Optional directory to save the KEGG database file if not using cache.

- bfc:

  BiocFileCache object for caching KEGG database files.

## Value

A data frame with KEGG IDs and names.

## Details

The valid KEGG database names are: kegg \| pathway \| brite \| module \|
ko \| genes \| \| vg \| vp \| ag \| genome \| ligand \| compound \|
glycan \| reaction \| rclass \| enzyme \| network \| variant \| disease
\| drug \| dgroup

If neither 'bfc' nor 'directory' is provided, the KEGG database will be
downloaded but not saved. It will be returned as a data frame.

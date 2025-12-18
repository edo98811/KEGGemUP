# Get KEGG db with caching.

Get KEGG db with caching.

## Usage

``` r
get_kegg_db(bfc, db_name = "compound")
```

## Arguments

- bfc:

  A BiocFileCache object for caching.

- db_name:

  KEGG database name (e.g., 'compound', 'glycan').

## Value

A data frame with KEGG IDs and names.

## Details

The valid KEGG database names are: kegg \| pathway \| brite \| module \|
ko \| genes \| \| vg \| vp \| ag \| genome \| ligand \| compound \|
glycan \| reaction \| rclass \| enzyme \| network \| variant \| disease
\| drug \| dgroup

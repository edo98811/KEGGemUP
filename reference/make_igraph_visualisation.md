# Plot igraph with improved layout and outlier handling

Plot igraph with improved layout and outlier handling

## Usage

``` r
make_igraph_visualisation(
  g,
  eliminate_distance_outliers = TRUE,
  size_multiplier = 0.3,
  text_dist = 1.1,
  text_cex = 0.5
)
```

## Arguments

- g:

  An igraph object to plot.

- eliminate_distance_outliers:

  Logical, if TRUE, replaces outlier node positions with mean positions
  to improve layout visualization (default: TRUE).

- size_multiplier:

  Numeric multiplier to adjust node sizes in the plot (default: 0.3).

- text_dist:

  Numeric distance for node labels from the nodes (default: 1.1).

- text_cex:

  Numeric scaling factor for node label text size (default

## Value

A plot of the igraph object with improved layout.

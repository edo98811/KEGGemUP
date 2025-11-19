# Combine multiple differential expression results into a single data frame

Combine multiple differential expression results into a single data
frame

## Usage

``` r
combine_results_in_dataframe(results_list)
```

## Arguments

- results_list:

  A named list where each element is a differential expression result
  containing a data frame (de_table), value column name (value_column),
  and feature column name (feature_column)

## Value

A combined data frame with columns: KEGG, value, source

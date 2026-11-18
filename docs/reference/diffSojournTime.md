# Difference of sojourn times

This function computes the difference between two sojourn times.

## Usage

``` r
diffSojournTime(x1, x2)
```

## Arguments

- x1:

  a sojournTime object, called with `draws = TRUE` and
  `returndraws = TRUE`, computed for a specified seed.

- x2:

  a sojournTime object, called with `draws = TRUE` and
  `returndraws = TRUE`, computed with the same seed as `x1`.

## Value

returns the median, the 2.5% and 97.5% quantiles, the mean, the standard
deviation of the difference between the two sojourn times and the number
of removed draws (eventually due to computational issues)

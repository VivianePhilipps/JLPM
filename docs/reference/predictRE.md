# Random effects prediction

The random effects of a `jointLPM` model are predicted as the mode of
the conditional distribution of the random effects given the subject's
longitudinal and survival observations.

## Usage

``` r
predictRE(x, data, control = NULL, variance = FALSE)
```

## Arguments

- x:

  a `jointLPM` object

- data:

  a data.frame object containing the data from the predictions are to be
  computed

- control:

  optional list of arguments to be passed to the `marqLevAlg`
  optimization function.

- variance:

  logical indicating whether the variance of the estimated predictions
  should be returned. Default ot FALSE.

## Value

a matrix containing the predicted random effects (in columns) for each
subject (in rows). If `variance = TRUE`, returns a list containing the
previous matrix and a vector containing the upper part of each variance
matrix.

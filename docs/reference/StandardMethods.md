# Standard methods for estimated models

coef and vcov methods for estimated jointLPM models.

## Usage

``` r
# S3 method for class 'jointLPM'
coef(object, ...)
```

## Arguments

- object:

  an object of class `jointLPM`

- ...:

  other arguments. There are ignored in these functions.

## Value

For `coef`, the vector of the estimated parameters.

For `vcov`, the variance-covariance matrix of the estimated parameters.

## Author

Viviane Philipps

# Summary of a joint latent process model

This function provides a summary of model estimated with the `jointLPM`
function.

## Usage

``` r
# S3 method for class 'jointLPM'
summary(object, ...)
```

## Arguments

- object:

  an object inheriting from class `jointLPM` for a joint latent process
  model.

- ...:

  further arguments to be passed to or from other methods. They are
  ignored in this function.

## Value

The function is mainly used for its side effects. It returns invisibily
a list of two matrices containing the estimates, their standard errors,
Wald statistics and associated p-values for the survival submodel (first
element of the list) and for the mixed model's fixed effects (second
element).

## Author

Viviane Philipps, Tiphaine Saulnier and Cecile Proust-Lima

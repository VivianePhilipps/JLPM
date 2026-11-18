# createX0

Create the matrix containing all covariates appearing in the formula

## Usage

``` r
createX0(form, data)
```

## Arguments

- form:

  a list of formula

- data:

  a data frame containing all variables mentionned in form

## Value

a list containing

- X0:

  the pooled (with cbind) model matrix. The number of lines is the same
  as in 'data'

- idform:

  a list of the same length as 'form' specifying which term of X0 appear
  in each formula

## Details

Caution 1 : 'data' should not contain any missing values. There is no
check for that within the function. Caution 2 : if interaction terms are
present, the names are sorted. For example x:age becomes age:x in X0.
Caution 3 : intercept is not systematically in the first column. It will
only if form\[\[1\]\] contains an intercept.

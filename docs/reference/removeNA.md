# removeNA

Remove missing values

## Usage

``` r
removeNA(form, data)
```

## Arguments

- form:

  a list of formula

- data:

  a data frame containing all variables mentionned in form

## Value

a list containing

- newdata:

  a data frame created from data, containing only the variables used in
  the formula, and without any missing values

- Y:

  the names of the response variables

- YlinesNA:

  the lines removed from each response variables

- nmes:

  the number of measures for each response variable

- X:

  the names of the covariates

- XlinesNA:

  the lines removed from data because of missing values in the
  covariates

## Details

Caution : if several response variables are included in 'form', then the
resulting data frame will contain one column (named outcome) grouping
all these response variables. The irst lines will refer to the first
response variables, a second block to the second, etc. All covariates
will be repeated for each response variable.

## Examples

``` r
 d <- data.frame(y1=c(3,6,4,2), y2=c(6.5,NA,9.5,4.5), x=c(0,0,1,0))
 removeNA(form=list(y1+y2~x), data=d)
#> $newdata
#>    outcome x
#> 1      3.0 0
#> 2      6.0 0
#> 3      4.0 1
#> 4      2.0 0
#> 11     6.5 0
#> 31     9.5 1
#> 41     4.5 0
#> 
#> $Y
#> [1] "y1" "y2"
#> 
#> $YlinesNA
#> $YlinesNA[[1]]
#> integer(0)
#> 
#> $YlinesNA[[2]]
#> [1] 2
#> 
#> 
#> $nmes
#> [1] 4 3
#> 
#> $X
#> [1] "x"
#> 
#> $XlinesNA
#> NULL
#> 
#> $terms
#> $terms[[1]]
#> ~x
#> attr(,"variables")
#> list(x)
#> attr(,"factors")
#>   x
#> x 1
#> attr(,"term.labels")
#> [1] "x"
#> attr(,"order")
#> [1] 1
#> attr(,"intercept")
#> [1] 1
#> attr(,"response")
#> [1] 0
#> attr(,".Environment")
#> <environment: 0x574d453d6f80>
#> attr(,"predvars")
#> list(x)
#> attr(,"dataClasses")
#>         x 
#> "numeric" 
#> 
#> 
## 1      3.0 0   | first block referring to y1
## 2      6.0 0   |
## 3      4.0 1   |
## 4      2.0 0   |
## 11     6.5 0     | second block referring to y2
## 31     9.5 1     |
## 41     4.5 0     |
```

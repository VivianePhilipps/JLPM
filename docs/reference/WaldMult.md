# WaldMult Performs multivariate Wald tests.

Given a vector of estimated parameters (coef = (theta_1, ... theta_n)')
and their associated variance matrix (vcov = Var(coef)), the WaldMult
function performs either the test of the null hypothesis :

## Usage

``` r
WaldMult(
  Mod,
  coef,
  vcov,
  pos = NULL,
  contrasts = NULL,
  A = NULL,
  name = NULL,
  value = 0
)
```

## Arguments

- Mod:

  an object for which the `coef` and `vcov` methods are diefined.

- coef:

  if `Mod` is missing, the vector of estimated parameters

- vcov:

  if `Mod` is missing, the matrix of variance of coef

- pos:

  vector containing the indices of the parameters involved in the test.
  Required if A is not specified.

- contrasts:

  optional vector containing the the c1, ..., cm value in case 2 (see
  details). Should be of same lengtht as pos.

- A:

  matrix A of case 3 (see details). A should have as many columns as
  parameters in coef.

- name:

  characters string giving the name of the test printed in the output
  (the row names of the output). By default, the name's test is the null
  hypothesis.

- value:

  the value against which to test. By default, value=0.

## Value

If contrasts is NULL, the function returns a matrix with 1 row and 2
columns containing the value of the Wald test's statistic and the
associated p-value.

If contrasts is not NULL, the function returns a matrix with 1 row and 4
columns containing the value of the coefficient (dot product of pos and
contrasts), his standard deviation, the value of the Wald test's
statistic and the associated p-value.

If A is specified, the function returns a matrix with m rows and 2
columns containing the value of the Wald test's statistic and the
associated p-value.

## Details

1\. H0 : theta_m = (theta_j1, ..., theta_jm)' = 0 ie, tests if a subset
of m parameters are simultaneously equal to zero. The Wald statistic is
w = theta_m' \* Var(theta_m)^(-1) \* theta_m and the associated p-value
is obtained by p = P(chi2\_(dof=m) \> w)

2\. H0 : theta = c1\*theta_j1 + ... + cm\*theta_jm = 0 ie, tests if a
combination of m parameters is equal to zero. The Wald statistic is w =
theta^2 / Var(theta) = (c1\*theta_j1 + ... + cm\*theta_jm)^2 /
Var(c1\*theta_j1 + ... + cm\*theta_jm) and the associated p-value is
obtained by p = P(chi2\_(dof=1) \> w)

3\. H0 : Theta = A\*coef = 0 with A a matrix with m lines and n columns.
ie, tests if multiple combinations of parameters are simultaneously
equal to zero. The Wald statistic is w = Theta \* Var(Theta)^(-1) \*
Theta' = (A\*coef)' \* (A\*vcov\*A')^(-1) \* (A\*coef) and the
associated p-value is obtained by p = P(chi2\_(dof=m) \> w)

The function can be applied to any object having a `coef`

## Examples

``` r
require(lcmm)
m <- hlme(Y ~ Time + X2, random = ~ 1, subject = "ID", data = data_hlme)
summary(m)
#> Heterogenous linear mixed model 
#>      fitted by maximum likelihood method 
#>  
#> hlme(fixed = Y ~ Time + X2, random = ~1, subject = "ID", data = data_hlme)
#>  
#> Statistical Model: 
#>      Dataset: data_hlme 
#>      Number of subjects: 100 
#>      Number of observations: 326 
#>      Number of latent classes: 1 
#>      Number of parameters: 5  
#>  
#> Iteration process: 
#>      Convergence criteria satisfied 
#>      Number of iterations:  11 
#>      Convergence criteria: parameters= 3.4e-07 
#>                          : likelihood= 6.3e-07 
#>                          : second derivatives= 4.2e-14 
#>  
#> Goodness-of-fit statistics: 
#>      maximum log-likelihood: -932.2  
#>      AIC: 1874.4  
#>      BIC: 1887.42  
#>  
#>  
#> Maximum Likelihood Estimates: 
#>  
#> Fixed effects in the longitudinal model:
#> 
#>               coef      Se   Wald p-value
#> intercept 25.87421 0.87824 29.462 0.00000
#> Time      -0.46375 0.07725 -6.003 0.00000
#> X2         2.80425 1.59293  1.760 0.07833
#> 
#> 
#> Variance-covariance matrix of the random-effects:
#>           intercept
#> intercept  49.09354
#> 
#>                              coef      Se
#> Residual standard error:  2.57999 0.12135
#> 

mm <- hlme(Y ~ -1 + Time + factor(X2), random = ~ 1, subject = "ID", data = data_hlme)
summary(mm)
#> Heterogenous linear mixed model 
#>      fitted by maximum likelihood method 
#>  
#> hlme(fixed = Y ~ -1 + Time + factor(X2), random = ~1, subject = "ID", 
#>     data = data_hlme)
#>  
#> Statistical Model: 
#>      Dataset: data_hlme 
#>      Number of subjects: 100 
#>      Number of observations: 326 
#>      Number of latent classes: 1 
#>      Number of parameters: 5  
#>  
#> Iteration process: 
#>      Convergence criteria satisfied 
#>      Number of iterations:  14 
#>      Convergence criteria: parameters= 1.4e-07 
#>                          : likelihood= 2.1e-07 
#>                          : second derivatives= 1.7e-14 
#>  
#> Goodness-of-fit statistics: 
#>      maximum log-likelihood: -932.2  
#>      AIC: 1874.4  
#>      BIC: 1887.42  
#>  
#>  
#> Maximum Likelihood Estimates: 
#>  
#> Fixed effects in the longitudinal model:
#> 
#>                 coef      Se   Wald p-value
#> Time        -0.46375 0.07725 -6.003 0.00000
#> factor(X2)0 25.87421 0.87824 29.462 0.00000
#> factor(X2)1 28.67846 1.37732 20.822 0.00000
#> 
#> 
#> Variance-covariance matrix of the random-effects:
#>           intercept
#> intercept  49.09354
#> 
#>                              coef      Se
#> Residual standard error:  2.57999 0.12135
#> 

## Retrieve from model mm the difference between X2 levels as in model m :
WaldMult(Mod = mm, pos = c(2, 3), contrasts = c(-1, 1))
#>                                            coef      Se Wald Test p_value
#> factor(X2)0 * -1 + factor(X2)1 * 1 =  0 2.80425 1.59293   3.09914 0.07833
## or
WaldMult(coef = coef(mm), vcov = vcov(mm), A = matrix(c(0, -1, 1, 0, 0), 1, 5))
#>               Wald Test p_value
#> A*coef = (0)'   3.09914 0.07833
```

# Computation of sojourn time under specific condition

This function computes posterior computations from a joint shared
random-effect model with ordinal longitudinal outcomes, aka a joint item
response theory model. Specifically, the function computes the expected
time spent before reaching a given level of impairment specified by one
or multiple items for a covariate profile.

## Usage

``` r
sojournTime(
  x,
  maxState,
  condState = NULL,
  newdata,
  var.time,
  startTime = 0,
  nMC = 1000,
  upperTime = 150,
  subdivisions = 100L,
  rel.tol = .Machine$double.eps^0.25,
  draws = FALSE,
  ndraws = 2000,
  returndraws = FALSE,
  cl = NULL,
  seed = NULL
)
```

## Arguments

- x:

  an object of class `jointLPM` representing of joint shared random
  effects model with binary or ordinal longitudinal outcome(s)

- maxState:

  a list specifying the items and the corresponding levels defining the
  maximum state for the computation of the sojourn time. For instance
  `maxState=list(Y = 3)` will compute the expected sojourn time
  corresponding to an impairment of Y lower or equal to 3.

- condState:

  an optional list specifying the initial state at start time (argument
  `startTime`) from which to compute the residual sojourn time. The
  state can be defined either as a single value or as an interval. For
  example `condState=list(Y = 2)` means that the state at start time is
  equal to 2, whereas `condState=list(Y = c(1, 2))` means that the state
  at start time is between 1 and 2 (both included).

- newdata:

  a one line data frame specifying the covariate profile for which the
  sojourn time is computed.

- var.time:

  a character string specifying the name of the time variable in the
  longitudinal submodel. Note that this time covariate should not be
  included in newdata.

- startTime:

  a numeric value specifying the time from which the residual sojourn
  time is computed (the lower bound of the integral over time). Default
  to 0 for the expected sojourn time over lifetime.

- nMC:

  an integer giving the number of Monte Carlo simulations used to
  compute the integral over the random effects.

- upperTime:

  a numeric specifying the upper bound of the integral over time.
  Default to 150 (150 years as an approximation of infinity).

- subdivisions:

  passed to the `integrate` function.

- rel.tol:

  passed to the `integrate` function.

- draws:

  logical indicating if 95% confidence interval should be computed.
  Default to FALSE.

- ndraws:

  integer giving the number of draws to be used to compute the 95%
  confidence interval. Default to 2000.

- returndraws:

  logical indicating if the `ndraws` results should be returned. Default
  to FALSE.

- cl:

  either a cluster created with `makeCluster` or an integer specifying
  the number of cores that should be used for computation. Only used
  with draws = TRUE.

- seed:

  integer only used with `draws = TRUE` to set the random seed.

## Value

if `draws = FALSE`, returns a single value. If `draws = TRUE` and
`returndraws = FALSE`, returns the median, the 2.5% and 97.5% quantiles,
the mean, the standard deviation and the number of removed draws
(eventually due to computational issues). If `draws = TRUE` and
`returndraws = TRUE`, returns the `ndraws` values.

## Details

1\. Lifetime expected sojourn time: the expected time before reaching
level k at item(s) Y is the integral over time t, from 0 to infinity, of
P(Y(t) \<= k, T \> t), i.e., int_0^infty P(Y(t) \<= k, T \> t) dt.

2\. Residual expected sojourn time from a time s: conditionally on being
below m at item(s) Y at time s and being alive at time s, the sojourn
time below level k at item(s) Y is computed as: int_s^infty P(Y(t) \<=
k, T \> t \| Y(s) \<= m, T \> s) dt = int_s^infty P(Y(t) \<= k, T \> t,
Y(s) \<= m) dt / P(Y(s) \<= m, T \> s)

## Author

Viviane Philipps and Cecile Proust-Lima

## Examples

``` r
library(lcmm)
paq <- paquid[which(paquid$age_init < paquid$agedem), ]
paq$age65 <- (paq$age - 65) / 10
paq$ageinit65 <- (paq$age_init - 65) / 10
paq$agedem65 <- (paq$agedem - 65) / 10

#### Estimation of the joint model with one ordinal longitudinal item
if (FALSE) { # \dontrun{
M2 <- jointLPM(fixed = HIER ~ age65 * male,
                random = ~ age65,
                subject = "ID", 
                link = "thresholds",
                survival = Surv(ageinit65, agedem65, dem) ~ male,
                sharedtype = 'RE',
                var.time = "age65",
                data = paq, 
                methInteg = "QMC", 
                nMC = 1000,
                B = c(0.6, 2.399, -0.409, 2.076, 6.338, 0.994, -0.223, -0.005,
                     -0.299, 0.174, 0.523, 1.044, 1.064, 0.506))

#### Computation of the expected lifetime sojourn time with HIER impairment
#### up to 2 (HIER = 2)
#### (=int_0^150 P(HIER(t) <= 2, T > t) dt)
sojournTime(M2, list(HIER = 2), newdata = data.frame(male = 0), 
var.time = "age65")

#### Computation of the expected residual time with maximum HIER impairment
#### of 2 (HIER = 2), given the impairment was 1 at time 0.5
#### (=int_0.5^150 P(HIER(t) <= 2, T > t | HIER(0.5) <= 1, T > 0.5) dt)
sojournTime(M2, list(HIER = 2), condState = list(HIER = 1), startTime = 0.5,
newdata = data.frame(male = 0), var.time = "age65")
} # }
```

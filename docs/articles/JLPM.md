# Joint Latent Process Models

The R package JLPM implements in the jointLPM function the estimation of
a joint shared random effects model.

 

The longitudinal data ($y_{1},y_{2},\ldots,y_{K}$) can be continuous or
ordinal and are modelled using a mixed model.

 

For the continuous case, a transformation $H_{k}$ is estimated for
outcome $y_{k}$ ($k \in 1,\ldots,K$) :

$$H_{k}\left( y_{k}\left( t_{ijk} \right),\eta_{k} \right) = X_{i}\left( t_{ijk} \right)\beta + Z_{i}\left( t_{ijk} \right)u_{i} + \varepsilon_{ijk}$$

For an ordinal outcome, mixed models are combined to the Item Response
Theory:

$$\left. {\mathbb{P}}\left( y_{k}\left( t_{ijk} \right) = m \right)\Leftrightarrow{\mathbb{P}}\left( \eta_{k,m - 1} < X_{i}\left( t_{ijk} \right)\beta + Z_{i}\left( t_{ijk} \right)u_{i} + \varepsilon_{ijk} < \eta_{k,m} \right) \right.$$

where $X_{i}\left( t_{ijk} \right)$ and $Z_{i}\left( t_{ijk} \right)$
are vectors of covariates measured at time $t_{ijk}$ for subject $i$,
$\beta$ is the vector of fixed effects, $u_{i} \sim \mathcal{N}(0,B)$
are the random effects,
$\varepsilon_{ijk} \sim \mathcal{N}\left( 0,\sigma_{k}^{2} \right)$ the
measurement errors, $\eta_{k}$ are the parameters of the link function
$H_{k}$ or the thresholds associated to the outcome $y_{k}$.

Note that even with multiple longitudinal outcomes, a univariate mixed
model is estimated.

 

The time-to-event data are modelled in a proportional hazard model where
different associations between the longitudinal outcome and the event
can be included.

 

With an association through the random effects of the longitudinal
model:

$$\alpha_{i}(t) = \alpha_{0}(t,\omega)\exp\left( {\widetilde{X}}_{i}\gamma + \delta u_{i} \right)$$

With an association through the current level of the longitudinal model:

$$\alpha_{i}(t) = \alpha_{0}(t,\omega)\exp\left( {\widetilde{X}}_{i}\gamma + \delta\left( X_{i}(t)\beta + Z_{i}(t)u_{i} \right) \right)$$

With an association through the current slope of the longitudinal model:

$$\alpha_{i}(t) = \alpha_{0}(t,\omega)\exp\left( {\widetilde{X}}_{i}\gamma + \delta\left( \frac{d}{dt}X_{i}(t)\beta + \frac{d}{dt}Z_{i}(t)u_{i} \right) \right)$$

with $\alpha_{0}(t,\omega)$ the baseline hazard function at time $t$,
parameterized with $\omega$, ${\widetilde{X}}_{i}$ a vector of
time-fixed covariates, $\gamma$ the fixed effects, and $\delta$ the
association parameter.

 

The parameters $\beta,B,\sigma_{k},\eta_{k},\omega,\gamma,\delta$ are
estimated using a Marquardt-Levenberg algorithm.

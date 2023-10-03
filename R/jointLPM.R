#' Estimation of latent process joint models for multivariate longitudinal 
#' outcomes and time-to-event data.
#' 
#' This function fits extended joint models with shared random effects.
#' The longitudinal submodel handles multiple continuous longitudinal outcomes
#' (Gaussian or  non-Gaussian, curvilinear) as well as 
#' ordinal longitudinal outcomes in a mixed effects framework. The model assumes that all the outcomes 
#' measure the same underlying latent process defined as their common
#' factor, and each outcome is related to this latent common factor by
#' outcome-specific measurement models whose nature depends on the type of the 
#' associated outcome (linear model for Gaussian outcome, curvilinear model for 
#' non-Gaussian outcome, cumulative probit model for ordinal outcome).
#' At the latent process level, the model estimates a standard linear mixed model.
#' The survival submodel handles right-censored (possibly left-truncated) time-to-events with competing risks.
#' The association between the longitudinal and the survival data is captured by including 
#' the random effect from the mixed model or the predicted current level 
#' of the underlying process as a linear predictor in the proportional hazard survival model.
#' Parameters of the measurement models, of the latent process mixed model and of the 
#' survival model are estimated simultaneously using a maximum likelihood method,
#' through a Marquardt-Levenberg algorithm.
#' 
#' 
#' 
#' A. THE MEASUREMENT MODELS
#' 
#' \code{jointLPM} function estimates one measurement model per outcome to link
#' each outcome Y_k(t) with the underlying latent common factor L(t) they measure.
#' To fix the latent process dimension, we chose to constrain at the latent process 
#' level the intercept of the mixed model at 0 and the 
#' standard error of the first random effect at 1. The nature of each measurment
#' model adapts to the type of the outcome it models. 
#' 
#' 1. For continuous Gaussian outcomes, linear models are used and required 2 parameters for 
#' the transformation (Y(t) - b1)/b2
#' 
#' 2. For continuous non-Gaussian outcomes, curvilinear models use 
#' parametrized link function to link outcomes to the latent process. 
#' With the "beta" link function, 4 parameters are required for the
#' following transformation: [ h(Y(t)',b1,b2) - b3]/b4 where h is the Beta CDF
#' with canonical parameters c1 and c2 that can be derived from b1 and b2 as
#' c1=exp(b1)/[exp(b2)*(1+exp(b1))] and c2=1/[exp(b2)*(1+exp(b1))], and Y(t)'
#' is the rescaled outcome i.e. Y(t)'= [ Y(t) - min(Y(t)) + epsY ] / [
#' max(Y(t)) - min(Y(t)) +2*epsY ].
#' With the "splines" link function, n+2 parameters are required for the
#' following transformation b_1 + b_2*I_1(Y(t)) + ... + b_{n+2}*I_{n+1}(Y(t)),
#' where I_1,...,I_{n+1} is the basis of quadratic I-splines. To constraint the
#' parameters to be positive, except for b_1, the program estimates b_k^* (for
#' k=2,...,n+2) so that b_k=(b_k^*)^2.
#' 
#' 3. For discrete ordinal outcomes, cumulative probit models are used. For a
#' (n+1)-level outcome, the model consist of determining n thresholds t_k in the 
#' latent process scale which correspond to the outcome level changes. Then,
#' Y(t) = n' <=> t_n' < L(t) + e <= t_(n'+1) with e the standard error of the outcome.
#' To ensure that t_1 < t_2 < ... < t_n, the program estimates t'_1, t'_2, ..., t'_n
#' such that t_1=t'_1, t_2=t_1+(t'_2)^2, t_3=t_2+(t'_3)^2, ...
#' 
#' 
#' B. THE SURVIVAL MODEL
#' 
#' a. BASELINE RISK FUNCTIONS
#' 
#' For the baseline risk functions, the following parameterizations were considered. 
#' 
#' 1. With the "Weibull" function: 2 parameters are necessary w_1 and w_2 so that 
#' the baseline risk function a_0(t) = w_1^2*w_2^2*(w_1^2*t)^(w_2^2-1) if logscale=FALSE  
#' and \cr
#' a_0(t) = exp(w_1)*exp(w_2)(t)^(exp(w_2)-1) if logscale=TRUE.
#' 
#' 2. with the "piecewise" step function and nz nodes (y_1,...y_nz), nz-1 parameters 
#' are necesssary p_1,...p_nz-1 so that the baseline risk function a_0(t) = p_j^2 
#' for y_j < t =< y_j+1 if logscale=FALSE and a_0(t) = exp(p_j) for y_j < t =< y_j+1 if logscale=TRUE.
#' 
#' 3. with the "splines" function and nz nodes (y_1,...y_nz), nz+2 parameters 
#' are necessary s_1,...s_nz+2 so that the baseline risk function a_0(t) = sum_j s_j^2 M_j(t) 
#' if logscale=FALSE and a_0(t) = sum_j exp(s_j) M_j(t) if logscale=TRUE where M_j is the basis of cubic M-splines.
#' Two parametrizations of the baseline risk function are proposed (logscale=TRUE or FALSE) 
#' because in some cases, especially when the instantaneous risks are very close to 0, 
#' some convergence problems may appear with one parameterization or the other. 
#' As a consequence, we recommend to try the alternative parameterization (changing logscale option) 
#' when a model does not converge (maximum number of iterations reached) and
#' where convergence criteria based on the parameters and likelihood are small.
#' 
#' 
#' b. ASSOCIATION BETWEEN LONGITUDINAL AND SURVIVAL DATA
#' 
#' The association between the longitudinal and the survival data is captured by including 
#' a function of the elements from the latent process mixed model as a predictor in the survival model.
#'  We implement two association structures,
#' that should be specified through \code{sharedtype} argument.
#' 
#' 1. the random effect from the latent process linear mixed model (\code{sharedtype='RE'}) :
#' the q random effects modeling the individual deviation in the longitudinal model are also included
#' in the survival model, so that a q-vector of parameters measures the association
#' between the risk of event and the longitudinal outcome(s).
#' 
#' 2. the predicted current level of the underlying process (\code{sharedtype='CL'}) :
#' the predicted latent process defined by the mixed model appears as
#' time-dependent covariate in the survival model.
#' The association between the longitudinal process and the risk of event
#' is then quantified by a unique parameter.
#' 
#'  
#' C. THE VECTOR OF PARAMETERS B
#' 
#' The parameters in the vector of initial values \code{B} or in the vector of
#' maximum likelihood estimates \code{best} are included in the following
#' order:
#' (1) parameters for the baseline risk function: 2 parameters for each Weibull, 
#' nz-1 for each piecewise constant risk and nz+2 for each splines risk. In the 
#' presence of competing events, the number of parameters should be adapted to 
#' the number of causes of event;
#' (2) for all covariates in survival, one parameter 
#' is required. Covariates parameters should be included in the same order as in survival.
#' In the presence of cause-specific effects, the number of parameters should be
#' multiplied by the number of causes;
#' (3) parameter(s) of association between the longitudinal 
#' and the survival process: for \code{sharedtype='RE'}, one parameter per random effect
#' and per cause of event is 
#' required; for \code{sharedtype='CL'}, one parameter per cause of event is required;
#' (4) for all covariates  in fixed, one parameter is required. Parameters should
#' be included in the same  order as in fixed;
#' (5)for all covariates included with \code{contrast()} in \code{fixed}, one
#' supplementary parameter per outcome is required excepted for the last
#' outcome for which the parameter is not estimated but deduced from the others;
#' (6) the variance of each random-effect specified in random 
#' (excepted the intercept which is constrained to 1) 
#' if idiag=TRUE and the inferior triangular variance-covariance matrix of all 
#' the random-effects if idiag=FALSE;
#' (7) if \code{cor} is specified, the standard error of the Brownian motion or 
#' the standard error and the correlation parameter of the autoregressive process;
#' (8) parameters of each measurement model: 2 for "linear", 4 for "beta", 
#' n+2 for "splines" with n nodes, n for "thresholds" for a (n+1)-level outcome; 
#' (9) if \code{randomY=TRUE}, the standard 
#' error of the outcome-specific random intercept (one per outcome); 
#' (10) the outcome-specific standard errors (one per outcome)
#' 
#' 
#' 
#' C. CAUTIONS REGARDING THE USE OF THE PROGRAM
#' 
#' Some caution should be made when using the program. Convergence criteria are
#' very strict as they are based on the derivatives of the log-likelihood in
#' addition to the parameter and log-likelihood stability. In some cases, the
#' program may not converge and reach the maximum number of iterations fixed at
#' 100 by default. In this case, the user should check that parameter estimates at the
#' last iteration are not on the boundaries of the parameter space.
#' 
#' If the parameters are on the boundaries of the parameter space, the
#' identifiability of the model is critical. This may happen especially with
#' splines parameters that may be too close to 0 (lower boundary). When
#' identifiability of some parameters is suspected, the program can be run
#' again from the former estimates by fixing the suspected parameters to their
#' value with option posfix. This usually solves the problem. An alternative is
#' to remove the parameters of the Beta of Splines link function from the
#' inverse of the Hessian with option partialH. If not, the program should be 
#' run again with other initial values, with a higher maximum number of iterations 
#' or less strict convergence tolerances.
#' 
#' To reduce the computation time, this program can be carried out in parallel mode,
#' ie. using multiple cores which number can be specified with argument \code{nproc}.
#' 
#' 
#' 
#' @param fixed a two-sided linear formula object for specifying the
#' fixed-effects in the linear mixed model at the latent process level. The
#' response outcomes are separated by \code{+} on the left of \code{~} and the
#' covariates are separated by \code{+} on the right of the \code{~}. For
#' identifiability purposes, the intercept specified by default should not be
#' removed by a \code{-1}. Variables on which a contrast above the different
#' outcomes should also be estimated are included with \code{contrast()}.
#' @param random a one-sided formula for the random-effects in the
#' latent process mixed model. Covariates with a random-effect are separated
#' by \code{+}. An intercept should always be included for identifiability.
#' @param subject name of the covariate representing the grouping structure.
#' @param idiag optional logical for the variance-covariance structure of the
#' random-effects. If \code{FALSE}, a non structured matrix of
#' variance-covariance is considered (by default). If \code{TRUE} a diagonal
#' matrix of variance-covariance is considered.
#' @param randomY optional logical for including an outcome-specific random
#' intercept. If \code{FALSE} no outcome-specific random intercept is added
#' (default). If \code{TRUE} independent outcome-specific random intercepts
#' with parameterized variance are included.
#' @param link optional vector of families of parameterized link functions defining
#' the measurement models (one by outcome). Option "linear" (by default) specifies a linear
#' link function. Other possibilities include "beta" for estimating a link
#' function from the family of Beta cumulative distribution functions, "Splines" 
#' for approximating the link function by I-splines and "thresholds" for ordinal
#' outcomes modelled by cumulative probit models. For splines case, the number of 
#' nodes and the nodes location should be also specified. The number of nodes is 
#' first entered followed by \code{-}, then the location is specified with "equi", 
#' "quant" or "manual" for respectively equidistant
#' nodes, nodes at quantiles of the marker distribution or interior nodes
#' entered manually in argument \code{intnodes}. It is followed by \code{-splines}.
#' For example, "7-equi-splines" means
#' I-splines with 7 equidistant nodes, "6-quant-splines" means I-splines with 6
#' nodes located at the quantiles of the marker distribution and
#' "9-manual-splines" means I-splines with 9 nodes, the vector of 7 interior
#' nodes being entered in the argument \code{intnodes}.
#' @param intnodes optional vector of interior nodes. This argument is only
#' required for a I-splines link function with nodes entered manually.
#' @param epsY optional positive real used to rescale the marker in
#' (0,1) when the beta link function is used. By default, epsY=0.5.
#' @param cor optional indicator for inclusion of an autocorrelated Gaussian
#' process in the linear mixed model at the latent process level. Option
#' \code{BM(time)} indicates a brownian motion with parameterized variance while
#' option \code{AR(time)} specifies an autoregressive process in continuous time
#' with parameterized variance and correlation intensity. In both cases, \code{time}
#' is the variable representing the measurement times. By default, no autocorrelated
#' Gaussian process is added.
#' @param survival two-sided formula object. The left side of the formula corresponds 
#' to a Surv() object for right-censored (\code{Surv(EntryTime,Time,Indicator)})
#' and possibly left-truncated (\code{Surv(EntryTime,Time,Indicator)}).
#' Multiple causes of event can be considered 
#' in the Indicator (0 for censored, k for event of cause k). The right side of the 
#' formula specifies the covariates to include in the survival model.
#' Cause-specific covariate effect are specified with \code{cause()}.For example,  
#' Surv(Time,Indicator) ~ X1 + cause(X2) indicates a common effect of X1 and a cause-specific effect of X2. 
#' @param hazard optional family of hazard function assumed for the survival model. 
#' By default, "Weibull" specifies a Weibull baseline risk function. Other possibilities 
#' are "piecewise" for a piecewise constant risk function or "splines" for a cubic M-splines 
#' baseline risk function. For these two latter families, the number of nodes and the 
#' location of the nodes should be specified as well, separated by -. The number of 
#' nodes is entered first followed by -, then the location is specified with "equi", 
#' "quant" or "manual" for respectively equidistant nodes, nodes at quantiles of the 
#' times of event distribution or interior nodes entered manually in argument hazardnodes. 
#' It is followed by - and finally "piecewise" or "splines" indicates the family of 
#' baseline risk function considered. Examples include "5-equi-splines" for M-splines 
#' with 5 equidistant nodes, "6-quant-piecewise" for piecewise constant risk over 5 
#' intervals and nodes defined at the quantiles of the times of events distribution 
#' and "9-manual-splines" for M-splines risk function with 9 nodes, the vector of 7 
#' interior nodes being entered in the argument hazardnodes. In the presence of competing 
#' events, a vector of hazards should be provided such as hazard=c("Weibull","5-quant-splines") 
#' with 2 causes of event, the first one modelled by a Weibull baseline cause-specific 
#' risk function and the second one by splines.
#' @param hazardnodes optional vector containing interior nodes if splines or piecewise 
#' is specified for the baseline hazard function in hazard.
#' @param TimeDepVar optional vector containing an intermediate time
#' corresponding to a change in the risk of event. This time-dependent
#' covariate can only take the form of a time variable with the assumption that
#' there is no effect on the risk before this time and a constant effect on the
#' risk of event after this time (example: initiation of a treatment to account
#' for).
#' @param logscale optional boolean indicating whether an exponential
#' (logscale=TRUE) or a square (logscale=FALSE -by default) transformation is
#' used to ensure positivity of parameters in the baseline risk functions. See
#' details section
#' @param startWeibull optional numeric with Weibull hazard functions only.
#' Indicates the shift in the Weibull distribution.
#' @param hazardrange optional vector indicating the range of the survival times 
#' (that is the minimum and maximum). By default, the range is defined according 
#' to the minimum and maximum observed values of the survival times. The option 
#' should be used only for piecewise constant and Splines hazard functions.
#' @param data data frame containing all variables named in \code{fixed},
#'  \code{random}, \code{cor}, \code{survival} and \code{subject}.
#' @param B optional specification for the initial values for the parameters.
#' Initial values should be entered in the order detailed in \code{details} section.
#' @param convB optional threshold for the convergence criterion based on the
#' parameter stability. By default, convB=0.0001.
#' @param convL optional threshold for the convergence criterion based on the
#' log-likelihood stability. By default, convL=0.0001.
#' @param convG optional threshold for the convergence criterion based on the
#' derivatives. By default, convG=0.0001.
#' @param maxiter optional maximum number of iterations for the Marquardt
#' iterative algorithm. By default, maxiter=100.
#' @param nsim number of points used to plot the estimated link functions. By
#' default, nsim=100.
#' @param range optional vector indicating the range of the outcomes (that is
#' the minimum and maximum). By default, the range is defined according to the
#' minimum and maximum observed values of the outcome. The option should be
#' used only for Beta and Splines transformations.
#' @param subset optional vector giving the subset of observations in
#' \code{data} to use. By default, all lines.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param posfix Optional vector giving the indices in vector B of the
#' parameters that should not be estimated. Default to NULL, all parameters are
#' estimated.
#' @param partialH optional vector giving the indices in vector B of parameters that
#' can be dropped from the Hessian matrix to define convergence criteria.
#' @param verbose logical indicating if information about computation should be
#' reported. Default to TRUE.
#' @param returndata logical indicating if data used for computation should be
#' returned. Default to FALSE, data are not returned.
#' @param methInteg character indicating the type of integration to compute the
#' log-likelihood. 'MCO' for ordinary Monte Carlo, 'MCA' for antithetic Monte Carlo,
#' 'QMC' for quasi Monte Carlo. Default to "QMC".
#' @param nMC integer, number of Monte Carlo simulations. Default to 1000.
#' @param sharedtype indicator of shared random function type. \code{'RE'} indicates
#' an association through the random effects included in the linear mixed model.
#' \code{'CL'} defines a association through the predicted current level of the latent process.
#' @param var.time name of the variable representing the measurement times.
#' @param nproc number of cores for parallel computation.
#' @param clustertype the type of cluster that should internally be created.
#' See \code{parallel::makeCluster} for possible values.
#' 
#' @return An object of class "jointLPM" is returned containing some internal information
#' used in related functions. Users may investigate the following elements : 
#' \item{ns}{number of grouping units in the dataset} 
#' \item{loglik}{log-likelihood of the model} 
#' \item{best}{vector of parameter estimates in the same order as specified in 
#' \code{B} and detailed in section \code{details}}
#' \item{V}{vector containing the upper triangle matrix of variance-covariance
#' estimates of \code{best} with exception for variance-covariance parameters
#' of the random-effects for which \code{V} contains the variance-covariance
#' estimates of the Cholesky transformed parameters displayed in \code{cholesky}} 
#' \item{gconv}{vector of convergence criteria: 1. on the parameters, 2. on the 
#' likelihood, 3. on the derivatives} 
#' \item{conv}{status of convergence: =1 if the convergence criteria were satisfied, 
#' =2 if the maximum number of iterations was reached, =4 or 5 if a problem occured 
#' during optimisation} 
#' \item{call}{the matched call} 
#' \item{niter}{number of Marquardt iterations} 
#' \item{nevent}{number of occured event}
#' \item{pred}{table of individual predictions and residuals in the underlying
#' latent process scale; it includes marginal predictions (pred_m), marginal
#' residuals (resid_m), subject-specific predictions (pred_ss) and
#' subject-specific residuals (resid_ss) and finally the transformed
#' observations in the latent process scale (obs).}
#' \item{predRE}{table containing individual predictions of the random-effects : 
#' a column per random-effect, a line per subject.}
#' \item{predRE_Y}{table containing individual predictions of the outcome-specific
#' random intercept}
#' \item{predSurv}{table containing the predicted baseline risk function and
#' the predicted cumulative baseline risk function }
#' \item{cholesky}{vector containing the estimates of the Cholesky transformed
#' parameters of the variance-covariance matrix of the random-effects}
#' \item{estimlink}{table containing the simulated values of each outcome and
#' the corresponding estimated link function} 
#' \item{epsY}{definite positive reals used to rescale the markers in (0,1) when 
#' the beta link function is used. By default, epsY=0.5.} 
#' \item{AIC}{the Akaike's information criterion}
#' \item{BIC}{the Bayesian information criterion}
#' \item{CPUtime}{the runtime in seconds}
#' \item{data}{the original data set (if returndata is TRUE)}
#' @author Viviane Philipps, Tiphaine Saulnier and Cecile Proust-Lima
#' 
#' @references
#' Saulnier, Philipps, Meissner, Rascol, Pavy-Le-Traon, Foubert-Samier, Proust-Lima (2021).
#' Joint models for the longitudinal analysis of measurement scales in the presence 
#' of informative dropout   arXiv:2110.02612
#' 
#' Philipps, Hejblum, Prague, Commenges, Proust-Lima (2021).
#' Robust and efficient optimization using a Marquardt-Levenberg algorithm with 
#' R package marqLevAlg, The R Journal 13:2.  
#' 
#' @examples
#' #### Examples with paquid data from R-package lcmm
#' library(lcmm)
#' paq <- paquid[which(paquid$age_init<paquid$agedem),]
#' paq$age65 <- (paq$age-65)/10
#'
#' #### Example with one Gaussian marker :
#' ## We model the cognitive test IST according to age, sexe and eduction level. We assume
#' ## a Weibull distribution for the time to dementia and link the longitudinal and survival
#' ## data using the random effects.
#' ## We provide here the call to the jointLPM function without optimization (maxiter=0). The
#' ## results should therefore not be interpreted.
#' M0 <- jointLPM(fixed = IST~age65*(male+CEP),
#'                 random=~age65,
#'                 idiag=FALSE,
#'                 subject="ID",
#'                 link="linear",
#'                 survival=Surv(age_init,agedem,dem)~male,
#'                 sharedtype='RE',
#'                 hazard="Weibull",
#'                 data=paq,
#'                 var.time="age65",
#'                 maxiter=0)
#' M0$best ## these are the initial values of each of the 15 parameters
#' 
#' \donttest{ 
#' ## Estimation with one Gaussian marker
#' ## We remove the maxiter=0 option to estimate the model. We specify initial values
#' ## to reduce the runtime, but this can take several minutes.
#' binit1 <- c(0.1039, 5.306, -0.1887, -1.0355, -4.3817, -1.0543, -0.1161, 0.8588,
#' 0.0538, -0.1722, -0.2224, 0.3296, 30.7768, 4.6169, 0.7396)
#' M1 <- jointLPM(fixed = IST~age65*(male+CEP),
#'                 random=~age65,
#'                 idiag=FALSE,
#'                 subject="ID",
#'                 link="linear",
#'                 survival=Surv(age_init,agedem,dem)~male,
#'                 sharedtype='RE',
#'                 hazard="Weibull",
#'                 data=paq,
#'                 var.time="age65",
#'                 B=binit1)
#' ## Optimized the parameters to be interpreted :
#' summary(M1)
#'
#' 
#' #### Estimation with one ordinal marker :
#' ## We consider here the 4-level hierarchical scale of dependence HIER and use "thresholds"
#' ## to model it as an ordinal outcome. We assume an association between the current level
#' ## of dependency and the risk of dementia through the option sharedtype="CL".
#' ## We use a parallel optimization on 2 cores to reduce computation time.
#' binit2 <- c(0.0821, 2.4492, 0.1223, 1.7864, 0.0799, -0.2864, 0.0055, -0.0327, 0.0017,
#' 0.3313, 0.9763, 0.9918, -0.4402)
#' M2 <- jointLPM(fixed = HIER~I(age-65)*male,
#'                 random = ~I(age-65),
#'                 subject = "ID", 
#'                 link = "thresholds",
#'                 survival = Surv(age_init,agedem,dem)~male,
#'                 sharedtype = 'CL',
#'                 var.time = "age",
#'                 data = paq, 
#'                 methInteg = "QMC", 
#'                 nMC = 1000,
#'                 B=binit2,
#'                 nproc=2)
#' summary(M2)
#' 
#' }
#' 
#' @export
#' 
jointLPM <- function(fixed,random,subject,idiag=FALSE,cor=NULL,link="linear",intnodes=NULL,epsY=0.5,randomY=FALSE, var.time,
                survival=NULL,hazard="Weibull",hazardrange=NULL,hazardnodes=NULL,TimeDepVar=NULL,logscale=FALSE,startWeibull=0, sharedtype='RE',
                methInteg="QMC",nMC=1000,data,subset=NULL,na.action=1,
                B,posfix=NULL,maxiter=100,convB=0.0001,convL=0.0001,convG=0.0001,partialH=NULL,
                nsim=100,range=NULL,verbose=TRUE,returndata=FALSE,
                nproc=1, clustertype=NULL)
{
    ptm <- proc.time()
    
    cl <- match.call()
    
    nom.subject <- as.character(subject)
    
    if(missing(random)) stop("At least one random effect is required")
    if(random==~-1) stop("At least one random effect is required")
    if(missing(fixed)) stop("The argument fixed must be specified in any model")
    if(!inherits(fixed,"formula")) stop("The argument fixed must be a formula")
    if(!inherits(random,"formula")) stop("The argument random must be a formula")
    if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
    if(nrow(data)==0) stop("Data should not be empty")
    if(missing(subject)){ stop("The argument subject must be specified")}
    if(!is.numeric(data[,subject])) stop("The argument subject must be numeric")
    if(all(link %in% c("linear","beta","thresholds")) & !is.null(intnodes)) stop("Intnodes should only be specified with splines links")
    if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")
    if(!(sharedtype%in%c('RE','CL')))stop("The value of argument sharedtype must be 'RE' (for shared random effects) or 'CL' (for shared latent process current level)")
    if(missing(var.time) | length(var.time)!=1)stop("The argument var.time is missing or is not of length 1")
    if(!(var.time %in% colnames(data))) stop("Unable to find variable 'var.time' in 'data'")
    if(sharedtype == 'CL' & missing(cor)==FALSE) print("WARNING : wi not computed on current level prediction 'cause considered not shared")
    if(sharedtype == 'CL' & missing(TimeDepVar)==FALSE) stop("model with sharedtype='CL' not yet programmed with time dependent effect on survival")
    
    #    if(length(posfix) & missing(B)) stop("A set of initial parameters must be specified if some parameters are not estimated")
    
    
    ## garder data tel quel pour le renvoyer
    if(returndata==TRUE)
    {
        datareturn <- data
    }
    else
    {
        datareturn <- NULL
    }
    
    ## test de l'argument cor
    ncor <- 0
    cor.type <- cl$cor[1]
    cor.time <- cl$cor[2]
    cor <- paste(cor.type,cor.time,sep="-")
    if (!isTRUE(all.equal(cor,character(0))))
    {
        if (substr(cor,1,2)=="AR") { ncor <- 2 }
        else if (substr(cor,1,2)=="BM") { ncor <- 1  }
        else { stop("The argument cor must be of type AR or BM") }
        
        if(!(strsplit(cor,"-")[[1]][2] %in% colnames(data))) stop("Unable to find time variable from argument 'cor' in 'data'")
        else { cor.var.time <- strsplit(cor,"-")[[1]][2] }
    }
    
    ##pour acces aux attributs des formules
    afixed <- terms(fixed, specials=c("factor","contrast"))
    if(attr(afixed,"intercept")==0) stop("An intercept should appear in fixed for identifiability purposes")
    
    arandom <- terms(random, specials=c("factor"))
    ##fixed sans contrast
    fixed2 <- gsub("contrast","",fixed)
    fixed2 <- formula(paste(fixed2[2],fixed2[3],sep="~"))   
    afixed2 <- terms(fixed2)
    
    ##verifier si toutes les varialbes sont dans data
    variables <- c(attr(afixed,"variables"),attr(arandom,"variables"))
    variables <- unlist(lapply(variables,all.vars))  
    if(!all(variables %in% colnames(data))) stop(paste("Data should contain the variables",paste(unique(variables),collapse=" ")))
    
    
    ##contrast
    contr <- ~-1
    if(!is.null(attr(afixed,"specials")$contrast))
    {
        vcontr <- attr(afixed,"term.labels")[setdiff(attr(afixed,"specials")$contrast-1,untangle.specials(afixed,"contrast",2)$terms)]
        vcontr <- gsub("contrast","",vcontr)
        contr <- as.formula(paste("~-1+",paste(vcontr,collapse="+")))
    }
    acontr <- terms(contr)
    
    
    ##liste des outcomes
    nomsY <- as.character(attr(afixed,"variables")[2])
    nomsY <- strsplit(nomsY,split=" + ",fixed=TRUE)
    nomsY <- as.vector(nomsY[[1]])
    ny <- length(nomsY)
    
    ##pas de contrast ni randomY si un seul Y
    if(ny<2 & length(attr(afixed,"specials")$contrast)) stop("No contrast can be included with less than two outcomes")
    if(ny<2 & randomY==TRUE) stop("With less than 2 outcomes randomY should be FALSE")
    
    ##liste des variables utilisees  (sans les interactions et sans les Y)
    ttesLesVar <- colnames(get_all_vars(afixed,data=data[1,]))
    ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(arandom,data=data[1,])))
    if (ncor>0) ttesLesVar <- unique(c(ttesLesVar,cor.var.time))
    else ttesLesVar <- unique(ttesLesVar)
    ttesLesVar <- setdiff(ttesLesVar, nomsY)
    
    ## argument subset
    form1 <- paste(c(nom.subject,nomsY,ttesLesVar),collapse="+")
    if(!isTRUE(all.equal(as.character(cl$subset),character(0))))
    {
        cc <- cl
        cc <- cc[c(1,which(names(cl)=="subset"))]
        cc[[1]] <- as.name("model.frame")
        cc$formula <- formula(paste("~",form1))
        cc$data <- data
        cc$na.action <- na.pass
        data <- eval(cc)
    }
    
    attributes(data)$terms <- NULL
    
    ## si subject est un factor
    if(is.factor(data[,nom.subject]))
    {
        data[,nom.subject] <- as.numeric(data[,nom.subject])
    }
    
    ## partie survie
    if(is.null(survival))
    {
        nbevt <- 0
        idtrunc <- 0
        nprisq <- 0
        nrisqtot <- 0
        nvarxevt <- 0
        nvarxevt2 <- 0
        typrisq <- 0
        noms.surv <- NULL
        nom.timedepvar <- NULL
        form.commun <- ~-1
        form.cause <- ~-1
        survival <- ~-1
        nz <- 0
        zi <- 0
        minT <- 0
        maxT <- 0
    }
    else
    {
        ## objet Surv
        surv <- cl$survival[[2]]
        
        if(length(surv)==3) #censure droite sans troncature gauche
        {
            idtrunc <- 0 
            
            Tevent <- getElement(object=data,name=as.character(surv[2]))
            Event <- getElement(object=data,name=as.character(surv[3]))  
            Tentry <- rep(0,length(Tevent)) #si pas de troncature, Tentry=0
            
            noms.surv <-  c(as.character(surv[2]),as.character(surv[3]))
            
            surv <- do.call("Surv",list(time=Tevent,event=Event,type="mstate"))  
        }
        
        if(length(surv)==4) #censure droite et troncature
        {
            idtrunc <- 1
            
            Tentry <- getElement(object=data,name=as.character(surv[2]))
            Tevent <- getElement(object=data,name=as.character(surv[3]))
            Event <- getElement(object=data,name=as.character(surv[4]))  
            
            noms.surv <-  c(as.character(surv[2]),as.character(surv[3]),as.character(surv[4]))   
            
            surv <- do.call("Surv",list(time=Tentry,time2=Tevent,event=Event,type="mstate"))   
        }  
        
        ## nombre d'evenement concurrents
        nbevt <- length(attr(surv,"states"))
        if(nbevt<1) stop("No observed event in the data")
        
        
        ## pour la formule pour survivial, creer 3 formules : 
        ## une pour les covariables en mixture, une pour les covariables avec effet specifique a la cause, et une pour les effets communs.  
        form.surv <- cl$survival[3]
        
        noms.form.surv <- all.vars(attr(terms(formula(paste("~",form.surv))),"variables"))
        if(length(noms.form.surv)==0)
        {
            form.cause <- ~-1
            form.commun <- ~-1
            asurv <- terms(~-1)
        }
        else
        {
            ##creer la formula pour cause
            form1 <- gsub("mixture","",form.surv)       #TS: pas de classes latentes
            form1 <- formula(paste("~",form1))
            asurv1 <- terms(form1,specials="cause")  
            ind.cause <- attr(asurv1,"specials")$cause
            if(length(ind.cause))
            {
                form.cause <- paste(labels(asurv1)[ind.cause],collapse="+")
                form.cause <- gsub("cause","",form.cause)
                form.cause <- formula(paste("~",form.cause))
            }
            else
            {
                form.cause <- ~-1 
            }
            
            
            ## creer la formule pour ni cause ni mixture
            asurv <- terms(formula(paste("~",form.surv)),specials=c("cause"))
            ind.commun <- setdiff(1:length(labels(asurv)),unlist(attr(asurv,"specials")))
            if(length(ind.commun))
            {
                form.commun <- paste(labels(asurv)[ind.commun],collapse="+")
                form.commun <- gsub("cause","",form.commun)   # si X1:cause(X2)
                form.commun <- formula(paste("~",form.commun))  
            }
            else
            {
                form.commun <- ~-1 
            }
        }
    }
    
    nom.timedepvar <- NULL
    if(!missing(TimeDepVar))
    {
        if(!is.null(TimeDepVar))
        {  
            nom.timedepvar <- as.character(TimeDepVar)
            if(!isTRUE(nom.timedepvar %in% colnames(data))) stop(paste("Data should contain variable",nom.timedepvar))
        }
    }
    
    
    ##verifier si toutes les variables sont dans data
    varSurv <- unique(all.vars(terms(survival)))
    if(!is.null(nom.timedepvar)){if(!(nom.timedepvar %in% all.vars(terms(survival)))) stop("Variable in 'TimeDepVar' should also appear as a covariate in the 'survival' argument")}  
    if(!all(varSurv %in% colnames(data))) stop(paste("Data should contain the variables",paste(varSurv,collapse=" ")))
    ttesLesVar <- unique(c(ttesLesVar,varSurv))
    
    
    ##subset de data avec les variables utilisees
    newdata <- data[,c(nom.subject,nomsY,noms.surv,ttesLesVar),drop=FALSE]
    
    
    ## remplacer les NA de TimeDepVar par Tevent
    Tint <- Tevent
    nvdepsurv <- 0  
    if(!is.null(nom.timedepvar))
    {
        Tint <- newdata[,nom.timedepvar]
        Tint[(is.na(Tint))] <- Tevent[(is.na(Tint))]
        Tint[Tint>Tevent] <- Tevent[Tint>Tevent]
        Tint[Tint<Tentry] <- Tentry[Tint<Tentry]
        nvdepsurv <- 1
        if (length(Tint[Tint<Tevent])==0)
        {
            stop("TimeDepVar is always greater than Time of Event. \n")  
            nvdepsurv <- 0
        }
        if (length(Tint[Tint>Tentry])==0)
        {
            Tint <- Tevent
            stop("TimeDepVar is always lower than Time of Entry (0 by default). \n")
            nvdepsurv  <- 0
        }
        
        newdata[,nom.timedepvar] <- Tint 
    }
    
    dataSurv <- NULL
    if((nbevt>0))
    {
        dataSurv <- data.frame(getElement(object=data,name=nom.subject),Tentry,Tevent,Event,Tint)
    }
    
    
    ##un data frame par outcome et creation Y0
    dataY <- paste("data",nomsY,sep=".")
    Y0 <- NULL
    IND <- NULL
    outcome <- NULL
    data0 <- NULL
    nayk <- vector("list",ny)
    for (k in 1:ny)
    {
        dtemp <- newdata[,c(nom.subject,nomsY[k],noms.surv,ttesLesVar)]
        ##enlever les NA
        linesNA <- apply(dtemp,2,function(v) which(is.na(v)))
        linesNA <- unique(unlist(linesNA))
        if(length(linesNA)) nayk[[k]] <- linesNA
        if(na.action==1 & length(linesNA)>0) dtemp <- dtemp[-linesNA,]
        if(na.action==2 & length(linesNA)>0) stop("Data contains missing values")
        assign(dataY[k],dtemp)
        Y0 <- c(Y0, dtemp[,nomsY[k]])
        IND <- c(IND, dtemp[,nom.subject])
        outcome <- c(outcome,rep(nomsY[k],nrow(dtemp)))
        data0 <- rbind(data0, dtemp[,setdiff(colnames(dtemp),nomsY[k]),drop=FALSE])   #dataset sans NA avec les covariables utilisees; obs ordonnees par outcome
    }
    
    ##creation de X0 (ttes les var + interactions)
    Xfixed <- model.matrix(fixed2[-2], data=data0)
    Xrandom <- model.matrix(random, data=data0)
    Xcontr <- model.matrix(contr,data=data0)
    Xsurv <- model.matrix(form.commun,data=data0)
    Xsurvcause <- model.matrix(form.cause,data=data0)
    
    
    z.fixed <- strsplit(colnames(Xfixed),split=":",fixed=TRUE)
    z.fixed <- lapply(z.fixed,sort)
    
    z.random <- strsplit(colnames(Xrandom),split=":",fixed=TRUE)
    z.random <- lapply(z.random,sort)
    
    if(contr != ~-1)
    {
        z.contr <- strsplit(colnames(Xcontr),split=":",fixed=TRUE)
        z.contr <- lapply(z.contr,sort)
    }
    else
    {
        z.contr <- list()
    }
    
    if(form.commun != ~-1)
    {
        z.surv <- strsplit(colnames(Xsurv),split=":",fixed=TRUE)
        z.surv <- lapply(z.surv,sort)
    }
    else
    {
        z.surv <- list() 
    }
    
    if(form.cause != ~-1)
    {
        z.survcause <- strsplit(colnames(Xsurvcause),split=":",fixed=TRUE)
        z.survcause <- lapply(z.survcause,sort)
    }
    else
    {
        z.survcause <- list() 
    }
    
    
    if(!all(z.contr %in% z.fixed))  stop("The covariates in contrast should also appear in fixed")
    
    X0 <- cbind(Xfixed, Xrandom, Xsurv, Xsurvcause)
    
    nom.unique <- unique(colnames(X0))
    X0 <- X0[,nom.unique,drop=FALSE]
    
    form.cor <- ~-1
    if (ncor>0)
    {
        if(!(cor.var.time %in% colnames(X0)))
        {
            X0 <- cbind(X0, data0[,cor.var.time])
            colnames(X0) <- c(nom.unique, cor.var.time)
            nom.unique <- c(nom.unique,cor.var.time)
            form.cor <- formula(paste("~-1+",cor.var.time))
        }
    }
    
    X0 <- as.matrix(X0)
    ##X0 fini
    
    
    
    ##test de link
    if (length(link)!=1 & length(link)!=ny) stop("One link per outcome should be specified")
    if(any(link %in% c("splines","Splines")))
    {
        link[which(link %in% c("splines","Splines"))] <- "5-quant-splines"
    }
    if(length(link)==1 & ny>1)
    {
        link <- rep(link, ny)
    }
    
    idlink <- rep(2,ny)
    idlink[which(link=="linear")] <- 0
    idlink[which(link=="beta")] <- 1
    idlink[which(link=="thresholds")] <- 3
    
    spl <- strsplit(link[which(idlink==2)],"-")
    if(any(sapply(spl,length)!=3)) stop("Invalid argument 'link'")
    
    nySPL <- length(spl)
    nybeta <- sum(idlink==1)
    nyORD <- sum(idlink==3)
    ##remplir range si pas specifie
    if(!is.null(range) & length(range)!=2*(nySPL+nybeta)) stop("Length of vector range is not correct.")
    if((length(range)==2*(nySPL+nybeta)) & (nySPL+nybeta>0))
    {
        ind12 <- which(idlink==1 | idlink==2)
        for (i in 1:(nySPL+nybeta))
        {
            rg <- range(get(dataY[ind12[i]])[,nomsY[ind12[i]]])
            if(rg[1]<range[2*(i-1)+1] | rg[2]>range[2*(i-1)+2]) stop("The range specified do not cover the entire range of the data")
        }
    }
    if((is.null(range) & (nybeta+nySPL)>0) | length(range)!=2*(nySPL+nybeta))
    {
        range <- NULL
        for(k in which(idlink!=0))
        {
            min1 <- min(get(dataY[k])[,nomsY[k]])
            min2 <- round(min1,3)
            if(min1<min2) min2 <- min2-0.001
            
            max1 <- max(get(dataY[k])[,nomsY[k]])
            max2 <- round(max1,3)
            if(max1>max2) max2 <- max2+0.001
            
            range <- c(range, min2, max2)
        }
    }
    
    
    ## epsY
    if (any(idlink==1))
    {
        if (any(epsY<=0))
        {
            stop("Argument 'epsY' should be positive.")
        }
        
        if(length(epsY)==1) epsY <- rep(epsY,nybeta)
        
        if(length(epsY)!=nybeta) stop(paste("Argument 'epsY' should be of length",nybeta))
        if(nybeta!=ny)
        {
            epsY2 <- rep(0,ny)
            epsY2[which(idlink==1)] <- epsY
            epsY <- epsY2
        }
        
    }
    
    
    nbzitr <- rep(2,ny) #nbzitr = nb de noeuds si splines, 2 sinon
    nbnodes <- NULL  #que pour les splines
    spltype <- NULL
    if(nySPL>0)
    {
        for (i in 1:nySPL)
        {
            nbnodes <- c(nbnodes, spl[[i]][1])
            spltype <- c(spltype, spl[[i]][2])
            if(spl[[i]][3] != "splines") stop("Invalid argument link")
        }
    }
    nbnodes <- as.numeric(nbnodes)
    nbzitr[which(idlink==2)] <- nbnodes
    
    ##test splines
    if(!(all(spltype %in% c("equi","quant","manual")))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")
    
    ##tester longueur de intnodes
    if(!is.null(intnodes))
    {
        if(length(intnodes) != sum(nbnodes[which(spltype=="manual")]-2)) stop(paste("Vector intnodes should be of length",sum(nbnodes[which(spltype=="manual")]-2)))
    }
    
    ##intnodes2 : contient tous les noeuds interieurs (pas seulement ceux de manual)
    intnodes2 <- rep(NA,sum(nbnodes-2))
    nb <- 0
    nbspl <- 0
    for (k in 1:ny)
    {
        if (idlink[k]!=2) next
        else
        {
            nbspl <- nbspl+1
            
            if(spltype[nbspl]=="manual")
            {
                nodes <- intnodes[(nb+1):(nb+nbnodes[nbspl]-2)]
                if(!length(nodes)) stop("The length of intnodes is not correct")
                intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <-  nodes
                nb <- nb+nbnodes[nbspl]-2
                
                idrg <- length(which(idlink[1:k] != 0))
                if(any(nodes <= range[2*(idrg-1)+1]) | any(nodes >= range[2*idrg])) stop("Interior nodes must be in the range of the outcome")
            }
            
            if(spltype[nbspl]=="equi")
            {
                nodes <- seq(range[2*(nbspl-1)+1], range[2*nbspl], length.out=nbnodes[nbspl])
                nodes <- nodes[-nbnodes[nbspl]]
                nodes <- nodes[-1]
                intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <- nodes
            }
            
            if(spltype[nbspl]=="quant")
            {
                nodes <- quantile(get(dataY[k])[,nomsY[k]], probs=seq(0,1,length.out=nbnodes[nbspl]))
                if(length(unique(nodes)) != length(nodes)) stop(paste("Some nodes are equal for link number",k,"; Please try to reduce the number of nodes or use manual location."))
                nodes <- nodes[-nbnodes[nbspl]]
                nodes <- nodes[-1]
                intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <- as.vector(nodes)
            }
        }
    }
    
    if(nb != length(intnodes)) stop(paste("The vector intnodes should be of length",nb))
    
    ##remplir zitr
    m <- 0
    if(nySPL>0) m <- max(nbnodes)
    zitr <- matrix(0,max(m,2),ny)
    nb12 <- 0
    nbspl <- 0
    for (k in 1:ny)
    {
        if((idlink[k]==0) | (idlink[k]==3)) zitr[1:2,k] <- c(min(get(dataY[k])[,nomsY[k]]),max(get(dataY[k])[,nomsY[k]]))
        
        if(idlink[k]==1)
        {
            nb12 <- nb12 + 1
            zitr[1:2,k] <- range[2*(nb12-1)+1:2]
        }
        
        if(idlink[k]==2)
        {
            nb12 <- nb12+1
            nbspl <- nbspl+1
            zitr[2:(nbzitr[k]-1),k] <- intnodes2[ifelse(nbspl==1,0,1)*sum(nbnodes[1:(nbspl-1)]-2) + 1:(nbnodes[nbspl]-2)]
            zitr[1,k] <- range[2*(nb12-1)+1]
            zitr[nbnodes[nbspl],k]  <- range[2*nb12]
            
            ##verifier s'il y a des obs entre les noeuds
            hcounts <- hist(get(dataY[k])[,nomsY[k]],breaks=zitr[1:nbnodes[nbspl],k],plot=FALSE,include.lowest=TRUE,right=TRUE)$counts
            if(any(hcounts==0)) stop(paste("Link function number",k,"can not be estimated. Please try other nodes such that there are observations in each interval."))
        }
    }
    
    ##uniqueY0 et indiceY0
    uniqueY0 <- NULL
    indiceY0 <- NULL
    nvalSPLORD <- rep(0,ny)
    nbmod <- rep(0,ny)
    modalites <- vector("list",ny)
    nb <- 0
    for (k in 1:ny)
    {
        if((idlink[k]!=2) & (idlink[k]!=3))
        {
            indiceY0 <- c(indiceY0, rep(0,length(get(dataY[k])[,nomsY[k]])))
            next
        }
        
        yk <- get(dataY[k])[,nomsY[k]]
        uniqueTemp <- sort(unique(yk))
        permut <- order(order(yk))  # sort(y)[order(order(y))] = y
        if(length(as.vector(table(yk)))==length(uniqueTemp))
        {
            indice <- rep(1:length(uniqueTemp), as.vector(table(yk)))
            if(idlink[k]==2)
            {
                indiceTemp <- nb + indice[permut]
            }
            else
            {
                indiceTemp <- indice[permut]
            }
            
            nb <- nb + length(uniqueTemp)
            
            uniqueY0 <- c(uniqueY0, uniqueTemp)
            indiceY0 <- c(indiceY0, indiceTemp)
            nvalSPLORD[k] <- length(uniqueTemp)
        }
        else
        {
            uniqueY0 <- c(uniqueY0, yk)
            indiceY0 <- c(indiceY0, ifelse(idlink[k]==2,nb,0)+c(1:length(yk)))
            nb <- nb + length(yk)
            nvalSPLORD[k] <- length(yk)
        }
        if(idlink[k]==3)
        {
            nbmod[k] <- length(na.omit(uniqueTemp))
            modalites[[k]] <- uniqueTemp
        }
    }
    
    
    ##ordonner les mesures par individu
    matYX <- cbind(IND,Y0,indiceY0,outcome,X0)
    matYXord <- matYX[order(IND),]
    Y0 <- as.numeric(matYXord[,2])
    X0 <- apply(matYXord[,-c(1:4),drop=FALSE],2,as.numeric)
    IND <- matYXord[,1]
    outcome <- matYXord[,4]
    indiceY0 <- as.numeric(matYXord[,3])
    
    dataSurv <- dataSurv[which(dataSurv[,1] %in% IND),]
    dataSurv <- dataSurv[order(dataSurv[,1]),]
    nmes <- as.vector(table(dataSurv[,1]))
    data.surv <- apply(dataSurv[cumsum(nmes),-1],2,as.numeric)
    tsurv0 <- data.surv[,1]
    tsurv <- data.surv[,2]
    devt <- data.surv[,3]
    tsurvint <- data.surv[,4]
    ind_survint <- (tsurvint<tsurv) + 0
    
    ## test de hazard
    arghaz <- hazard
    hazard <- rep(hazard,length.out=nbevt)
    if(any(hazard %in% c("splines","Splines")))
    {
        hazard[which(hazard %in% c("splines","Splines"))] <- "5-quant-splines"
    }
    if(any(hazard %in% c("piecewise","Piecewise")))
    {
        hazard[which(hazard %in% c("piecewise","Piecewise"))] <- "5-quant-piecewise"
    }
    
    haz13 <- strsplit(hazard[which(!(hazard=="Weibull"))],"-")
    if(any(sapply(haz13,length)!=3)) stop("Invalid argument hazard")
    
    nz <- rep(2,nbevt)
    locnodes <- NULL
    typrisq <- rep(2,nbevt)
    nprisq <- rep(2,nbevt)
    
    nznodes <- 0 #longueur de hazardnodes
    ii <- 0
    if(any(hazard!="Weibull"))
    {
        
        for (i in 1:nbevt)
        {
            if(hazard[i]=="Weibull") next;
            
            ii <- ii+1
            
            nz[i] <- as.numeric(haz13[[ii]][1])
            if(nz[i]<3) stop("At least 3 nodes are required")
            typrisq[i] <- ifelse(haz13[[ii]][3] %in% c("splines","Splines"),3,1)
            nprisq[i] <- ifelse(haz13[[ii]][3] %in% c("splines","Splines"),nz[i]+2,nz[i]-1)
            locnodes <- c(locnodes, haz13[[ii]][2])
            if(!(haz13[[ii]][3] %in% c("splines","Splines","piecewise","Piecewise"))) stop("Invalid argument hazard")
            
            if((haz13[[ii]][2]=="manual"))
            {
                nznodes <- nznodes + nz[i]-2
            }
            
            if(!all(locnodes %in% c("equi","quant","manual"))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")
        }
        
        if(!is.null(hazardnodes))
        {
            if(!any(locnodes == "manual"))  stop("hazardnodes should be NULL if the nodes are not chosen manually")
            
            if(length(hazardnodes) != nznodes) stop(paste("Vector hazardnodes should be of length",nznodes))
        }
    }
    else
    {
        if(!is.null(hazardnodes)) stop("hazardnodes should be NULL if Weibull baseline risk functions are chosen")
    }
    
    
    if(nbevt>1 & length(arghaz)==1 & nznodes>0)
    {
        hazardnodes <- rep(hazardnodes,length.out=nznodes*nbevt)
    }
    
    nrisqtot <- sum(nprisq)
    
    zi <- matrix(0,nrow=max(nz),ncol=nbevt)
    nb <- 0
    
    minT1 <- 0
    maxT1 <- max(tsurv)
    tsurvevt <- tsurv
    
    if(idtrunc==1)
    {
        minT1 <- min(tsurv,tsurv0)
        maxT1 <- max(tsurv,tsurv0)
    }
    
    ## arrondir
    minT2 <- round(minT1,3)
    if(minT1<minT2) minT2 <- minT2-0.001
    minT <- minT2
    
    maxT2 <- round(maxT1,3)
    if(maxT1>maxT2) maxT2 <- maxT2+0.001
    maxT <- maxT2
    
    if(length(hazardrange)){
      if(hazardrange[1]>minT) stop(paste("hazardrange[1] should be <=",minT))
      if(hazardrange[2]>maxT) stop(paste("hazardrange[2] should be >=",maxT))
      minT <- hazardrange[1]
      maxT <- hazardrange[2]
    }
    
    startWeib <- rep(0,nbevt)
    startWeib[which(typrisq==2)] <- rep(startWeibull, length.out=length(which(typrisq==2)))
    ii <- 0
    for(i in 1:nbevt)
    {
        if(typrisq[i]==2)
        {
            if(minT < startWeib[i]) stop("Some entry or event times are bellow startWeibull")
            zi[1:2,i] <- c(startWeib[i],maxT)
        }
        else
        {
            ii <- ii+1
            
            if(locnodes[ii]=="manual")
            {
                zi[1:nz[i],i] <- c(minT,hazardnodes[nb+1:(nz[i]-2)],maxT)
                nb <- nb + nz[i]-2
            }
            if(locnodes[ii]=="equi")
            {
                zi[1:nz[i],i] <- seq(minT,maxT,length.out=nz[i])
            }
            if(locnodes[ii]=="quant")
            {
                pi <- c(1:(nz[i]-2))/(nz[i]-1)
                qi <- quantile(tsurvevt,prob=pi)
                zi[1,i] <- minT
                zi[2:(nz[i]-1),i] <- qi
                zi[nz[i],i] <- maxT
            }
        }
    }
    
    
    ###TS: cas ou gi(bi,t)=niv.courant -> construction des matrices design pr predire lambda a chq pnt de quadrature GK
    if(sharedtype == 'RE'){  #gi(bi,t)=bi
        Xpred <- 0
        Xpred_Ti <- 0
        nbXpred <- 0
    }
    if(sharedtype == 'CL'){   #gi(bi,t)=niv.courant(t)
        
        # /!\ si varexp dependante du tps (autre que var.time), prediction impossible
        nom.var <- c(attr(terms(fixed2[-2]), "term.labels"),attr(terms(random), "term.labels"))  # noms des covariables EF et EA
        nom.var <- nom.var[!str_detect(nom.var,var.time)] # nom.var[!(nom.var == var.time)] # noms des variables, autre que celle incluant var.time
        if(length(nom.var)>0){
            for(v in 1:length(nom.var)){ #verif: covariables (hors var.time) independantes du temps
                tmp <- unique(na.omit(data[,c(subject,nom.var[v])]))  #dataframe 2 colonnes : subject et var v, en supprimant les lignes doublons
                if(nrow(tmp) != length(unique(IND))) #var v dependante du temps
                    stop(paste(nom.var[v]," variable seems to be time dependant, can't use sharedtype='CL' due to impossibility to predict"))
            }  
        }
        
        ## matrices design
        
        # database contenant les variables des EFs et des EAs, 1 ligne par sujet
        data_tmp <- as.data.frame(unique(na.omit(data[,c(subject,nom.var)])))
        colnames(data_tmp) <- c(subject,nom.var)
        data_tmp$tmp_name <- NA
        colnames(data_tmp)[which(colnames(data_tmp) == "tmp_name")] <- var.time
        if(nrow(data_tmp)!=length(unique(IND))) stop("Make sure at least 1 visit per subject has complete data")  # si un sujet a un valeur NA a 1 covar a chq visit, alors il n'est plus pris en compte
        if(nrow(data_tmp)!=length(tsurv)) stop("pblm")
        
        # Xpredcl_vectTi (Ti event) vecteur par patient
        data_tmp[,var.time] <- tsurv # remplacer la variable de temps par T
        #matrices design
        mat_ef <- model.matrix(object = fixed2[-2], data = data_tmp) # EF sans contraste et sans outcome, 
        mat_ef <- as.data.frame(mat_ef[,-1])  # et sans intercept car non-estime
        mat_ea <- model.matrix(object = random, data = data_tmp) # EA
        # concatenation
        Xpredcl_vectTi <- cbind(data_tmp[,var.time],mat_ef,mat_ea)  # 1ere colonne = tps de quadrature
        #nb_ind_Xpred <- c(ncol(mat_ef),ncol(mat_ea))  #nb de var des EF, EA
        Xpredcl_vectTi <- as.matrix(Xpredcl_vectTi)
        Xpred_Ti <- Xpredcl_vectTi
        nbXpred <- ncol(Xpred_Ti)
        
        # noeuds de quadrature Gauss-Kronrod
        ptGK_1 <- 0.991455371120812639206854697526329
        ptGK_2 <- 0.949107912342758524526189684047851
        ptGK_3 <- 0.864864423359769072789712788640926
        ptGK_4 <- 0.741531185599394439863864773280788
        ptGK_5 <- 0.586087235467691130294144838258730
        ptGK_6 <- 0.405845151377397166906606412076961
        ptGK_7 <- 0.207784955007898467600689403773245
        ptGK_8 <- 0.000000000000000000000000000000000
        
        # Xpredcl (Ti event)
        data_tmp[,var.time] <- tsurv # remplacer la variable de temps par T
        data_tmp[,var.time] <- (data_tmp[,var.time] - 0) / 2  #centr = centre de l intervalle 0 -> Ti
        data_tmp$mid_interv_surv <- data_tmp[,var.time]  #hlgth = demi longueur de l intervalle 0 -> Ti
        data_tmp <- data_tmp[rep(1:nrow(data_tmp),each=15),]  #1ligne/point quadrature/patient
        data_tmp$ptGK <- c(ptGK_1,-ptGK_1,ptGK_2,-ptGK_2,ptGK_3,-ptGK_3,ptGK_4,-ptGK_4,ptGK_5,-ptGK_5,ptGK_6,-ptGK_6,ptGK_7,-ptGK_7,ptGK_8) #pnts de quadrature GK
        data_tmp[,var.time] <- data_tmp[,var.time] + data_tmp$mid_interv_surv * data_tmp$ptGK  #centr+absc ou absc=hlgth*pnt
        #matrices design
        mat_ef <- model.matrix(object = fixed2[-2], data = data_tmp) # EF sans contraste et sans outcome,
        mat_ef <- as.data.frame(mat_ef[,-1])  # et sans intercept car non-estime
        mat_ea <- model.matrix(object = random, data = data_tmp) # EA
        # concatenation
        Xpredcl <- cbind(data_tmp[,var.time],mat_ef,mat_ea)  # 1ere colonne = tps de quadrature
        Xpredcl <- as.matrix(Xpredcl)
        Xpred <- Xpredcl
        
        # Xpredcl0 (Ti0, tronc gche)
        if(idtrunc==1){
            data_tmp[,var.time] <- rep(tsurv0,each=15) # remplacer la variable de temps par T0
            data_tmp[,var.time] <- (data_tmp[,var.time] - 0) / 2  #centr = centre de l intervalle 0 -> T0i
            data_tmp$mid_interv_surv <- data_tmp[,var.time]  #hlgth = demi longueur de l intervalle 0 -> T0i
            data_tmp[,var.time] <- data_tmp[,var.time] + data_tmp$mid_interv_surv * data_tmp$ptGK  #centr+absc ou absc=hlgth*pnt
            #matrices design
            mat_ef <- model.matrix(object = fixed2[-2], data = data_tmp) # EF sans contraste et sans outcome,
            mat_ef <- as.data.frame(mat_ef[,-1])  # et sans intercept car non-estime
            mat_ea <- model.matrix(object = random, data = data_tmp) # EA
            # concatenation
            Xpredcl0 <- cbind(data_tmp[,var.time],mat_ef,mat_ea)
            Xpredcl0 <- as.matrix(Xpredcl0)
            Xpred <- cbind(Xpred,Xpredcl0)
        }else{
            Xpredcl0 <- NULL
        }
        
        nbXpred <- c(nbXpred, ncol(Xpred))
        
    }
    
    
    ##parametres pour Fortran
    ns <- length(unique(IND))
    nv <- dim(X0)[2]
    nobs <- length(Y0)
    idiag0 <- ifelse(idiag==TRUE,1,0)
    nalea <- ifelse(randomY==TRUE,ny,0)
    logspecif <- as.numeric(logscale)
    #loglik <- 0
    ni <- 0
    istop <- 0
    gconv <- rep(0,3)
    resid_m <- rep(0,nobs)
    resid_ss <- rep(0,nobs)
    Yobs <- rep(0,nobs)
    time <- seq(minT,maxT,length.out=nsim)
    risq_est <- matrix(0,nrow=nsim,ncol=nbevt)
    risqcum_est <- matrix(0,nrow=nsim,ncol=nbevt)
    predRE_Y <- rep(0,ns*nalea)
    rlindiv <- rep(0,ns)
    marker <- rep(0,nsim*ny)
    transfY <- rep(0,nsim*ny)
    
    
    
    ##nmes
    nmes <- matrix(0,ns,ny)
    for (k in 1:ny)
    {
        INDpresents <- which(unique(IND) %in% get(dataY[k])[,nom.subject])
        nmes[INDpresents,k] <- as.vector(table(get(dataY[k])[,nom.subject]))
    }
    maxmes <- max(apply(nmes,1,sum))
    
    
    
    ##remplir idprob, etc
    z.X0 <- strsplit(nom.unique,split=":",fixed=TRUE)
    z.X0 <- lapply(z.X0,sort)
    
    idea <- (z.X0 %in% z.random) + 0
    idg <- (z.X0 %in% z.fixed) + 0
    idcontr <- (z.X0 %in% z.contr) + 0
    
    if (ncor>0) idcor <- colnames(X0) %in% cor.var.time +0
    else idcor <- rep(0,nv)
    
    idsurv <- z.X0 %in% z.surv + 2*(z.X0 %in% z.survcause)
    idsurv[1] <- 0 # 0 pour l'intercept
    
    idtdv <- z.X0 %in% nom.timedepvar + 0
    
    ## Si pas TimeDepVar dans formule survival
    if(length(nom.timedepvar) & all(idtdv==0))
    {
        stop("Variable in 'TimeDepVar' should also appear as a covariate in the 'survival' argument")
    }
    
    ## nb coef ppour survie
    nvarxevt <- sum(idsurv==1) + nbevt*sum(idsurv==2)
    
    nea <- sum(idea)
    predRE <- rep(0,ns*nea)
    
    ## nb parametres d'association  #TS
    if(sharedtype == 'RE')
        nasso <- nea*nbevt
    if(sharedtype == 'CL')
        nasso <- nbevt
    
    ## prm partie long
    nef <- sum(idg==1)-1
    ncontr <- (ny-1)*sum(idcontr)
    nvc <- ifelse(idiag0==1,nea,nea*(nea+1)/2)-1
    ntr <- rep(0,ny)
    ntr[which(idlink==0)] <- 2
    ntr[which(idlink==1)] <- 4
    ntr[which(idlink==2)] <- nbzitr[which(idlink==2)]+2
    ntr[which(idlink==3)] <- nbmod[which(idlink==3)]-1
    ntrtot <- sum(ntr)
    
    ##nombre total de parametres
    NPM <- nrisqtot + nvarxevt + nasso +
        nef + ncontr + nvc + ncor + ntrtot + nalea + ny
    
    
    V <- rep(0, NPM*(NPM+1)/2)  #pr variance des parametres
    
    ## parametres MC
    methInteg <- switch(methInteg,"MCO"=1,"MCA"=2,"QMC"=3)
    seqMC <- 0
    dimMC <- 0
    if(methInteg==3)
    {
        dimMC <- nea+nalea
        if(ncor>0) dimMC <- dimMC+maxmes
        if(dimMC>0) seqMC <- randtoolbox::sobol(n=nMC,dim=dimMC,normal=TRUE,scrambling=1)
    }
    
    

    ###valeurs initiales
    if(!(missing(B)))
    {
        if(!is.vector(B)) stop("B should be a vector")
        
        if (length(B)==NPM) b <- B
        else stop(paste("Vector B should be of length",NPM))
        
    }
    else ## B missing
    {
        b <- rep(0,NPM)
        
        if(nbevt>0){  #TS : si Weibull et prms au carre pr positivite -> valeurs par defaut = 1
            if(any(hazard!="Weibull")==FALSE & isFALSE(logscale)==TRUE){
                for(i in 1:nbevt){
                    if(typrisq[i]==2){
                        b[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- 1
                        if(idtrunc==1 & any(Tentry==0)) #TS : si entree retardee et au moins un temps vaut 0
                          b[sum(nprisq[1:i])-nprisq[i]+nprisq[i]] <- 1.25 #sinon pblm lors de recherche prms, risq instant tend vers infini qd w2 < 1 et t=0
                    }
                }
            }
            if(any(hazard %in% c("splines","Splines")) & idtrunc==1 & any(Tentry==0)){
                for(i in 1:nbevt){
                    if(typrisq[i]==3){
                        b[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- 10^-7 #sinon pgrm crash
                    }
                }
            }
        }
        
        if (nvc>0)
        {
            if(idiag==1) b[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- rep(1,nvc)
            if(idiag==0)
            {
                init.nvc <- diag(nea)
                init.nvc <- init.nvc[upper.tri(init.nvc, diag=TRUE)]
                b[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- init.nvc[-1]
            }
        }
        
        if(ncor>0) b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor] <- 1
        
        if(nalea>0) b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+1:nalea] <- 1
        
        b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+1:ny] <-  1
        
        for(k in 1:ny)
        {
            if(idlink[k]==0)
            {
                b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-1] <- mean(get(dataY[k])[,nomsY[k]])
                b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])] <- 1
            }
            if(idlink[k]==1)
            {
                b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-3] <- 0
                b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-2] <- -log(2)
                b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-1] <- 0.7
                b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])] <- 0.1
            }
            if(idlink[k]==2)
            {
                b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-ntr[k]+1] <- -2
                b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-ntr[k]+2:ntr[k]] <- 0.1
            }
            if(idlink[k]==3)
            {
                b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-ntr[k]+1] <- 0
                if(ntr[k]>1) b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-ntr[k]+2:ntr[k]] <- 0.1
            }
        }
    }
    
    
    
    
    ## ## faire wRandom et b0Random
    ## nef2 <- sum(idg!=0)-1 + (ny-1)*sum(idcontr)
    ## NPM2 <- nef2+ nvc+ncor+nalea+ny+ntrtot
    
    ## wRandom <- rep(0,NPM)
    ## b0Random <- rep(0,ng-1)
    
    ## l <- 0
    ## t <- 0
    ## m <- 0
    ## for (i in 1:nv)
    ## {
    ##     if(idg[i]==1)
    ##     {
    ##         if(i==1) next
    ##         l <- l+1
    ##         t <- t+1
    ##         wRandom[nprob+t] <- l
    ##     }
    ##     if(idg[i]==2)
    ##     {
    ##         if(i==1)
    ##         {
    ##             t <- t+ng-1
    ##             b0Random <- c(b0Random,rep(0,ng-1))
    ##             next
    ##         }
    ##         l <- l+1
    ##         for (g in 1:ng)
    ##         {
    ##             t <- t+1
    ##             wRandom[nprob+t] <- l
    ##         }
    ##     }
    ##     if(idcontr[i]==1)
    ##     {
    ##         wRandom[nef-ncontr+m+1:(ny-1)] <- nef2-ncontr+m+1:(ny-1)
    ##         m <- m+ny-1
    ##     }
    ## }
    
    ## if(nvc>0)
    ## {
    ##     wRandom[nef+1:nvc] <- nef2+1:nvc
    ## }
    ## if(nw>0)
    ## {
    ##     b0Random <- c(b0Random,rep(1,ng-1))
    ## }
    
    ## if(ncor>0) wRandom[nef+nvc+nw+1:ncor] <- nef2+nvc+1:ncor
    
    ## wRandom[nef+nvc+nw+ncor+1:ny] <- nef2+nvc+ncor+1:ny
    
    ## if(nalea>0)
    ## {
    ##     wRandom[nef+nvc+nw+ncor+ny+1:nalea] <- nef2+nvc+ncor+ny+1:nalea
    ## }
    
    ## wRandom[nef+nvc+nw+ncor+ny+nalea+1:ntrtot] <- nef2+nvc+ncor+ny+nalea+1:ntrtot
    ## ## wRandom et b0Random ok.
    
    
    ##------------------------------------------
    ##------nom au vecteur best
    ##--------------------------------------------
    
    nom.X0 <- colnames(X0)
    nom.X0[nom.X0=="(Intercept)"] <- "intercept"
    
    if(nbevt>0)
    {
        ##prm fct de risque
        if(isTRUE(logscale))
        {
            for(i in 1:nbevt)
            {
                nom1 <- rep(paste("event",i,sep=""),nprisq[i])
                if(typrisq[i]==2)
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:2]," log(Weibull",1:2,")",sep="")
                }
                if(typrisq[i]==1)
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," log(piecewise",1:(nz[i]-1),")",sep="")
                }
                if(typrisq[i]==3)
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," log(splines",1:(nz[i]+2),")",sep="")
                }
            }
        }
        else
        {
            for(i in 1:nbevt)
            {
                nom1 <- rep(paste("event",i,sep=""),nprisq[i])
                if(typrisq[i]==2)
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:2]," +/-sqrt(Weibull",1:2,")",sep="")
                }
                if(typrisq[i]==1)
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(piecewise",1:(nz[i]-1),")",sep="")
                }
                if(typrisq[i]==3)
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(splines",1:(nz[i]+2),")",sep="")
                }
            }
        }
        
        
        ##prm covariables survival
        nom1 <- NULL
        for(j in 1:nv)
        {
            if(idsurv[j]==1) #X
            {
                if(idtdv[j]==1)
                {
                    nom1 <- c(nom1,paste("I(t>",nom.timedepvar,")",sep=""))
                }
                else
                {
                    nom1 <- c(nom1,nom.X0[j])
                }
            }
            
            if(idsurv[j]==2) #cause(X)
            {
                if(idtdv[j]==1)
                {
                    nom1 <- c(nom1,paste("I(t>",nom.timedepvar,") event",1:nbevt,sep=""))
                }
                else
                {
                    nom1 <- c(nom1,paste(nom.X0[j],paste("event",1:nbevt,sep="")))
                }
            }
            
        }
        
        if(nvarxevt>0) names(b)[nrisqtot+1:nvarxevt] <- nom1
        
        
        if(sharedtype == 'RE'){   #TS
            for(i in 1:nbevt)
            {
                names(b)[nrisqtot+nvarxevt+(nbevt-1)*nea+1:nea] <- paste("event",i," asso",1:nea,sep="")
            }
        }
        if(sharedtype == 'CL')
            names(b)[nrisqtot+nvarxevt+1:nbevt] <- paste("event",1:nbevt," asso",sep="")
        
        
    }
    
    
    names(b)[nrisqtot+nvarxevt+nasso+1:nef] <- nom.X0[-1][idg[-1]!=0]
    if(ncontr!=0) names(b)[nrisqtot+nvarxevt+nasso+nef+1:ncontr] <- paste("contrast",paste(rep(1:sum(idcontr),each=ny-1),rep(1:(ny-1),sum(idcontr)),sep=""),sep="")
    
    if(idlink[1]==0) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntr[1]]<- c("Linear 1","Linear 2")
    if(idlink[1]==1) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntr[1]]<- paste("Beta",c(1:ntr[1]),sep="")
    if(idlink[1]==2) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntr[1]]<- paste("I-splines",c(1:ntr[1]),sep="")
    if(idlink[1]==3) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntr[1]]<- paste("Thresh",c(1:ntr[1]),sep="")
    if(ny>1)
    {
        for (yk in 2:ny)
        {
            if(idlink[yk]==0) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:(yk-1)])+1:ntr[yk]]<- c("Linear 1","Linear 2")
            if(idlink[yk]==1) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:(yk-1)])+1:ntr[yk]]<- paste("Beta",c(1:ntr[yk]),sep="")
            if(idlink[yk]==2) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:(yk-1)])+1:ntr[yk]]<- paste("I-splines",c(1:ntr[yk]),sep="")
            if(idlink[yk]==3) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:(yk-1)])+1:ntr[yk]]<- paste("Thresh",c(1:ntr[yk]),sep="")
        }
    }
    if(nvc!=0)names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- paste("varcov",c(1:nvc))
    
    if(ncor>0) {names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+1:ncor] <- paste("cor",1:ncor,sep="")}
    if(nalea!=0) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+1:nalea] <- paste("std.randomY",1:ny,sep="")
    
    names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+1:ny] <- paste("std.err",1:ny)
    namesb <- names(b)
    
    
    
    ## prm fixes
    fix <- rep(0,NPM)
    if(length(posfix))
    {
        if(any(!(posfix %in% 1:NPM))) stop("Indexes in posfix are not correct")
        
        fix[posfix] <- 1
    }
    if(length(posfix)==NPM) stop("No parameter to estimate")
    
    if(!all(partialH %in% 1:NPM)) stop(paste("partialH should contain indices between 1 and",NPM))
    
    # varcov RE
    if(nvc!=0)
    {
        mvc <- matrix(0,nea,nea)
        if(idiag0==0)
        {
            mvc[upper.tri(mvc, diag=TRUE)] <- c(1,b[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc])
            mvc <- t(mvc)
            mvc[upper.tri(mvc, diag=TRUE)] <- c(1,b[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc])
            ch <- chol(mvc)
            
            b[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- ch[upper.tri(ch, diag=TRUE)][-1]
        }
        else
        {
            b[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- sqrt(b[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc])
        }            
    }
    
    # fixed prms
    nfix <- sum(fix)
    bfix <- 0
    if(nfix>0)
    {
        bfix <- b[which(fix==1)]
        b <- b[which(fix==0)]
        NPM <- NPM-nfix
    }
    
    
    
    #browser()
    ###estimation
    
    conv3 <- c(convB, convL, convG) #TS: pr reduire le nb d arguments ds appel fct Fortran
    
    sharedtype <- ifelse(sharedtype == 'RE',1,2) #recodage 1 pr RE, 2 pr CL  #TS
    
    
    
    ###############
    ###   MLA   ###
    ###############
    
    if(maxiter==0)
    {
        vrais <- loglik(b,Y0,X0,tsurv0,tsurv,devt,ind_survint,idea,idg,idcor,idcontr,
                        idsurv,idtdv,typrisq,nz,zi,nbevt,idtrunc,logspecif,ny,ns,nv,
                        nobs,nmes,idiag0,ncor,nalea,NPM,nfix,bfix,epsY,
                        idlink,nbzitr,zitr,uniqueY0,indiceY0,nvalSPLORD,fix,
                        methInteg,nMC,dimMC,seqMC,sharedtype,nbXpred,Xpred_Ti,Xpred)
        
        out <- list(conv=2, V=rep(NA, length(b)), best=b, predRE=NA, predRE_Y=NA,
                    Yobs=NA, resid_m=NA, resid_ss=NA, risqcum_est=NA, risq_est=NA,
                    marker=NA, transfY=NA, gconv=rep(NA,3), niter=0, loglik=vrais)
    }
    else
    {   
        res <- mla(b=b, m=length(b), fn=loglik,
                   clustertype=clustertype,.packages=NULL,
                   epsa=convB,epsb=convL,epsd=convG,
                   digits=8,print.info=verbose,blinding=FALSE,
                   multipleTry=25,file="",partialH=partialH,
                   nproc=nproc,maxiter=maxiter,minimize=FALSE,
                   Y0=Y0,X0=X0,Tentr0=tsurv0,Tevt0=tsurv,Devt0=devt,
                   ind_survint0=ind_survint,idea0=idea,idg0=idg,idcor0=idcor,
                   idcontr0=idcontr,idsurv0=idsurv,idtdv0=idtdv,typrisq0=typrisq,
                   nz0=nz,zi0=zi,nbevt0=nbevt,idtrunc0=idtrunc,logspecif0=logspecif,
                   ny0=ny,ns0=ns,nv0=nv,nobs0=nobs,nmes0=nmes,idiag0=idiag0,
                   ncor0=ncor,nalea0=nalea,npm0=NPM,nfix0=nfix,bfix0=bfix,
                   epsY0=epsY,idlink0=idlink,nbzitr0=nbzitr,zitr0=zitr,
                   uniqueY0=uniqueY0,indiceY0=indiceY0,nvalSPLORD0=nvalSPLORD,
                   fix0=fix,methInteg0=methInteg,nMC0=nMC,dimMC0=dimMC,seqMC0=seqMC,
                   idst0=sharedtype,nXcl0=nbXpred,Xcl_Ti0=Xpred_Ti,Xcl_GK0=Xpred)
        
        out <- list(conv=res$istop, V=res$v, best=res$b, predRE=NA, predRE_Y=NA, Yobs=NA,
                    resid_m=NA, resid_ss=NA, risqcum_est=NA, risq_est=NA, marker=NA,
                    transfY=NA, gconv=c(res$ca, res$cb, res$rdm), niter=res$ni,
                    loglik=res$fn.value)
    }
    
    
    
    ## creer best a partir de b et bfix
    best <- rep(NA,length(fix))
    best[which(fix==0)] <- out$best
    best[which(fix==1)] <- bfix
    out$best <- best
    NPM <- NPM+nfix
        

   
    ## mettre NA pour les variances et covariances non calculees et  0 pr les prm fixes
    if(length(posfix))
    {
        mr <- NPM-length(posfix)
        Vr <- matrix(0,mr,mr)
        Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
        Vr <- t(Vr)
        Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
        V <- matrix(0,NPM,NPM)
        V[setdiff(1:NPM,posfix),setdiff(1:NPM,posfix)] <- Vr
        V <- V[upper.tri(V,diag=TRUE)]
    }
    else
    {
        V <- out$V
    }
    
    
    ## Creation du vecteur cholesky
    Cholesky <- rep(0,(nea*(nea+1)/2))
    if(idiag0==0 & nvc>0)
    {
        Cholesky[1:(nvc+1)] <- c(1,out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc])
        ## Construction de la matrice U
        U <- matrix(0,nrow=nea,ncol=nea)
        U[upper.tri(U,diag=TRUE)] <- Cholesky[1:(nvc+1)]
        z <- t(U) %*% U
        out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- z[upper.tri(z,diag=TRUE)][-1]
    }
    if(idiag0==1 & nvc>0)
    {
        id <- 1:nea
        indice <- rep(id+id*(id-1)/2)
        Cholesky[indice] <- c(1,out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc])
        out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc]**2
    }
    
    ##predictions
    predRE <- matrix(out$predRE,ncol=nea,byrow=T)
    predRE <- data.frame(unique(IND),predRE)
    colnames(predRE) <- c(nom.subject,nom.X0[idea!=0])
    
    if (nalea!=0)
    {
        predRE_Y <- matrix(out$predRE_Y,ncol=ny,byrow=TRUE)
        predRE_Y <- data.frame(unique(IND),predRE_Y)
        colnames(predRE_Y)  <- c(nom.subject,nomsY)
    }
    else
    {
        predRE_Y <- rep(NA,nalea*ns)
    }
    
    
    ##pred
    pred_m <- out$Yobs-out$resid_m
    pred_ss <- out$Yobs - out$resid_ss
    pred <- data.frame(IND,outcome,pred_m,out$resid_m,pred_ss,out$resid_ss,out$Yobs)
    
    colnames(pred)<-c(nom.subject,"Yname","pred_m","resid_m","pred_ss","resid_ss","obs")
    rownames(pred) <- NULL
    
    
    ## risques
    if(nbevt>0)
    {
        risqcum_est <- matrix(out$risqcum_est,nrow=nsim,ncol=nbevt)
        risq_est <- matrix(out$risq_est,nrow=nsim,ncol=nbevt)
        predSurv <- cbind(time,risq_est,risqcum_est)
        
        temp <- paste("event",1:nbevt,".RiskFct",sep="")
        temp1 <- paste("event",1:nbevt,".CumRiskFct",sep="")
        colnames(predSurv) <- c("time",temp,temp1)
        rownames(predSurv) <- 1:nsim
    }
    else
    {
        predSurv <- NA
    }
    
    ##estimlink
    ysim <- matrix(out$marker,nsim,ny)
    transfo <- matrix(out$transfY,nsim,ny)
    estimlink <- as.vector(rbind(ysim,transfo))
    estimlink <- matrix(estimlink,nsim,2*ny)
    colnames(estimlink) <- paste(c("","transf"),rep(nomsY, each=2),sep="")
    
    if(any(idlink==3))
    {
        if(any(2*nbmod[which(idlink==3)] > nsim))
        {
            nsim2 <- max(2*nbmod[which(idlink==3)])
            
            estimlink2 <- matrix(NA, nrow=nsim2, ncol=2*ny)
            estimlink2[1:nsim,] <- estimlink
            
            colnames(estimlink2) <- colnames(estimlink)
            estimlink <- estimlink2
        }
        
        seuils <- function(x)
        {
            n <- length(x)
            if(n>1) res <- c(x[1],x[1]+cumsum(x[-1]^2))
            if(n==1) res <- x[1]
            return(res)
        }
        
        sumntr <- 0
        for (k in 1:ny)
        {
            if(idlink[k]==3)
            {
                estimlink[,2*(k-1)+1] <- NA
                estimlink[,2*(k-1)+2] <- NA
                
                nb <- nbmod[k]-1
                
                marker <- rep(modalites[[k]], each=2)
                seuilsk <- seuils(out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1:nb])
                transfY <- rep(seuilsk, each=2)
                transfY <- c(-Inf,transfY,Inf)
                
                estimlink[1:length(marker), 2*(k-1)+1] <- marker
                estimlink[1:length(transfY), 2*(k-1)+2] <- transfY
            }
            
            sumntr <- sumntr + ntr[k]
        }
        
        ## remplacer -Inf et Inf
        hy <- estimlink[,2*(1:ny)]
        maxhy <- max(hy[is.finite(hy)])
        minhy <- min(hy[is.finite(hy)])
        for (k in 1:ny)
        {
            if(idlink[k]==3)
            {
                estimlink[1, 2*(k-1)+2] <- minhy - (maxhy - minhy)/10
                estimlink[2*nbmod[k], 2*(k-1)+2] <- maxhy + (maxhy - minhy)/10
            }
        }
        
    }
    
    N <- rep(NA,12)
    N[1] <- 0 #nprob
    N[2] <- nrisqtot
    N[3] <- nvarxevt
    N[4] <- nasso
    N[5] <- nef
    N[6] <- ncontr
    N[7] <- nvc
    N[8] <- 0 #nw
    N[9] <- ncor
    N[10] <- ntrtot
    N[11] <- nalea
    N[12] <- ny
    N[13] <- nobs
    N[14] <- nbevt
    
    nevent <- rep(0,nbevt)
    for(ke in 1:nbevt)
    {
        nevent[ke] <- length(which(devt==ke))
    }
    
    Nprm <- c(nprisq,ntr)
    
    
    nom.X0[nom.X0=="(Intercept)"] <- "Intercept"
    
    ## noms des variables
    Names <- list(Xnames=nom.X0,Ynames=nomsY,
                  ID=nom.subject,Tnames=noms.surv,
                  TimeDepVar.name=nom.timedepvar,
                  Xvar=setdiff(ttesLesVar,noms.surv))
    
    names(modalites) <- nomsY
    
    form <- list(fixed=fixed2[-2], random=random, contr=contr,
                 form.commun=form.commun, form.cause=form.cause,
                 form.cor=form.cor)
    
    cost <- proc.time()-ptm
    
    res <-list(ns=ns,idg=idg,idcontr=idcontr,idea=idea,idcor=idcor,
               idsurv=idsurv,idtdv=idtdv,
               loglik=out$loglik,best=out$best,V=V,gconv=out$gconv,conv=out$conv,
               call=cl,niter=out$niter,N=N,nevent=nevent,Nprm=Nprm,
               idiag=idiag0,pred=pred,predRE=predRE,
               predRE_Y=predRE_Y,Names=Names,form=form,cholesky=Cholesky,
               logspecif=logspecif,predSurv=predSurv,typrisq=typrisq,hazardnodes=zi,nz=nz,
               estimlink=estimlink,epsY=epsY,linktype=idlink,linknodes=zitr,nbnodes=nbnodes,nbmod=nbmod,mod=modalites,
               na.action=nayk,AIC=2*(length(out$best)-length(posfix)-out$loglik),BIC=(length(out$best)-length(posfix))*log(ns)-2*out$loglik,data=datareturn,
               #wRandom=wRandom,b0Random=b0Random,
               posfix=posfix,CPUtime=cost[3])
    
    names(res$best) <- namesb
    class(res) <-c("jointLPM")
    
    
    return(res)
}

loglik <- function(b0,Y0,X0,Tentr0,Tevt0,Devt0,ind_survint0,idea0,idg0,idcor0,idcontr0,
                   idsurv0,idtdv0,typrisq0,nz0,zi0,nbevt0,idtrunc0,logspecif0,ny0,ns0,
                   nv0,nobs0,nmes0,idiag0,ncor0,nalea0,npm0,nfix0,bfix0,
                   epsY0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0,nvalSPLORD0,fix0,
                   methInteg0,nMC0,dimMC0,seqMC0,idst0,nXcl0,Xcl_Ti0,Xcl_GK0)
{
    res <- 0
    .Fortran(C_loglik,as.double(Y0),as.double(X0),as.double(Tentr0),as.double(Tevt0),as.integer(Devt0),as.integer(ind_survint0),as.integer(idea0),as.integer(idg0),as.integer(idcor0),as.integer(idcontr0),as.integer(idsurv0),as.integer(idtdv0),as.integer(typrisq0),as.integer(nz0),as.double(zi0),as.integer(nbevt0),as.integer(idtrunc0),as.integer(logspecif0),as.integer(ny0),as.integer(ns0),as.integer(nv0),as.integer(nobs0),as.integer(nmes0),as.integer(idiag0),as.integer(ncor0),as.integer(nalea0),as.integer(npm0),as.double(b0),as.integer(nfix0),as.double(bfix0),as.double(epsY0),as.integer(idlink0),as.integer(nbzitr0),as.double(zitr0),as.double(uniqueY0),as.integer(indiceY0),as.integer(nvalSPLORD0),as.integer(fix0),as.integer(methInteg0),as.integer(nMC0),as.integer(dimMC0),as.double(seqMC0),as.integer(idst0),as.integer(nXcl0),as.double(Xcl_Ti0),as.double(Xcl_GK0),loglik_res=as.double(res))$loglik_res
}

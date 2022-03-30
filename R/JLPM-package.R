#' Estimation of joint latent process models
#' 
#' Functions for the estimation of joint latent process models (JLPM).
#' Continuous and ordinal outcomes are handled for the longitudinal part,
#' whereas the survival part considers multiple competing events.
#' The likelihood is computed using Monte-Carlo integration. Estimation is achieved
#' by maximizing the log-likelihood using a robust iterative algorithm.
#'
#' Please report to the JLPM-team any question or suggestion regarding the package
#' via github only (https://github.com/VivianePhilipps/JLPM/issues).
#' 
#' @name JLPM-package
#' @docType package
#' @author Cecile Proust-Lima, Viviane Philipps, Tiphaine Saulnier
#' 
#' @references
#' Saulnier, Philipps, Meissner, Rascol, Pavy-Le-Traon, Foubert-Samier, Proust-Lima (2021).
#' Joint models for the longitudinal analysis of measurement scales in the presence 
#' of informative dropout, arXiv:2110.02612.
#' 
#' Philipps, Hejblum, Prague, Commenges, Proust-Lima (2021).
#' Robust and efficient optimization using a Marquardt-Levenberg algorithm with 
#' R package marqLevAlg, The R Journal 13:2.
#'
#' @keywords package
#' @importFrom graphics axis hist lines matlines matplot mtext par plot points segments polygon
#' @importFrom grDevices rainbow rgb col2rgb n2mfrow
#' @importFrom stats as.formula formula get_all_vars integrate median model.frame model.matrix na.fail na.omit na.pass pchisq pnorm qnorm quantile rnorm sd terms residuals vcov fitted coef update
#' @importFrom survival Surv untangle.specials
#' @importFrom randtoolbox sobol
#' @importFrom stringr str_detect
#' @importFrom parallel clusterEvalQ clusterExport clusterSetRNGStream makeCluster parApply stopCluster
#' @importFrom marqLevAlg mla
#' @useDynLib JLPM, .registration=TRUE, .fixes="C_"
NULL










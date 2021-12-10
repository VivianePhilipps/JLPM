#' Joint latent process models
#' 
#' Functions for the estimation of joint latent process models (JLPM)
#' 
#' @name JLPM-package
#' @docType package
#' @author Cecile Proust-Lima, Viviane Philipps, Tiphaine Saulnier
#' 
#' \email{cecile.proust-lima@@inserm.fr}
#' @references
#' 
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










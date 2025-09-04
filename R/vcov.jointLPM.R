#' @export
vcov.jointLPM <- function(object, ...)
{
    res <- matrix(0, length(object$best), length(object$best))

    res[upper.tri(res, diag = TRUE)] <- object$V
    res <- t(res)
    res[upper.tri(res, diag = TRUE)] <- object$V

    noms <- sub("varcov", "cholesky", names(object$best))
    rownames(res) <- colnames(res) <- noms

    return(res)
}

fn_solve <- function(v, A, tol   = sqrt(.Machine$double.eps), ridge = 0) {
  v <- as.numeric(v)
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A), length(v) == nrow(A))

  A <- 0.5 * (A + t(A))
  if (ridge > 0) A <- A + ridge * diag(nrow(A))

  ch <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(ch)) {
    y <- forwardsolve(t(ch), v)
    x <- backsolve(ch, y)          # x = A^{-1} v
    return(drop(t(v) %*% x))       # v' A^{-1} v
  }

  ee <- eigen(A, symmetric = TRUE)
  lam <- ee$values
  U   <- ee$vectors

  cutoff <- tol * max(1, max(abs(lam)))
  keep   <- abs(lam) > cutoff

  if (!any(keep)) {
    return(0)
  }

  proj <- drop(t(U[, keep, drop = FALSE]) %*% v)  # U^T v
  chi  <- sum( (proj * proj) / lam[keep] )        # v' A^{+} v

  return(chi)
}

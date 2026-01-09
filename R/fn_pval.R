fn_pval <- function(mhat, zva, norr, crn) {
  m1 <- mhat[, 1]
  m2 <- mhat[, 2]
  rho <- sum(m1 * m2) / sqrt(sum(m1^2) * sum(m2^2))
  cmat <- matrix(c(1, rho, rho, 1), 2, 2)
  svd_res <- svd(cmat)
  rsig <- svd_res$u %*% sqrt(diag(svd_res$d)) %*% t(svd_res$v)
  sa <- norr %*% rsig
  sa_max <- pmax(abs(sa[, 1]), abs(sa[, 2]))
  mxo <- max(abs(zva))
  pval <- mean(sa_max > mxo / crn)
  list(mxo = mxo, pval = pval)
}

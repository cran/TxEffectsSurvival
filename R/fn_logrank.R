fn_logrank <- function(y, d, z) {

  I <- order(y)
  oy <- y[I]
  od <- d[I]
  oz <- z[I]
  nn <- length(y)

  s1 <- rev(cumsum(rev(oz)))
  s0 <- rev(cumsum(rep(1, nn)))

  t <- sum(od * (oz - s1 / s0))
  v <- sum(od * (s1 / s0 - (s1 / s0)^2))
  lr <- t / sqrt(v)
  pvalr <- 2 * pnorm(abs(lr), lower.tail = FALSE)

  results <- list()
  results$stat_lr <- lr^2
  results$pval_lr <- pvalr

  results
}

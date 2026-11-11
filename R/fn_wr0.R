fn_wr0 <- function(yh, hcen, yd, dcen, zo) {

  loc <- order(zo)
  ydct <- yd[loc]
  dcenct <- dcen[loc]
  yhct <- yh[loc]
  hcenct <- hcen[loc]

  n <- length(zo)
  n2 <- sum(zo)
  n1 <- n - n2

  s1 <- matrix(1, nrow = n1, ncol = n2)
  s2 <- matrix(1, nrow = n1, ncol = n2)
  tdw <- s1; thw <- s1; tdl <- s1; thl <- s1

  for(ci in 1:n1){
    tdw[ci, ] <- dcenct[ci] * (ydct[ci] <= ydct[n1 + (1:n2)])
    thw[ci, ] <- hcenct[ci] * (yhct[ci] <= yhct[n1 + (1:n2)])
    tdl[ci, ] <- dcenct[n1 + (1:n2)] * (ydct[n1 + (1:n2)] <= ydct[ci])
    thl[ci, ] <- hcenct[n1 + (1:n2)] * (yhct[n1 + (1:n2)] <= yhct[ci])
  }

  s2 <- tdl + (1 - tdl) * (1 - tdw) * thl
  s1 <- tdw + (1 - tdl) * (1 - tdw) * thw

  s1i <- rowSums(s1)
  s2i <- rowSums(s2)
  s1j <- colSums(s1)
  s2j <- colSums(s2)

  win <- mean(s1)
  los <- mean(s2)

  wr <- los / win
  results <- list(wr = wr, ci_wr = rep(NA_real_, 2L),
                  zv_wr = NA_real_, pval_wr = NA_real_)
  if (win == 0 || los == 0) return(results)

  # Exclude products from the identical control-treatment pair.
  same_pair <- sum(s1 * s2)
  shared_control <- c(
    sum(s1i * (s1i - 1)),
    sum(s1i * s2i) - same_pair,
    sum(s2i * (s2i - 1))
  ) / (n1 * n2 * (n2 - 1))
  shared_treatment <- c(
    sum(s1j * (s1j - 1)),
    sum(s1j * s2j) - same_pair,
    sum(s2j * (s2j - 1))
  ) / (n1 * n2 * (n1 - 1))

  independent <- c(win^2, win * los, los^2)
  covariance <- (shared_control - independent) / n1 +
    (shared_treatment - independent) / n2
  log_variance <- sum(covariance * c(1 / win^2, -2 / (los * win), 1 / los^2))
  if (!is.finite(log_variance) || log_variance <= 0) return(results)

  se <- sqrt(log_variance)
  results$ci_wr <- wr * exp(c(-1, 1) * 1.96 * se)
  results$zv_wr <- log(wr) / se
  results$pval_wr <- 2 * pnorm(abs(results$zv_wr), lower.tail = FALSE)
  results

}


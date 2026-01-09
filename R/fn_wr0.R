fn_wr0 <- function(yh, hcen, yd, dcen, zo) {

  loc <- order(zo)
  zoct <- zo[loc]
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

  pwxyyp <- sum(s1i * (s1i - 1)) / (n1 * n2 * (n1 - 1))
  pwxxpy <- sum(s1j * (s1j - 1)) / (n1 * n2 * (n2 - 1))
  plxyyp <- sum(s2i * (s2i - 1)) / (n1 * n2 * (n1 - 1))
  plxxpy <- sum(s2j * (s2j - 1)) / (n1 * n2 * (n2 - 1))
  pmxyyp <- sum(s1i * (s2i - 1)) / (n1 * n2 * (n1 - 1))
  pmxxpy <- sum(s1j * (s2j - 1)) / (n1 * n2 * (n2 - 1))

  ts <- n2 / n
  cs <- 1 - ts

  xi1110 <- pwxyyp - win^2
  xi1210 <- pmxyyp - win * los
  xi2210 <- plxyyp - los^2

  xi1101 <- pwxxpy - win^2
  xi1201 <- pmxxpy - win * los
  xi2201 <- plxxpy - los^2

  sig11 <- xi1110 / ts + xi1101 / cs
  sig12 <- xi1210 / ts + xi1201 / cs
  sig22 <- xi2210 / ts + xi2201 / cs

  wr <- los / win
  siglog <- sig11 / win^2 + sig22 / los^2 - 2 * sig12 / (los * win)

  wci <- c(wr * exp(-1.96 * sqrt(siglog / n)), wr * exp(1.96 * sqrt(siglog / n)))

  zvawr <- log(wr) / sqrt(siglog / n)
  pvawr <- 2 * pnorm(abs(zvawr), lower.tail = FALSE)

  results <- list()
  results$wr <- wr
  results$ci_wr <- wci
  results$zv_wr <- zvawr
  results$pval_wr <- pvawr

  results

}


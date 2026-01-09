fn_tntinf <- function(yh, hcen, yd, dcen, zo, lin, norr, trp, trpn, ak) {

  n <- length(yh)
  crn <- 1 + 15 / max(c(n, 100))

  # yh, hcen
  yor <- yh
  dor <- hcen

  rk <- sapply(1:n, function(i) sum(yor[i] >= yor))
  ginfh <- fn_ginf(yor, dor, zo, trp, trpn, ak)

  mhat01 <- ginfh$mhat0[rk]
  mhatp101 <- ginfh$mhatp10[rk]
  mhatg01 <- ginfh$mhatg0[rk]
  mhatrm1 <- ginfh$mhatrm0[rk]
  wrall1 <- ginfh$wrall
  stall1 <- ginfh$stall
  zva1 <- ginfh$zva
  zvap11 <- ginfh$zvap1
  zval1 <- ginfh$zval
  zvarm1 <- ginfh$zvarm
  zvat1 <- ginfh$zvat
  zvalt1 <- ginfh$zvalt
  zvap1t1 <- ginfh$zvap1t
  zvarmt1 <- ginfh$zvarmt

  wrl <- ginfh$wrl
  stndl <- ginfh$stndl
  nwr <- ginfh$nwr
  stndrm <- ginfh$stndrm
  tqt <- ginfh$tqt
  stnd <- ginfh$stnd
  tdif <- ginfh$tdif
  tqtp1 <- ginfh$tqtp1
  stndp1 <- ginfh$stndp1
  tdifp1 <- ginfh$tdifp1
  tqtl <- ginfh$tqtl
  stndl <- ginfh$stndl
  tdifl <- ginfh$tdifl
  tqtn <- ginfh$tqtn
  stndrm <- ginfh$stndrm
  tdifn <- ginfh$tdifn

  ci1l <- c(wrl - 1.96 * stndl, wrl + 1.96 * stndl)
  ci1rm <- c(nwr - 1.96 * stndrm, nwr + 1.96 * stndrm)

  lef <- tqt - 1.96 * stnd * tdif * crn
  rit <- tqt + 1.96 * stnd * tdif * crn
  ci1t <- log(2 / (1 - exp(-exp(c(lef, rit)))) - 1)

  lef <- tqtp1 - 1.96 * stndp1 * tdifp1 * crn
  rit <- tqtp1 + 1.96 * stndp1 * tdifp1 * crn
  ci1p1t <- log(2 / (1 - exp(-exp(c(lef, rit)))) - 1)

  lef <- tqtl - 1.96 * stndl * tdifl * crn
  rit <- tqtl + 1.96 * stndl * tdifl * crn
  ci1lt <- log(2 / (1 - exp(-exp(c(lef, rit)))) - 1)

  lef <- tqtn - 1.96 * stndrm * tdifn * crn
  rit <- tqtn + 1.96 * stndrm * tdifn * crn
  ci1rmt <- log(2 / (1 - exp(-exp(c(lef, rit)))) - 1)

  # yd, dcen
  yor <- yd
  dor <- dcen
  rk <- sapply(1:n, function(i) sum(yor[i] >= yor))
  ginfd <- fn_ginf(yor, dor, zo, trp, trpn, ak)

  mhat02 <- ginfd$mhat0[rk]
  mhatp102 <- ginfd$mhatp10[rk]
  mhatg02 <- ginfd$mhatg0[rk]
  mhatrm2 <- ginfd$mhatrm0[rk]
  wrall2 <- ginfd$wrall
  stall2 <- ginfd$stall
  zva2 <- ginfd$zva
  zvap12 <- ginfd$zvap1
  zval2 <- ginfd$zval
  zvarm2 <- ginfd$zvarm
  zvat2 <- ginfd$zvat
  zvalt2 <- ginfd$zvalt
  zvap1t2 <- ginfd$zvap1t
  zvarmt2 <- ginfd$zvarmt

  wrl <- ginfd$wrl
  stndl <- ginfd$stndl
  nwr <- ginfd$nwr
  stndrm <- ginfd$stndrm
  tqt <- ginfd$tqt
  stnd <- ginfd$stnd
  tdif <- ginfd$tdif
  tqtp1 <- ginfd$tqtp1
  stndp1 <- ginfd$stndp1
  tdifp1 <- ginfd$tdifp1
  tqtl <- ginfd$tqtl
  stndl <- ginfd$stndl
  tdifl <- ginfd$tdifl
  tqtn <- ginfd$tqtn
  stndrm <- ginfd$stndrm
  tdifn <- ginfd$tdifn

  ci2l <- c(wrl - 1.96 * stndl, wrl + 1.96 * stndl)
  ci2rm <- c(nwr - 1.96 * stndrm, nwr + 1.96 * stndrm)

  lef <- tqt - 1.96 * stnd * tdif * crn
  rit <- tqt + 1.96 * stnd * tdif * crn
  ci2t <- log(2 / (1 - exp(-exp(c(lef, rit)))) - 1)

  lef <- tqtp1 - 1.96 * stndp1 * tdifp1 * crn
  rit <- tqtp1 + 1.96 * stndp1 * tdifp1 * crn
  ci2p1t <- log(2 / (1 - exp(-exp(c(lef, rit)))) - 1)

  lef <- tqtl - 1.96 * stndl * tdifl * crn
  rit <- tqtl + 1.96 * stndl * tdifl * crn
  ci2lt <- log(2 / (1 - exp(-exp(c(lef, rit)))) - 1)

  lef <- tqtn - 1.96 * stndrm * tdifn * crn
  rit <- tqtn + 1.96 * stndrm * tdifn * crn
  ci2rmt <- log(2 / (1 - exp(-exp(c(lef, rit)))) - 1)

  # all
  wralla <- rbind(wrall1, wrall2)
  stalla <- rbind(stall1, stall2)

  mhatall <- cbind(mhat01, mhat02)
  mhatallp1 <- cbind(mhatp101, mhatp102)
  mhatalll <- cbind(mhatg01, mhatg02)
  mhatallrm <- cbind(mhatrm1, mhatrm2)

  # pvalues
  m1 <- mhat01
  m2 <- mhat02
  z1 <- zva1
  z2 <- zva2

  calc_pval <- function(m1, m2, z1, z2) {
    rho <- sum(m1 * m2) / prod(sqrt(colSums(cbind(m1, m2)^2)))
    cmat <- matrix(c(1, rho, rho, 1), 2, 2)
    svd_res <- svd(cmat)
    rsig <- svd_res$u %*% sqrt(diag(svd_res$d)) %*% t(svd_res$v)
    sa <- norr %*% rsig
    mxo <- max(abs(c(z1, z2)))
    pval <- mean(apply(abs(sa), 1, max) > mxo/crn)

    list(mxo = mxo, pval = pval)
  }

  res.mx <- calc_pval(mhat01, mhat02, zva1, zva2)
  res.mxt <- calc_pval(mhat01, mhat02, zvat1, zvat2)
  res.mxp1 <- calc_pval(mhatp101, mhatp102, zvap11, zvap12)
  res.mxp1t <- calc_pval(mhatp101, mhatp102, zvap1t1, zvap1t2)
  res.mxl <- calc_pval(mhatg01, mhatg02, zval1, zval2)
  res.mxlt <- calc_pval(mhatg01, mhatg02, zvalt1, zvalt2)
  res.mxrm <- calc_pval(mhatrm1, mhatrm2, zvarm1, zvarm2)
  res.mxrmt <- calc_pval(mhatrm1, mhatrm2, zvarmt1, zvarmt2)

  mx <-res.mx$mxo
  mxt <-res.mxt$mxo
  mxp1 <- res.mxp1$mxo
  mxp1t <- res.mxp1t$mxo
  mxl <- res.mxl$mxo
  mxlt <- res.mxlt$mxo
  mxrm <- res.mxrm$mxo
  mxrmt <- res.mxrmt$mxo

  pvalmx <-res.mx$pval
  pvalmxt <-res.mxt$pval
  pvalmxp1 <- res.mxp1$pval
  pvalmxp1t <- res.mxp1t$pval
  pvalmxl <- res.mxl$pval
  pvalmxlt <- res.mxlt$pval
  pvalmxrm <- res.mxrm$pval
  pvalmxrmt <- res.mxrmt$pval

  wrlin <- as.numeric(lin %*% wralla)
  stlin <- sqrt(mean((mhatall %*% lin)^2)) / sqrt(n)
  stlinp1 <- sqrt(mean((mhatallp1 %*% lin)^2)) / sqrt(n)
  stlinl <- sqrt(mean((mhatalll %*% lin)^2)) / sqrt(n)
  stlinrm <- sqrt(mean((mhatallrm %*% lin)^2)) / sqrt(n)
  stlinall <- c(stlin, stlinp1, stlinl, stlinrm)

  linci <- cbind(wrlin - 1.96 * stlinall, wrlin + 1.96 * stlinall)
  zvalin <- (wrlin - 1) / stlinall

  pvalin <- mean(abs(norr[,1]) > abs(zvalin[1]))
  pvalinp1 <- mean(abs(norr[,1]) > abs(zvalin[2]))
  pvalinn <- mean(abs(norr[,1]) > abs(zvalin[3]))
  pvalinrm <- mean(abs(norr[,1]) > abs(zvalin[4]))

  pvalinall <- c(pvalin, pvalinp1, pvalinn, pvalinrm)
  pvamxall <- c(pvalmxt, pvalmxp1t, pvalmxlt, pvalmxrmt)

  ci1 <- rbind(ci1l, ci1rm, ci1t, ci1p1t, ci1lt, ci1rmt)
  ci2 <- rbind(ci2l, ci2rm, ci2t, ci2p1t, ci2lt, ci2rmt)
  cil <- linci

  ## data-driven test
  rho <- sum(mhat01 * mhat02) / prod(sqrt(colSums(cbind(mhat01, mhat02)^2)))
  ard <- stall1[1]^2 + stall2[1]^2 - 2 * stall1[1] * stall2[1] * rho

  if (ard == 0) {
    ar <- 0.5
  } else {
    ar <- (stall2[1]^2 - stall1[1] * stall2[1] * rho) / ard
  }

  ar <- min(max(ar, 0), 1) #
  lin_ar <- c(ar, 1 - ar)
  wrlin_ar <- as.numeric(lin_ar %*% wralla[,1])
  stlin_ar <- sqrt(mean((mhatall %*% lin_ar)^2)) / sqrt(n)

  crnld <- 1 + 10 / max(c(n, 100))

  lincileswr <- c(
    wrlin_ar * exp(-1.96 * crnld * stlin_ar / wrlin_ar),
    wrlin_ar * exp(1.96 * crnld * stlin_ar / wrlin_ar)
  )

  zvalinleswr <- wrlin_ar * log(wrlin_ar) / (stlin_ar * crnld)
  peswrlog <- mean(abs(ak) > abs(zvalinleswr))

  results <- list()

  results$wralla <- wralla #nonterminal & terminal ratios
  results$ci1 <- ci1
  results$ci2 <- ci2
  results$cil <- cil
  results$mx <- c(mxt, mxp1t, mxlt, mxrmt)
  results$pvalmx <- c(pvalmxt, pvalmxp1t, pvalmxlt, pvalmxrmt)
  results$zvalin <- zvalin
  results$wrlin <- wrlin
  results$stlinall <- stlinall
  results$pvalinall <- pvalinall
  results$lin_ar <- lin_ar
  results$wrlin_ar <- wrlin_ar
  results$lincileswr <- lincileswr
  results$zvalinleswr <- zvalinleswr
  results$peswrlog <- peswrlog

  results


}


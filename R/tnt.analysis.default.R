tnt.analysis.default <- function(nonterm_time, nonterm_event, term_time, term_event, z,
                                 lin = c(0.5, 0.5), alpha = 0.05, repnum = 1E6, tau_r = 0.9, ...) {

  match.call()

  yh <- as.numeric(nonterm_time)
  hcen <- as.numeric(nonterm_event)
  yd <- as.numeric(term_time)
  dcen <- as.numeric(term_event)
  z <- as.numeric(z)
  lin_raw <- lin
  z_alpha <- qnorm(1 - alpha/2)
  n <- length(z)

  trp <- max(c(yh, yd)) + 0.001
  trpn <- trp * tau_r
  crn <- 1 + 15 / max(c(n, 100))

  ## data summary
  sum_table <- fn_ds(hcen, dcen, z)

  ## Estimates
  # each
  ginfh <- fn_ginf(yh, hcen, z, trp, trpn, z_alpha)
  ginfd <- fn_ginf(yd, dcen, z, trp, trpn, z_alpha)
  ci1 <- ginfh$ci
  ci2 <- ginfd$ci

  # all
  wralla <- rbind(h = ginfh$wrall, d = ginfd$wrall)
  stalla <- rbind(h = ginfh$stall, d = ginfd$stall)

  mhatall <- cbind(h = ginfh$mhat0, d = ginfd$mhat0)
  mhatallp1 <- cbind(h = ginfh$mhatp10, d = ginfd$mhatp10)
  mhatalll <- cbind(h = ginfh$mhatg0, d = ginfd$mhatg0)
  mhatallrm <- cbind(h = ginfh$mhatrm, d = ginfd$mhatrm)

  rho <- sum(mhatall[,1] * mhatall[,2]) / sqrt(sum(mhatall[,1]^2) * sum(mhatall[,2]^2))

  zvat <- c(ginfh$zvat, ginfd$zvat)
  zvap1t <- c(ginfh$zvap1t, ginfd$zvap1t)
  zvalt <- c(ginfh$zvalt, ginfd$zvalt)
  zvarmt <- c(ginfh$zvarmt, ginfd$zvarmt)

  mhat_list <- list(mhatall, mhatallp1, mhatalll, mhatallrm)
  z_list <- list(zvat, zvap1t, zvalt, zvarmt)

  norr <- matrix(rnorm(repnum * 2), ncol = 2)
  res_mx <- Map(fn_pval, mhat_list, z_list, MoreArgs = list(norr = norr, crn = crn))
  mx <- sapply(res_mx, `[[`, "mxo")
  pvalmx <- sapply(res_mx, `[[`, "pval")

  stlinall <- sapply(mhat_list, function(m) sqrt(mean((m %*% lin)^2)) / sqrt(n))

  wrlin <- as.numeric(lin %*% wralla)
  cil <- cbind(wrlin - z_alpha * stlinall, wrlin + z_alpha * stlinall)
  zvalin <- (wrlin - 1) / stlinall
  pvalinall <- 2 * pnorm(abs(zvalin), lower.tail = F)

  # data-driven test
  s1 <- stalla[1,1]; s2 <- stalla[2,1]
  ard <- s1^2 + s2^2 - 2 * s1 * s2 * rho
  ar <- if (ard == 0) 0.5 else (s2^2 - s1 * s2 * rho) / ard
  ar <- min(max(ar, 0), 1)  # clip between 0 and 1

  lin_ar <- c(ar, 1 - ar)
  wrlin_ar <- sum(lin_ar * wralla[,1])
  stlin_ar <- sqrt(mean((mhatall %*% lin_ar)^2)) / sqrt(n)

  crnld <- 1 + 10 / max(n, 100)
  log_term <- z_alpha * crnld * stlin_ar / wrlin_ar
  lincileswr <- wrlin_ar * exp(c(-log_term, log_term))
  zvalinleswr <- wrlin_ar * log(wrlin_ar) / (stlin_ar * crnld)
  peswrlog <- 2*pnorm(abs(zvalinleswr), lower.tail = F)

  ## Chi-square test
  vecr <- zvat/crn
  cmat <- matrix(c(1, rho, rho, 1), 2, 2)
  # chi <- as.numeric(t(vecr) %*% solve(cmat) %*% vecr)
  chi <- fn_solve(vecr, cmat)
  pvachi <- pchisq(chi, df = 2, lower.tail = FALSE)

  ## tnt.analysis for first obs
  y_first <- pmin(yh, yd)
  d_first <- dcen + (1 - dcen) * hcen
  logrank <- fn_logrank(y_first, d_first, z)
  wrtest0 <- fn_wr0(yh, hcen, yd, dcen, z)

  ## Test of Proportional Hazards test
  mhatph01 <- ginfh$mhatph
  mhatph02 <- ginfd$mhatph
  zvaph1 <- ginfh$zvaph
  zvaph2 <- ginfd$zvaph

  rhoph <- sum(mhatph01 * mhatph02)/prod(sqrt(colSums(cbind(mhatph01, mhatph02)^2)))
  zph12 <- c(zvaph1, zvaph2)
  cmath <- matrix(c(1, rhoph, rhoph, 1), 2, 2)

  svd_cmath <- svd(cmath)
  rsigh <- svd_cmath$u %*% diag(sqrt(svd_cmath$d)) %*% t(svd_cmath$v)
  sah <- norr %*% rsigh;
  mxph <- max(abs(c(zvaph1, zvaph2)));
  pvalph <- mean(pmax(abs(sah[,1]), abs(sah[,2])) > mxph)

  # Test of Equal Hazard Ratios
  lin_eq <- c(-1, 1)
  wrnon <- wralla[1,1]
  wrter <- wralla[2,1]

  qtn <- 1 - exp(-wrnon^2)
  qtt <- 1 - exp(-wrter^2)
  dwrl <- log(-log(qtn)) - log(-log(qtt))
  dtn <- 2 * wrnon * exp(-wrnon^2)/qtn/log(qtn)
  dtt <- 2 * wrter * exp(-wrter^2)/qtt/log(qtt)
  dtm <- diag(c(dtn, dtt))
  stl <- sqrt(mean((mhatall %*% dtm %*% lin_eq)^2))/sqrt(n)

  zvaephl <- dwrl/stl
  pvaephl <- 2*pnorm(abs(zvaephl), lower.tail = F)

  ci1 <- ci1[-c(1:2),]
  ci2 <- ci2[-c(1:2),]

  names(zvalin) <- names(pvalinall) <- c("ESWR", "LRGRE", "RICH", "RITCH")
  names(mx) <- names(pvalmx) <- c("ESWR", "LRGRE", "RICH", "RITCH")
  names(zvalinleswr) <- names(peswrlog) <- "ESWR"
  colnames(wralla) <- c("ESWR", "LRGRE", "RICH", "RITCH")
  rownames(ci1) <- rownames(ci2) <- rownames(cil) <- c("ESWR", "LRGRE", "RICH", "RITCH")
  colnames(ci1) <- colnames(ci2) <- colnames(cil) <- c("lower", "upper")
  names(wrlin) <- c("ESWR", "LRGRE", "RICH", "RITCH")
  names(wrlin_ar) <- "ESWR"
  names(lincileswr) <- c("lower", "upper")

  ## Output
  results <- list()

  results$alpha <- alpha
  results$lin <- lin
  results$rho <- rho
  results$lin_ar <- lin_ar
  results$sum_table <- sum_table
  results$global.tests$lincomb$stat <- zvalin
  results$global.tests$lincomb$p <- pvalinall
  results$global.tests$lincombw$stat <- zvalinleswr
  results$global.tests$lincombw$p <- peswrlog
  results$global.tests$max$stat <- mx
  results$global.tests$max$p <- pvalmx
  results$global.tests$chisq$stat <- chi
  results$global.tests$chisq$p <- pvachi
  results$logrank <- logrank
  results$wrtest0 <- wrtest0
  results$ci.results$indep$est <- wralla
  results$ci.results$indep$nonterm <- ci1
  results$ci.results$indep$term <- ci2
  results$ci.results$lincomb$est <- wrlin
  results$ci.results$lincomb$ci <- cil
  results$ci.results$lincombw$est <- wrlin_ar
  results$ci.results$lincombw$ci <- lincileswr
  results$phtest$est <- mxph
  results$phtest$p <- pvalph
  results$eqtest$est <- zvaephl
  results$eqtest$p <- pvaephl

  class(results) <- "tnt.analysis"

  results

}


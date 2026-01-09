print.tnt.analysis <- function(x, ...) {

  fmt  <- function(v) sprintf("%.3f", as.numeric(v))
  fmtp <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("= %.4f", as.numeric(p)))
  len <- 81
  ci_level <- round((1 - x$alpha) * 100, 2)

  zvalinf     <- fmt(x$global.tests$lincomb$stat)
  pvalinallf  <- fmtp(x$global.tests$lincomb$p)
  zvalinleswrf <- fmt(x$global.tests$lincombw$stat)
  peswrlogf   <- fmtp(x$global.tests$lincombw$p)
  mxf         <- fmt(x$global.tests$max$stat)
  pvalmxf     <- fmtp(x$global.tests$max$p)
  chif        <- fmt(x$global.tests$chisq$stat)
  pvachif     <- fmtp(x$global.tests$chisq$p)

  logrank_statf <- fmt(x$logrank$stat)
  logrank_pvalf <- fmtp(x$logrank$pval)

  wr0f        <- fmt(x$wrtest0$wr)
  zv_wr0f     <- fmt(x$wrtest0$zv_wr)
  pval_wr0f   <- fmtp(x$wrtest0$pval_wr)

  wralla <- x$ci.results$indep$est
  ci1 <- x$ci.results$indep$nonterm
  ci2 <- x$ci.results$indep$term
  wrlin <- x$ci.results$lincomb$est
  cil <- x$ci.results$lincomb$ci
  wrlin_ar <- x$ci.results$lincombw$est
  lincileswr <- x$ci.results$lincombw$ci

  mxphf       <- fmt(x$phtest$est)
  zvaephlf    <- fmt(x$phtest$p)

  pvalphf     <- fmtp(x$eqtest$est)
  pvaephlf    <- fmtp(x$eqtest$p)

  ## OUTPUT
  sum_table <- x$sum_table
  colnames(sum_table) <- c("Treatment", "Control", "All")
  rownames(sum_table) <- c("   Non-terminal event only", "   Terminal event only", "   Non-terminal + Terminal")
  # rownames(sum_table) <- c("   H only", "   H & D", "   D only")

  cat(strrep("=", len), "\n")
  cat("I. -------------- < Summary of Observed Events > -------------------------------\n")
  print(sum_table)
  cat(strrep("-", len), "\n\n")

  cat("II. ------------- < Tests of the Global Null Hypothesis > -----------------------\n\n")
  cat(sprintf(">> Linear combination test [l = (%.1f, %.1f)]\n", x$lin[1], x$lin[2]))
  for (i in 1:4) {
    cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n",
                c("ESWR", "LRGRE", "RICH", "RITCH")[i],
                zvalinf[i], pvalinallf[i]))
  }
  cat(sprintf("\n>> Linear combination test with data-driven weights [l = (%.2f, %.2f)]\n", x$lin_ar[1], x$lin_ar[2]))
  cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n", "ESWR", zvalinleswrf, peswrlogf))
  cat("\n>> Maximum test\n")
  for (i in 1:4) {
    cat(sprintf("  [%-5s]  statistic = %s,  p-value %s\n",
                c("ESWR", "LRGRE", "RICH", "RITCH")[i],
                mxf[i], pvalmxf[i]))
  }
  cat(paste("\n>> Chi-squred test: statistic = ", chif, ", p-value ", pvachif, sep = ""), "\n")
  cat(paste(">> Win ratio test: Chisq.statistic = ", zv_wr0f, ", p-value ", pval_wr0f, sep = ""), "\n")
  cat(paste(">> Log-rank test using time to first event: statistic = ", logrank_statf, ", p-value ", logrank_pvalf, sep = ""), "\n")
  cat(strrep("-", len), "\n\n")

  cat(paste("III. ----------- < Confidence Intervals (", ci_level, "% CI) > ------------------------------", sep = ""), "\n")
  cat(">> Non-terminal event\n")
  for (i in 1:4) {
    cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n",
                c("ESWR", "LRGRE", "RICH", "RITCH")[i],
                wralla[1, i], ci1[i, 1], ci1[i, 2]))
  }
  cat("\n>> Terminal event\n")
  for (i in 1:4) {
    cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n",
                c("ESWR", "LRGRE", "RICH", "RITCH")[i],
                wralla[2, i], ci2[i, 1], ci2[i, 2]))
  }
  cat(sprintf("\n>> Linear combinations [l = (%.1f, %.1f)]\n", x$lin[1], x$lin[2]))
  for (i in 1:4) {
    cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n",
                c("ESWR", "LRGRE", "RICH", "RITCH")[i],
                wrlin[i], cil[i, 1], cil[i, 2]))
  }
  cat(sprintf("\n>> Data-driven combinations [l = (%.2f, %.2f)]\n", x$lin_ar[1], x$lin_ar[2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ESWR",
              wrlin_ar, lincileswr[1], lincileswr[2]))
  cat(strrep("-", len), "\n\n")

  cat("IV. ------------ < Test of Proportional Hazards > -------------------------------\n")
  cat(paste(">> statistic = ", mxphf, ", p-value ", pvalphf, sep = ""), "\n")
  cat(strrep("-", len), "\n\n")

  cat("V. ------------- < Test of Equal Hazard Ratios > --------------------------------\n")
  cat(paste(">> statistic = ", zvaephlf, ", p-value ", pvaephlf, sep = ""), "\n")
  cat(strrep("-", len), "\n")
  cat(strrep("=", len))
}

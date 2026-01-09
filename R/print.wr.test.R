print.tnt.analysis <- function(x, ...) {

  fmt  <- function(v) sprintf("%.3f", as.numeric(v))
  fmtp <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("= %.4f", as.numeric(p)))
  len <- 81
  ci_level <- round((1 - x$alpha) * 100, 2)

  wrlinf      <- fmt(x$wrlin)
  mxf         <- fmt(x$mx)
  chif        <- fmt(x$chi)
  wr0f        <- fmt(x$wrtest0$wr)
  zv_wr0f     <- fmt(x$wrtest0$zv_wr)
  wrlin_arf   <- fmt(x$wrlin_ar)
  logrank_statf <- fmt(x$logrank$stat)
  mxphf       <- fmt(x$mxph)
  zvaephlf    <- fmt(x$zvaephl)
  zvalinf     <- fmt(x$zvalin)
  zvalinleswrf <- fmt(x$zvalinleswr)
  pvalinallf  <- fmtp(x$pvalinall)
  peswrlogf   <- fmtp(x$peswrlog)
  pvalmxf     <- fmtp(x$pvalmx)
  pvachif     <- fmtp(x$pvachi)
  pvalphf     <- fmtp(x$pvalph)
  pvaephlf    <- fmtp(x$pvaephl)
  pval_wr0f   <- fmtp(x$wrtest0$pval_wr)
  logrank_pvalf <- fmtp(x$logrank$pval)

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
                x$wralla[1, i], x$ci1[i + 2, 1], x$ci1[i + 2, 2]))
  }
  cat("\n>> Terminal event\n")
  for (i in 1:4) {
    cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n",
                c("ESWR", "LRGRE", "RICH", "RITCH")[i],
                x$wralla[2, i], x$ci2[i + 2, 1], x$ci2[i + 2, 2]))
  }
  cat(sprintf("\n>> Linear combinations [l = (%.1f, %.1f)]\n", x$lin[1], x$lin[2]))
  for (i in 1:4) {
    cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n",
                c("ESWR", "LRGRE", "RICH", "RITCH")[i],
                x$wrlin[i], x$cil[i, 1], x$cil[i, 2]))
  }
  cat(sprintf("\n>> Data-driven combinations [l = (%.2f, %.2f)]\n", x$lin_ar[1], x$lin_ar[2]))
  cat(sprintf("  [%-5s]  estimate = %.3f,  CI = (%.3f, %.3f)\n", "ESWR",
              x$wrlin_ar, x$lincileswr[1], x$lincileswr[2]))
  cat(strrep("-", len), "\n\n")

  cat("IV. ------------ < Test of Proportional Hazards > -------------------------------\n")
  cat(paste(">> statistic = ", mxphf, ", p-value ", pvalphf, sep = ""), "\n")
  cat(strrep("-", len), "\n\n")

  cat("V. ------------- < Test of Equal Hazard Ratios > --------------------------------\n")
  cat(paste(">> statistic = ", zvaephlf, ", p-value ", pvaephlf, sep = ""), "\n")
  cat(strrep("-", len), "\n")
  cat(strrep("=", len))
}

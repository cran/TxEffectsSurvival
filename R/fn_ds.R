fn_ds <- function(hcen, dcen, z) {

  h_only   <- (hcen == 1 & dcen == 0)
  h_and_d  <- (hcen == 1 & dcen == 1)
  d_only   <- (hcen == 0 & dcen == 1)

  count_group <- function(group) {
    z1  <- sum(group & z == 1)
    z0  <- sum(group & z == 0)
    all <- sum(group)
    c(`z=1` = z1, `z=0` = z0, all = all)
  }

  out <- rbind(
    h_only  = count_group(h_only),
    d_only  = count_group(d_only),
    h_and_d = count_group(h_and_d)
  )

  return(out)
}


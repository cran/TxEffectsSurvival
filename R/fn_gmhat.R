fn_gmhat <- function(d1, d2, zd1, zd2, ff1, ff2, yk1, yk2, n) {
  
  mhat <- rep(0, n);
  mhat[zd1] <- d1 * ff1 - cumsum(ff1 * d1 /yk1)
  mhat[zd2] <- d2 * ff2 - cumsum(ff2 * d2 /yk2)
  
  mhat
  
}
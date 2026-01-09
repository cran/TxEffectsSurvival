fn_glinsim <- function(lin, wralla, mhatall, mhatp1all, mhatp2all, ak, ctv, n, z_alpha) {
  
  wrlin <- as.numeric(lin %*% wralla);
  stlin <- sqrt(colMeans((mhatall %*% lin)^2))/sqrt(n)
  stlinp1 <- sqrt(colMeans((mhatp1all %*% lin)^2))/sqrt(n)
  stlinp2 <- sqrt(colMeans((mhatp2all %*% lin)^2))/sqrt(n)
  stlinall <- c(stlin,  stlinp1, stlinp2);
  
  lincil <- cbind(wrlin * exp(-z_alpha * stlinall/wrlin), wrlin * exp(z_alpha * stlinall/wrlin))
  zvalinl <- wrlin * log(wrlin)/stlinall;
  pvalinl <- mean(abs(ak) > abs(zvalinl[1]));
  
  linci <- cbind(wrlin - z_alpha * stlinall, wrlin + z_alpha * stlinall)
  zvalin <- (wrlin - ctv)/stlinall;
  pvalin <- mean(abs(ak) > abs(zvalin[1]));
  
  res <- list()
  res$stlin <- stlin
  res$lincil <- lincil
  res$zvalinl <- zvalinl
  res$pvalinl <- pvalinl
  res$linci <- linci
  res$zvalin <- zvalin
  res$pvalin <- pvalin
  res$wrlin <- wrlin
  
  res
  
}
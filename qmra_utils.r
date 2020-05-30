
qmra.exponential <- function(k,dose) {
  return(1 - exp(-k*dose))
}
qmra.bp <- function(alpha,N50,dose) {
  res <- 1 - ( 1 + ( dose / N50 ) * (2 ^ ( 1 / alpha) -1 ) ) ^ (-alpha)
  return(res)
}
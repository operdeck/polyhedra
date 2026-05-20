# shift elements from f up or down
# NB does not work as expected for n > |x|
shiftrotate <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

## TODO = seq_len
# safe sequence from 1 to n - returning empty vector if it would be backwards
safeseq <- function(n, by=1)
{
  if (n<1 & by>0) return(c())
  if (n>1 & by<0) return(c())
  return(seq(from=1, to=n, by=by))
}

# often used in polyhedra
phi <- (1+sqrt(5))/2


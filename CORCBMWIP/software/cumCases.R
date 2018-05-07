cumCases <- function(z, scores){
  zTable <- rev(table(cut(z, c(scores, Inf), right = FALSE)))
  return(cumsum(as.vector(zTable)))
}
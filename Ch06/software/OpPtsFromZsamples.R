# OpPtsFromZsamples.R
#z1 is the K1 array of ratings for non-diseased cases
#z2 is the K2 array of ratings for diseased cases
OpPtsFromZsamples <- function( z1, z2 ) {

  thresholds <- unique(c(z1, z2))
  z1Table <- rev(table(cut(z1, c(thresholds, Inf), right = FALSE)))
  z2Table <- rev(table(cut(z2, c(thresholds, Inf), right = FALSE)))
  
  FPF <- cumsum(as.vector(z1Table)) / length(z1)
  TPF <- cumsum(as.vector(z2Table)) / length(z2)
  
  FPF <- c(0, FPF)
  TPF <- c(0, TPF)
  
  return( list(
    FPF = FPF,
    TPF = TPF
  ) )
   
}

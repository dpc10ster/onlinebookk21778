desired_cdf <- function( r, param ) {
  f <- param$f
  r_1 <- param$r_1
  r_2 <- param$r_2
  
  zp <- array(dim = length(r))
  for (i in 1:length(r)){
    if( r[i] < r_1 ) {
      zp[i] <- max((f * r[i] / r_1)*100, 0)
    } else if( r[i] < r_2 ) {
      zp[i] <- f*100
    } else {
      zp[i] <- ( 1 - f ) * (r[i]-r_2)*100
    }
  }
  
  return( zp )
}

desired_cdf_inv <- function( zp, param ) {
  f <- param$f
  r_1 <- param$r_1
  r_2 <- param$r_2
  
  if( zp < f ) {
    r <- zp * r_1 / f
  } else {
    r <- r_2 + ( zp - f ) * ( 100 - r_2 ) / ( 1 - f )
  }
  
  return( r )
}

desired_cdf <- function(z,a1,a2,s1,s2) {
	y <- (1-f)*pnorm((z-a1)/s1)*100+f*pnorm((z-a2)/s2)*100

	return( y )
}



desired_cdf_inv_diff <- function( r, param ) {
	f <- param$f
	r_1 <- param$r_1
	r_2 <- param$r_2
	sigma <- param$sigma
	zp <- param$zp

	return( abs( zp - desired_cdf( r, param ) ) )
}



desired_cdf_inv <- function( zp, param ) {
	inv_param <- list(
		f = param$f,
		r_1 = param$r_1,
		r_2 = param$r_2,
		sigma = param$sigma,
		zp = zp
	)

	r <- 50		
	re <- optim( r, desired_cdf_inv_diff, param = inv_param, 
		method = "Brent", lower = 0, upper = 100 )

	ret <- re$par

	return( ret )
}

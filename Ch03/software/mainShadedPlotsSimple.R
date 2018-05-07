#mainShadedPlotsSimple.R
rm(list = ls())
mu <- 3;sigma <- 1;zeta <- 1;options(digits=3)
step <- 0.1

lowerLimit<- -1 # lower limit
upperLimit <- mu + 3*sigma # upper limit

seqNor <- seq(zeta,upperLimit,step)
cord.x <- c(zeta, seqNor,upperLimit)
# need two y-coords at each end point of range;
# one at zero and one at value of function
cord.y <- c(0,dnorm(seqNor),0)
curve(dnorm(x,0,1),xlim=c(lowerLimit,upperLimit),col='blue',
      ylab = "pdfs", xlab ="z")
polygon(cord.x,cord.y,col='blue')

curve(dnorm(x,mu,sigma),xlim=c(lowerLimit,upperLimit), add = TRUE, col = 'red')

crossing <- uniroot(function(x) dnorm(x) - dnorm(x,mu,sigma),
                    lower = 0, upper = 3)$root

crossing <- max(c(zeta, crossing))
seqAbn <- seq(crossing,upperLimit,step)
cord.x <- c(seqAbn, rev(seqAbn))
# reason for reverse
# we want to explicitly define the polygon
# we dont want R to close it

cord.y <- c()
for (i in seq(1,length(cord.x)/2)) {
  cord.y <- c(cord.y,dnorm(cord.x[i],mu, sigma))
}
for (i in seq(1,length(cord.x)/2)) {
  cord.y <- c(cord.y,dnorm(cord.x[length(cord.x)/2+i]))
}
polygon(cord.x,cord.y,lty = 0, col='red')

seqAbn <- seq(zeta,upperLimit,step)
seqAbn <- rep(seqAbn,each = 2)
for (i in seq(1,length(seqAbn), 2)) {
  # define xs and ys of two points, separated only along y-axis
  x <- c(seqAbn[i], seqAbn[i+1])
  y <- c(0,dnorm(seqAbn[i+1], mu, sigma))
  # draw vertical line between them
  lines(x,y, col = 'red', lty = 1, lwd = 2)
}

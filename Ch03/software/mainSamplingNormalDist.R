#MainSamplingNormalDist.R
rm(list = ls()); set.seed(1); options(digits=3) 
# get 2 sets of 10 random samples from the unit normal distribution
x <- rnorm(10);cat("1st set of 10 random samples are \n",x,"\n\n")
x <- rnorm(10);cat("2nd set of 10 random samples are \n",x,"\n\n")
# mean and standard deviation of 10,000 new samples from the unit normal distribution
cat("mean of 10,000 new samples =", mean(x1 <- rnorm(10000)),"\n")
cat("std. of above 10,000 samples =", sqrt(var(x1)),"\n")

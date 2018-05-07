rm( list = ls() ) #MainHelloWorld.R 
seed <- 1;set.seed(seed)
cat("Hello world!\n")
N <- 10000
samples <- runif(N)
hist(samples)

MeansOfSamples <- array(N)
for (i in 1:N) {
  samples <- runif(N)
  MeansOfSamples[i] <- mean(samples)
}

hist(MeansOfSamples, breaks = 100)

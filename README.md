# OR-2-HW

#Brownian motion is very easy to simulate. To start off, let's simulate a single
#instance of Brownian motion for 100 generations of discrete time in which the
#variance of the diffusion process is ??2 = 0.01 per generation.

t <- 0:100  # time
sig2 <- 0.01
## first, simulate a set of random deviates
x <- rnorm(n = length(t) - 1, sd = sqrt(sig2))
## now compute their cumulative sum
x <- c(0, cumsum(x))
plot(t, x, type = "l", ylim = c(-2, 2))

#The trick here is that we draw random normal deviates for the change over each
#t time intervals; and then to get the state of our chain at each time interval
#we just compute the cumulative sum of all the individual changes using cumsum.
#We can do this because the distribution of the changes under Brownian motion is
#invariant and does not depend on the state of the chain.

#We can also easily do a whole bunch of simulations like this at once, using the
#same conditions:
  
nsim <- 100
X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 
              1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-2, 2), type = "l")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)


# To see how the outcome depends on ??2, let's compate the result when we divide
# sig2 by 10:

X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2/10)), nsim, length(t) - 
1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-2, 2), type = "l")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)


# There are a number of different ways we could've done this. For example,
# instead of using apply, which economizes on code, we could have used a for
# loop as follows:

X <- matrix(0, nsim, length(t))
for (i in 1:nsim) X[i, ] <- c(0, cumsum(rnorm(n = length(t) - 1, sd = sqrt(sig2))))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-2, 2), type = "l")
for (i in 1:nsim) lines(t, X[i, ])


# As mentioned above, the expected variance under Brownian motion is just ??2 ×
# t. To see this easiest, we can just do the following. Here I will use 10,000
# simulations for 100 generations under the same conditions to "smooth out" our
# result:
  
nsim <- 10000
X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 
              1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
v <- apply(X, 2, var)
plot(t, v, type = "l", xlab = "time", ylab = "variance among simulations")

# The result will be Gaussian - due to the Central Limit Theorem. Here is an
# illustration of that:
  
t <- 0:100  # time
sig2 <- 0.01
nsim <- 1000
## we'll simulate the steps from a uniform distribution with limits set to
## have the same variance (0.01) as before
X <- matrix(runif(n = nsim * (length(t) - 1), min = -sqrt(3 * sig2), max = sqrt(3 * 
                                                                                  sig2)), nsim, length(t) - 1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-2, 2), type = "l")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)


var(X[, length(t)])


hist(X[, length(t)])

plot(density(X[, length(t)]))


# Brownian motion is generally assumed to be trendless; however it is possible
# to simulate and (under some conditions) fit a model of Brownian evolution with
# a trend. Here is an example of a simulation (using the same general approach
# as above) with a trend.

nsim = 100
X <- matrix(rnorm(mean = 0.02, n = nsim * (length(t) - 1), sd = sqrt(sig2/4)), 
            nsim, length(t) - 1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-1, 3), type = "l")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)

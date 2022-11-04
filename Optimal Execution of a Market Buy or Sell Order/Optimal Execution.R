# liquidity risk
# Optimal execution problem

# Problem:
# minimize: lambda_param * sqrt( sum(sigma^2 * x_t^2) ) + sum( 0.5 * abs(q_t) * p(q_t) )
# such that: sum(q_t) = S

### where:
### we sum over t=1 to t=n, where
### n            = number of time increments we have to execute all of a trade
### lambda_param = risk aversion to not executing immidiately
### sigma        = volatility of the asset you want to buy
### x_t          = quantity we hold of the asset at time t
### q_t          = quantity we trade at each time t
### p(q_t)       = bid-ask spread paid for an order of q_t
### S            = total quantity we want to sell

### We set S = 100, n = 5, sigma = 0.1


# objective function
tobeminimized <- function(q, lambda){
  # q == # of shares to buy each day
  # We substituted in the constraint (ie "such that: sum(q_t) = S") to have an n-1 dimensional problem: 
  ### for ex, once you know how much to trade on the first n-1 days, the # of trades on day n is also known
  q <- c(q, S-sum(q)) # substituted in the constraint here to have an n-1 dimensional problem
  
  # x == # of remaining shares you still have to execute on each day
  ### for ex, if S=100 and q=[15, 25, 18, 22, 20], we have x = [100, 85, 60, 42, 20]
  x <- c(S, S - cumsum(q[1:(n-1)]))
  
  # p = bid-ask spread paid for trading q_t shares at once
  p <- function(q_t){0.1 + 0.05*exp(0.03*q_t)}
  
  
  
  lambda*sqrt(sum(sigma2*x^2)) + 0.5*sum(abs(q)*p(q))
}


# parameters
S <- 100
n <- 5
sigma2 <- 0.1^2

# initial guess
# starting point is 15 1st day, 25 2nd day, ..., 22 4th day, and (100 - sum(q0)) 5th day
q0 <- c(15, 25, 18, 22)

# try with a specific lambda
q_opt <- optim(q0, tobeminimized, lambda = 1)$par
q_opt <- c(q_opt, 100-sum(q_opt))
q_opt


# try with many lambdas
lambdas <- c(0, 0.1, 1, 10, 100, 1000)

q <- sapply(lambdas, function(z)optim(q0, tobeminimized, lambda = z)$par)
q <- rbind(q, S-colSums(q))
q <- cbind(1:n, q)

q <- data.frame(q)
colnames(q) <- c('t',lambdas)

# plot the solutions
library(tidyverse)
ggplot(gather(q, lambda, 'q', 2:length(q)), aes(x = t, y = q)) + geom_line(aes(color = lambda)) +
  theme_bw() + ylab('# of shares to trade') + xlab('Day i') +
  ggtitle('Optimal execution varies with lambda (aversion to waiting)') +
  theme(plot.title = element_text(size=20, hjust = 0.5)) +
  theme(text=element_text(size=12)) +
  theme(axis.text=element_text(size=12))

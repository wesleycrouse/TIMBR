#load example data
data("mcv.data", package="TIMBR")

y <- mcv.data$y
prior.D <- mcv.data$prior.D
prior.D$fixed.diplo <- NULL

#priors over configurations in order of increasing runtime
#prior.M <- list(M.IDs="0,1,2,3,4,5,6,7", ln.probs=log(1))
prior.M <- list(M.IDs=c("0,1,0,1,1,0,0,0", "0,1,0,1,1,2,0,0", "0,1,2,3,4,5,6,7"), ln.probs=rep(log(0.5),2))
#prior.M <- mcv.data$prior.M$list
#prior.M <- ewenss.calc(J, list(type="gamma", shape=1, rate=2.333415))

rm(mcv.data)

####################
#define the stan_TIMBR function

stan_TIMBR <- function(y, prior.D, prior.M, p=0.01, prior.phi.v=2, chains=4, Z = NULL, W = NULL){
  #specify data for stan
  data.stan <- list()
  data.stan$y <- y
  data.stan$P <- prior.D$P
  data.stan$A <- prior.D$A
  data.stan$S <- ncol(data.stan$P)
  data.stan$J <- ncol(data.stan$A)
  data.stan$L <- array(numeric(),c(length(prior.M$M.IDs),data.stan$J,data.stan$J)) 
  for (i in 1:length(prior.M$M.IDs)){
    M <- TIMBR:::M.matrix.from.ID(prior.M$M.IDs[i])
    data.stan$L[i,,] <- t(chol((1-p)*M%*%t(M) + p*diag(data.stan$J)))
  }
  data.stan$M <- dim(data.stan$L)[1]
  data.stan$ln_prior_M = as.array(prior.M$ln.probs)
  data.stan$N <- length(y)
  data.stan$kappa <- 0.001*2
  data.stan$lambda <- 0.001*2
  data.stan$tau_mu <- sqrt(1000)
  data.stan$tau_delta <- sqrt(1000)
  data.stan$nu <- prior.phi.v
  data.stan$ones <- rep(1,data.stan$J)
  data.stan$zeros <- rep(0,data.stan$N)
  if (is.null(W)){
    data.stan$w_inv <- rep(1,data.stan$N)^(-1)
  } else {
    data.stan$w_inv <- W^(-1)
  }
  if (is.null(Z)){
    data.stan$Z <- matrix(0,data.stan$N,0)
  } else {
    data.stan$Z <- Z
  }
  data.stan$K <- ncol(data.stan$Z)
  
  #compile or load the stan model
  #stan.model <- rstan::stan_model(file="stan_TIMBR.stan")
  stan.model <- readRDS("stan_TIMBR.rds")
  
  if (data.stan$M==1){
    fit <- rstan::sampling(stan.model, data=data.stan, chains = chains, include=F, pars=c("ln_post_M"))
  } else {
    fit <- rstan::sampling(stan.model, data=data.stan, chains = chains)
  }
  
  fit
}

library(rstan)
options(mc.cores = parallel::detectCores())

fit <- stan_TIMBR(y, prior.D, prior.M)

posterior <- as.matrix(fit)

####################
#plot the posterior haplotype effects
library(ggplot2)
library(bayesplot)

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior[,grep("beta_h", colnames(posterior))],
           prob = 0.8) + plot_title

####################
#compute posterior configuration probabilities
post.M <- exp(posterior[,grep("ln_post_M", colnames(posterior)), drop = FALSE])
post.M <- colMeans(post.M)

results <- -sort(-post.M)
names(results) <- prior.M$M.IDs[order(-post.M)]

head(results)

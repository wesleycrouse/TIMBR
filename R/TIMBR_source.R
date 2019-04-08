#' @keywords internal
sumtozero.contrast <- function(K){
  u <- 1/((K-1)^(0.5))
  v <- (-1+(K^(0.5)))*((K-1)^(-1.5))
  w <- (K-2)*v+u
  
  C <- matrix(-v,K-1,K-1)
  diag(C) <- rep(w,K-1)
  rbind(C, rep(-u,K-1))
}

#' @keywords internal
m.rename <- function(m, string=T){
  unique.m <- unique(m)
  output <- sapply(m, function(y){match(x=y, unique.m)-1})
  ifelse(string, paste(output, collapse=","), output)
}

#' @keywords internal
ln.bell <- function(J){
  #exact calculation
  ln.b <- log(numbers::bell(J))
  
  if (ln.b==Inf){
    #approximation to avoid overflow
    lambertW.J <- VGAM::lambertW(J)
    ln.b <- -0.5*log(J) + (J+1/2)*(log(J)-log(lambertW.J)) + (J/lambertW.J - J - 1)
  }
  
  ln.b
}

#' @keywords internal
m.from.M.ID <- function(M.ID){
  as.numeric(unlist(strsplit(M.ID, ","))) + 1
}

#' @keywords internal
M.matrix.from.ID <- function(M.ID){
  m <- m.from.M.ID(M.ID)
  J <- length(m)
  M <- matrix(0, J, max(m))
  M[cbind(1:J, m)] <- 1
  M
}

#' @keywords internal
ln.beta.prior.marginalized <- function(beta, sigma.sq, prior.phi.b, prior.phi.a=0.5){
  #Abramowitz and Stegun 55 p.505 "Confluent Hypergeometric Functions"
  d <- length(beta)
  
  #parameters for hypergeometric U
  z.U <- 0.5*sum(beta^2)/sigma.sq
  a.U <- prior.phi.b + 0.5*d
  b.U <- prior.phi.a + 0.5*d
  
  (-0.5*d)*log(2) - 0.5*d*log(pi) - 0.5*d*log(sigma.sq) - lbeta(0.5, prior.phi.b) + lgamma(a.U) + log(gsl::hyperg_U(a.U, b.U, z.U))
}

#' Tree-based Inference of Multiallelism via Bayesian Regression
#'
#' Posterior samples and Bayes Factors using the TIMBR model
#'
#' @param y vector of phenotype values for each strain
#' @param prior.D list of inputs for the prior distribution of strain diplotype states; see data(mcv.data) for an example
#' @param prior.M list of inputs for the prior distribution of the allelic series model; see data(mcv.data) for examples
#' @param prior.phi.b shape parameter for the beta prime prior distribution on the variance component
#' @param samples number of samples to draw from the full posterior
#' @param samples.ml number of samples to draw from the condtiional posterior (if necessary)
#' @param Z design matrix for intercept and covariates; first column must be a vector of ones, which is the default
#' @param W vector of replicates for each strain; one replicate per strain by default
#' @param verbose optionally report function progress
#'
#' @return a list of input parameters, posterior samples and marginal densities, and the marginal likelihood
#' 
#' @examples
#' #example data
#' data(mcv.data)
#' str(mcv.data)
#' 
#' #call TIMBR using CRP
#' results <- TIMBR(mcv.data$y, mcv.data$prior.D, mcv.data$prior.M$crp)
#' 
#' #report the Bayes Factor
#' results$ln.BF
#' 
#' #report posterior probabilities for the top allelic series models
#' head(results$p.M.given.y)
#' 
#' #report mean posterior haplotype effects
#' colMeans(results$post.hap.effects)
#'
#' @export
TIMBR <- function(y, prior.D, prior.M, prior.phi.b=1, samples=10000, samples.ml=10000, Z=NULL, W=NULL, verbose=T){
  
  TIMBR.sampler <- function(iterations, calc.null.ml=T, update.M=T, update.alpha=T){
    
    nglm.hyperparameters.ml <- function(MC, calc.partial.ml=T){
      #compute hyperparameters for the normal-gamma linear model
      if (update.M){
        d <- ncol(MC)
        n.params <- p + d
      }
      
      if (update.M | !fixed.diplo){
        X <- cbind(Z, DA%*%MC)
        XtWy <- crossprod(X,Wy)
        ZtWDAMC <- ZtWDA%*%MC
        CtMtAtDtWDAMC <- crossprod(MC, AtDtWDA)%*%MC
        V.star.inv <- rbind(cbind(ZtWZ, ZtWDAMC),cbind(t(ZtWDAMC),CtMtAtDtWDAMC))
      }
      
      if (d!=0){
        V.star.inv[cbind((p+1):n.params, (p+1):n.params)] <- CtMtAtDtWDAMC[cbind(1:d, 1:d)] + rep(phi.sq^(-1),d)
      }
      
      if (n.params > 0){
        L.Vt.inv <- backsolve(chol(V.star.inv), diag(n.params))
        V.star <- tcrossprod(L.Vt.inv)
      } else {
        L.Vt.inv <- V.star.inv
        V.star <- V.star.inv
      }
      
      m.star <- V.star%*%XtWy
      psi.star <- ytWy - c(crossprod(m.star, V.star.inv)%*%m.star)
      
      return.list <- list("kappa.star"=n, "psi.star"=psi.star, "m.star"=m.star, "V.star"=V.star, "L.Vt.inv"=L.Vt.inv)
      
      if (calc.partial.ml){
        partial.ln.ml <- -0.5*kappa.star*log(psi.star) - 0.5*d*log(phi.sq) + 0.5*determinant(V.star)$modulus[1]
        return.list$partial.ln.ml <- partial.ln.ml
      }
      
      return.list
    }
    
    sample.crp.concentration <- function(){
      #sample latent variable zeta
      ln.zeta <- log(rbeta(1, alpha+1, J))
      
      #sample alpha conditional on zeta
      z <- (prior.alpha.shape+K-1)/(J*(prior.alpha.rate-ln.zeta))
      pi <- z/(1+z)
      
      if (rbinom(1,1,pi)==1){
        alpha <- rgamma(1, prior.alpha.shape+K, prior.alpha.rate-ln.zeta)
      } else {
        alpha <- rgamma(1, prior.alpha.shape+K-1, prior.alpha.rate-ln.zeta)
      }
      
      alpha
    }
    
    #precompute matrix products if model and diplotypes are fixed
    if (!update.M & fixed.diplo){
      DAMC <- DA%*%MC
      X <- cbind(Z, DAMC)
      XtWy <- crossprod(X,Wy)
      ZtWDAMC <- ZtWDA%*%MC
      CtMtAtDtWDAMC <- crossprod(MC, AtDtWDA)%*%MC
      V.star.inv <- rbind(cbind(ZtWZ, ZtWDAMC),cbind(t(ZtWDAMC),CtMtAtDtWDAMC))
    }
    
    #create objects to store results
    post.M <- matrix(NA, iterations, J)
    post.MCbeta <- matrix(NA, iterations, J)
    post.delta <- matrix(NA, iterations, p)
    post.phi.sq <- rep(NA, iterations)
    post.sigma.sq <- rep(NA, iterations)
    post.alpha <- rep(NA, iterations)
    post.hyperparameters <- vector("list", iterations)
    post.K <- rep(NA, iterations)
    post.y.hat <- matrix(NA, iterations, n)
    p.D.given.y <- matrix(0, n, ncol.P)
    
    #iterate sampler
    for (i in 1:iterations){
      if (i%%1000==0 & verbose){
        print(i)
      }
      
      #compute matrix quantitites that depend on D
      if (!fixed.diplo){
        DA <- A[D.list,,drop=F]
        ZtWDA <- ZtW%*%DA
        AtDtWDA <- crossprod(DA*sqrt.W)
        AtDtWy <- crossprod(DA, Wy)
      }
      
      #sample allelic series matrix M
      if (update.M){
        if (model.type=="crp"){
          j.order <- 1:J
        } else {
          #randomize order due to non-exchangeable prior
          j.order <- sample(1:J)
        }
        
        #sample each row of M conditional on the other rows
        for (j in j.order){
          
          if (j!=j.order[1]){
            M.current <- list(M=M, M.list=M.list, M.posteriors=M.posteriors[[M.indicator]], new.index=M.list[j])
          } else {
            M.current <- list()
          }
    
          #set current row to zero and update matrix columns if necessary
          M[j,] <- 0
          
          if (!(1 %in% M[,M.list[j]])){
            M.update.index <- M.list[-j] > M.list[j]
            M.list[-j][M.update.index] <- M.list[-j][M.update.index] - 1
            M <- M[,-M.list[j],drop=F]
            M.current$new.index <- ncol(M)+1
          }
          
          M.list[j] <- NA
          K <- ncol(M)
          C <- contrast.list[[K]]
          
          #enumerate all possible assignments of current row of M
          M.list.space <- lapply(1:(K+1), function(x){M.list[j] <- x; M.list})
          MC.space <- lapply(1:K, function(x){C[M.list.space[[x]],,drop=F]})
          MC.space[[K+1]] <- contrast.list[[K+1]][M.list.space[[K+1]],,drop=F]
          
          #calculate prior for all possible assignments of current row of M
          if (model.type=="crp"){
            #analytic form for exchangeable prior
            colsums.M <- matrixStats::colSums2(M)
            M.ln.prior <- log(c(colsums.M, alpha))
          } else if (model.type=="uniform"){
            #constant non-exchangeable prior
            M.ln.prior <- rep(0, K+1)
          } else if (model.type=="list" | model.type=="mixture"){
            #arbitrary non-exchangeable prior or a mixture of that prior with the CRP
            #requires looking up priors via hash table using a unique naming scheme for M
            M.space.vec <- lapply(1:(K+1), function(x){M.list[j] <- x; M.list})
            
            #compute unique names for possible values of M
            if (hash.names){
              #use a hash table to map string IDs of M to their unique names
              M.space.name <- sapply(M.space.vec, paste, collapse=",")
              M.space.key <- lapply(M.space.name, function(x){prior.M.names[[x]]})
              M.space.key.null <- which(sapply(M.space.key, is.null))
              
              #compute unique names for new string IDs and update hash table
              if (length(M.space.key.null) != 0){
                M.space.key[M.space.key.null] <- sapply(M.space.vec[M.space.key.null], m.rename)
                list2env(setNames(M.space.key[M.space.key.null], M.space.name[M.space.key.null]), envir = prior.M.names)
              }
            } else {
              #compute unique names for possible settings of M
              M.space.key <- lapply(M.space.vec, m.rename)
            }
            
            #use the unique names to look up priors for possible settings of M using hash table
            M.ln.prior <- lapply(M.space.key, function(x){prior.M.hash[[x]]})
            M.ln.prior.null <- which(sapply(M.ln.prior, is.null))
            
            if (length(M.ln.prior.null)!=0){
              if (model.type=="list"){
                #values of M that are not in the hash table have probability zero
                M.ln.prior[M.ln.prior.null] <- -Inf
              } else {
                #compute mixture prior for new M
                missing.weighted.input <- sapply(M.space.key[M.ln.prior.null], function(x){prior.M.input[[x]]})
                missing.weighted.input <- sapply(missing.weighted.input, function(x){ifelse(is.null(x), -Inf, x + prior.M.weight.ln)})
                missing.weighted.crp <- sapply(M.space.vec[M.ln.prior.null], dcrp, prior.alpha=list(type="gamma", shape=prior.alpha.shape, rate=prior.alpha.rate)) + prior.M.weight.ln.1minus
                missing.mixture <- sapply(1:length(M.ln.prior.null), function(x){matrixStats::logSumExp(c(missing.weighted.input[x], missing.weighted.crp[x]))})
                
                #update hash table with mixture prior for new M
                M.ln.prior[M.ln.prior.null] <- missing.mixture
                list2env(setNames(M.ln.prior[M.ln.prior.null], M.space.key[M.ln.prior.null]), envir = prior.M.hash)
              }
            }
            M.ln.prior <- unlist(M.ln.prior)
          }
          
          #calculate t-distributed likelihood for all possible assignments of current row of M

          if (j==j.order[1]){
            M.posteriors <- lapply(MC.space, nglm.hyperparameters.ml)
            #M.posteriors <- lapply(1:(K+1), function(x){ifelse(M.ln.prior[x]==-Inf, 0, nglm.hyperparameters.ml(MC.space[x]))})
          } else {
            M.posteriors <- vector("list", K+1)
            M.posteriors[[M.current$new.index]] <- M.current$M.posteriors
            M.posteriors[-M.current$new.index] <- lapply(MC.space[-M.current$new.index], nglm.hyperparameters.ml)
            #M.posteriors[-M.current$new.index] <- lapply((1:(K+1))[-M.current$new.index], function(x){ifelse(M.ln.prior[x]==-Inf, 0, nglm.hyperparameters.ml(MC.space[x]))})
          }
          
          M.ln.ml <- unlist(lapply(M.posteriors, function(x){x$partial.ln.ml}))
          
          #combine likelihoods with priors and scale by normalizing constant
          M.prob <- M.ln.ml + M.ln.prior
          M.prob <- exp(M.prob - matrixStats::logSumExp(M.prob))
          
          #sample assignment for current row of M from categorical distribution
          M.indicator <- match(rmultinom(1,1,M.prob), x=1)
          
          if (isTRUE(M.indicator==M.current$new.index) & isTRUE(M.current$M.list[j]!=M.current$new.index)){
            M.list <- M.current$M.list
            M <- M.current$M
            K <- K + 1
          } else {
            M.list[j] <- M.indicator
            
            #update M
            if (M.indicator > ncol(M)){
              M <- cbind(M,0)
              M[j, M.indicator] <- 1
              K <- K + 1
            } else {
              M[j, M.indicator] <- 1
            }
          }
        }
        #update quantities that depend on M
        C <- contrast.list[[K]]
        MC <- M%*%C
        AMC <- A%*%MC
        d <- ncol(C)
      }
      
      #sample concentration parameter if using CRP
      if (update.alpha){
        alpha <- sample.crp.concentration()
        if (prior.alpha.type=="beta.prime"){
          prior.alpha.rate <- rgamma(1, prior.alpha.b + prior.alpha.shape, prior.alpha.q + alpha)
        }
      }
      
      #sample error variance and linear coefficients from conjugate normal-gamma distribution 
      #note that coefficients are scaled as eta=beta/lambda
      #store hyperparameters for normal-gamma distribution 
      if (update.M){
        M.posteriors <- M.posteriors[[M.indicator]]
      } else {
        M.posteriors <- nglm.hyperparameters.ml(MC)
      }
      
      #sample error variance from inv-gamma distribution
      sigma.sq.inv <- rgamma(1, 0.5*M.posteriors$kappa.star, 0.5*M.posteriors$psi.star)
      sigma.sq <- sigma.sq.inv^(-1)
      sigma <- sqrt(sigma.sq)
      
      #sample linear coefficients from normal distribution and account for scaling
      theta <- c(M.posteriors$m.star + sigma*(M.posteriors$L.Vt.inv%*%rnorm(p+d)))
      
      if (p > 0){
        theta[-(1:p)] <- theta[-(1:p)]/lambda
        delta <- theta[1:p]
        eta <- theta[-(1:p)]
      } else {
        theta <- theta/lambda
        delta <- rep(0,p)
        eta <- theta
      }
      
      #update tau.sq from conjugate inv-gamma distribution
      tau.sq.shape <- prior.phi.b + 0.5*d
      tau.sq.rate <- 0.5 + 0.5*sigma.sq.inv*sum(eta^(2))
      tau.sq <- rgamma(1, tau.sq.shape, tau.sq.rate)^(-1)
      
      #update lambda from conjugate normal distribution
      Z.delta <- Z%*%delta
      y.prime <- c(y-Z.delta)
      MCeta <- MC%*%eta
      
      if (update.M | !fixed.diplo){
        lambda.var <- (sigma.sq.inv*c(crossprod(MCeta, AtDtWDA)%*%MCeta) + 1)^(-1)
        lambda.mean <- sum(sigma.sq.inv*y.prime*W*lambda.var*c(DA%*%MCeta))
      } else {
        lambda.var <- (sigma.sq.inv*c(crossprod(eta, CtMtAtDtWDAMC)%*%eta) + 1)^(-1)
        lambda.mean <- sum(sigma.sq.inv*y.prime*W*lambda.var*c(DAMC%*%eta))
      }
      
      lambda <- rnorm(1, lambda.mean, sqrt(lambda.var))
      
      #update quantities that depend on tau.sq and lambda
      phi.sq <- lambda^2*tau.sq
      beta <- lambda*eta
      MCbeta <- lambda*(MCeta)
      AMCbeta <- A%*%MCbeta
      
      #update dipltypes jointly from independent categorical distributions
      if (!fixed.diplo){
        #calculate independent normal likelihood for all possible assignments of each row of D
        D.ln.ml <- matrix(dnorm(rep(y.prime, each=ncol.P), rep(AMCbeta, n), rep(sigma*sqrt.W^(-1), each=ncol.P), log=T), n, ncol.P, byrow=T)
        
        #combine likelihood with prior and scale to prevent underflow of all probabilities
        D.prob <- D.ln.ml + ln.P
        D.prob <- exp(D.prob - matrixStats::rowMaxs(D.prob))
        
        #sample assignment for each row of D from independent categorical distributions
        #probabilities are normalized by rmultinom
        D <- t(apply(D.prob, 1, rmultinom, n=1, size=1))
        D.list <- apply(D, 1, match, x=1)
      }
      
      #store posterior samples and hyperparameters
      post.M[i,] <- M.list
      post.MCbeta[i,] <- MCbeta
      post.delta[i,] <- delta
      post.phi.sq[i] <- phi.sq
      post.sigma.sq[i] <- sigma.sq
      post.alpha[i] <- alpha
      post.hyperparameters[[i]] <- M.posteriors[1:4]
      post.K[i] <- K
      post.y.hat[i,] <- AMCbeta[D.list] + Z.delta
      p.D.given.y <- p.D.given.y + D
    }
    
    #report unique names for posterior samples of M
    if (update.M==F){
      post.M <- rep(m.rename(post.M[1,]), iterations)
    } else if (hash.names){
      post.M <- apply(post.M, 1, function(x){prior.M.names[[paste(x, collapse=",")]]})
    } else {
      post.M <- apply(post.M, 1, m.rename)
    }
    
    #calculate marginal posterior diplotype probabilities
    p.D.given.y <- p.D.given.y/iterations
    
    #update variable state for potential reduced run of the sampler
    D <<- D
    D.list <<- D.list
    lambda <<- lambda
    tau.sq <<- tau.sq
    phi.sq <<- phi.sq
    
    #return posterior samples and hyperparameters
    posterior.results <- list("post.M"=post.M, "post.MCbeta"=post.MCbeta, "post.delta"=post.delta, "post.sigma.sq"=post.sigma.sq, "post.phi.sq"=post.phi.sq, 
                              "p.D.given.y"=p.D.given.y, "post.hyperparameters"=post.hyperparameters, "post.K"=post.K, "post.y.hat"=post.y.hat)

    if (update.alpha){
      posterior.results$post.alpha <- post.alpha
    }
      
    if (calc.null.ml){
      update.M <- T
      posterior.results$ln.ml.null <- nglm.hyperparameters.ml(matrix(1,J,1)%*%sumtozero.contrast(1), calc.partial.ml=T)$partial.ln.ml
    }
    
    posterior.results
  }
  
  ####################
  #precompute invariant hyperparameters and matrix products
  #specify starting values for sampler
  n <- length(y)
  
  #default covariate matrix includes a mean
  if (is.null(Z)){
    Z <- matrix(1,n,1)
  }
  
  #default covariance matrix assumes no replicates
  if (is.null(W)){
    W <- rep(1,n)
  }
  
  model.type <- prior.M$model.type
  fixed.diplo <- prior.D$fixed.diplo
  
  P <- prior.D$P
  A <- prior.D$A
  
  p <- ncol(Z)
  J <- ncol(A)
  ncol.P <- ncol(P)
  sqrt.W <- sqrt(W)
  ZtWZ <- crossprod(Z*sqrt.W)
  ZtW <- t(Z*W)
  Wy <- W*y
  ytWy <- sum(y*Wy)
  kappa.star <- n
  
  phi.sq <- 0.5/prior.phi.b
  lambda <- 1
  tau.sq <- phi.sq/(lambda^2)
  
  D <- matrix(0, n, ncol.P)
  D.list <- apply(P, 1, which.max)
  D[cbind(1:n, D.list)] <- 1

  if (fixed.diplo){
    DA <- A[D.list,,drop=F]
    ZtWDA <- ZtW%*%DA
    AtDtWDA <- crossprod(DA*sqrt.W)
    AtDtWy <- crossprod(DA, Wy)
    P <- D
  } else {
    ln.P <- log(P)
  }
  
  hash.names <- F
  
  if (model.type=="fixed"){
    alpha <- NA
    M <- M.matrix.from.ID(prior.M$M.IDs)
    K <- ncol(M)
    C <- sumtozero.contrast(K)
    MC <- M%*%C
    AMC <- A%*%MC
    d <- ncol(C)
    n.params <- p + d
  } else {
    contrast.list <- lapply(1:J, sumtozero.contrast)
  }
  
  if (model.type=="crp"){
    prior.alpha.type <- prior.M$prior.alpha.type
    if (prior.alpha.type=="gamma"){
      prior.alpha.shape <- prior.M$prior.alpha.shape
      prior.alpha.rate <- prior.M$prior.alpha.rate
      alpha <- prior.alpha.shape/prior.alpha.rate
    } else if (prior.alpha.type=="fixed"){
      alpha <- prior.M$prior.alpha
    } else if (prior.alpha.type=="beta.prime"){
      prior.alpha.shape <- prior.M$prior.alpha.a
      prior.alpha.b <- prior.M$prior.alpha.b
      prior.alpha.q <- ifelse(is.null(prior.M$prior.alpha.q), 1, prior.M$prior.alpha.q)
      prior.alpha.rate <- prior.alpha.b/prior.alpha.q
      alpha <- prior.alpha.shape/prior.alpha.rate
    }
    M <- matrix(1, J, 1)
  } else if (model.type=="uniform"){
    alpha <- NA
    M <- matrix(1, J, 1)
  } else if (model.type=="list" | model.type=="mixture"){
    prior.M.hash <- new.env(hash = T)
    alpha <- NA
    M <- M.matrix.from.ID(prior.M$M.IDs[which.max(prior.M$ln.probs)])
    
    if (model.type=="list"){
      list2env(setNames(as.list(prior.M$ln.probs), prior.M$M.IDs), envir = prior.M.hash)
    } else {
      prior.M.input <- new.env(hash = T)
      list2env(setNames(as.list(prior.M$ln.probs), prior.M$M.IDs), envir = prior.M.input)
      prior.alpha.shape <- prior.M$prior.alpha.shape
      prior.alpha.rate <- prior.M$prior.alpha.rate
      prior.M.weight.ln <- log(prior.M$weight)
      prior.M.weight.ln.1minus <- log(1-prior.M$weight)
    }
    
    hash.names <- prior.M$hash.names
    if (hash.names){
      prior.M.names <- new.env(hash = T)
    }
  }
  
  M.list <- apply(M, 1, match, x=1)
  
  ####################
  #iterate TIMBR sampler
  if (verbose){
    print("Sampling from the full posterior")
  }
  
  if (model.type=="crp"){
    if (prior.alpha.type=="gamma" | prior.alpha.type=="beta.prime"){
      results <- TIMBR.sampler(samples)
    } else if (prior.alpha.type=="fixed"){
      results <- TIMBR.sampler(samples, update.alpha=F)
    }
  } else if (model.type=="fixed"){
    results <- TIMBR.sampler(samples, update.M=F, update.alpha=F)
  } else if (model.type=="uniform" | model.type=="list" | model.type=="mixture"){
    results <- TIMBR.sampler(samples, update.alpha=F)
  }
  
  #collect posterior probability of the null model
  ln.ml.null <- results$ln.ml.null
  
  post.M.ranked <- -sort(-table(results$post.M))/samples
  post.M.null <- post.M.ranked[names(post.M.ranked)==paste(rep(0, J), collapse=",")]
  
  if(length(post.M.null)==0){
    post.M.null <- 0
  }
  
  #if posterior probability of null model is relatively high, use this to calculate marginal likelihood
  if (names(post.M.ranked[1])==paste(rep(0, J), collapse=",") | post.M.null>=0.01){
    if (model.type=="crp"){
      if (prior.alpha.type=="gamma"){
        ln.ml <- ln.ml.null + dcrp(rep(0,J), list(type="gamma", shape=prior.alpha.shape, rate=prior.alpha.rate)) - log(post.M.null)
      } else if (prior.alpha.type=="fixed"){
        ln.ml <- ln.ml.null + dcrp(rep(0,J), list(type="fixed", alpha=prior.M$prior.alpha)) - log(post.M.null)
      } else if (prior.alpha.type=="beta.prime"){
        ln.ml <- ln.ml.null + dcrp(rep(0,J), list(type="beta.prime", a=prior.alpha.shape, b=prior.alpha.b, q=prior.alpha.q)) - log(post.M.null)
      }
    } else if (model.type=="fixed"){
      ln.ml <- ln.ml.null
    } else if (model.type=="uniform"){
      ln.ml <- ln.ml.null - ln.bell(J) - log(post.M.null)
    } else if (model.type=="list" | model.type=="mixture"){
      ln.ml <- ln.ml.null + prior.M.hash[[paste(rep(0,J), collapse=",")]] - log(post.M.null)
    }
  } else {
    #obtain samples of the normal-gamma hyperparameters conditional on a single model matrix M
    if (model.type=="fixed"){
      #these are given when the model.type is fixed
      nglm.hyperparameters <- results$post.hyperparameters
      post.phi.sq <- results$post.phi.sq
      p.D.given.y <- results$p.D.given.y
      samples.ml <- samples
    } else {
      #set M to the MAP and update dependent quantities
      M <- M.matrix.from.ID(names(post.M.ranked)[1])
      M.list <- apply(M, 1, match, x=1)
      K <- ncol(M)
      C <- contrast.list[[K]]
      MC <- M%*%C
      AMC <- A%*%MC
      d <- ncol(C)
      n.params <- p + d
      
      #iterate reduced sampler run with M fixed at MAP
      if (verbose){
        print("Sampling from the conditional posterior")
      }
      
      reduced.results <- TIMBR.sampler(samples.ml, update.M=F, calc.null.ml=F, update.alpha=F)
      nglm.hyperparameters <- reduced.results$post.hyperparameters
      post.phi.sq <- reduced.results$post.phi.sq
      p.D.given.y <- reduced.results$p.D.given.y
    }
    
    #set remaining variables to values with high posterior probability given M fixed at MAP
    kappa.star <- n
    psi.star <- sapply(nglm.hyperparameters, function(x){x$psi.star})
    m.star <- lapply(nglm.hyperparameters, function(x){x$m.star})
    V.star <- lapply(nglm.hyperparameters, function(x){x$V.star})
    
    theta <- matrixStats::rowMeans2(sapply(1:length(m.star), function(x){sapply(unlist(m.star[x]), function(x){x})}))
    
    if (p > 0){
      delta <- theta[1:p]
      beta <- theta[-(1:p)]
    } else {
      delta <- rep(0,p)
      beta <- theta
    }
    
    sigma.sq <- mean((kappa.star/psi.star)^(-1))
    
    if (!fixed.diplo){
      D <- matrix(0, n, ncol.P)
      D.list <- apply(p.D.given.y, 1, which.max)
      D[cbind(1:n, D.list)] <- 1
    }
    
    #calculate partial marginal likelihood at point of high posterior probability
    Zdelta <- Z%*%delta
    AMCbeta <- A%*%MC%*%beta
    
    p1 <- sum(dnorm(y, AMCbeta[D.list] + Zdelta, sqrt(sigma.sq*W^(-1)), log=T))
    p4 <- log(sigma.sq)
    p5 <- ln.beta.prior.marginalized(beta, sigma.sq, prior.phi.b)
    p7 <- matrixStats::logSumExp(dgamma(sigma.sq^(-1), 0.5*kappa.star, 0.5*psi.star, log=T))-log(samples.ml)
    p8 <- matrixStats::logSumExp(sapply(1:samples.ml, function(x){mvtnorm::dmvnorm(theta, m.star[[x]], V.star[[x]]*sigma.sq, log=T)}))-log(samples.ml)
    c <- -0.5*n*log(pi) + lgamma(0.5*kappa.star) - 0.5*sum(-log(W))
    
    if (!fixed.diplo){
      p2 <- sum(ln.P[cbind(1:n, D.list)])
      
      D.ln.ml <- matrix(dnorm(rep((y - Zdelta), each=ncol.P), rep(AMCbeta, n), rep(sqrt(sigma.sq)*sqrt.W^(-1), each=ncol.P), log=T), n, ncol.P, byrow=T)
      D.ln.prob <- D.ln.ml + ln.P
      D.ln.prob <-  D.ln.prob - matrixStats::rowLogSumExps(D.ln.prob)
      p9 <- sum(D.ln.prob[cbind(1:n, D.list)])
    } else {
      p2 <- 0
      p9 <- 0
    }
    
    if (model.type=="crp"){
      if (prior.alpha.type=="gamma"){
        p3 <- dcrp(apply(M, 1, match, x=1), list(type="gamma", shape=prior.alpha.shape, rate=prior.alpha.rate))
      } else if (prior.alpha.type=="fixed"){
        p3 <- dcrp(apply(M, 1, match, x=1), list(type="fixed", alpha=prior.M$prior.alpha))
      } else if (prior.alpha.type=="beta.prime"){
        p3 <- dcrp(apply(M, 1, match, x=1), list(type="beta.prime", a=prior.alpha.shape, b=prior.alpha.b, q=prior.alpha.q))
      }
      p6 <- log(post.M.ranked[1])
    } else if (model.type=="fixed"){
      p3 <- 0
      p6 <- 0
    } else if (model.type=="uniform"){
      p3 <- -ln.bell(J)
      p6 <- log(post.M.ranked[1])
    } else if (model.type=="list" | model.type=="mixture"){
      p3 <- prior.M.hash[[names(post.M.ranked[1])]]
      p6 <- log(post.M.ranked[1])
    }
    
    names(p6) <- NULL
    ln.ml <- p1 + p2 + p3 + p4 + p5 - p6 - p7 - p8 - p9 - c
  }
  
  #return posterior samples, marginal densities, and marginal likelihood
  output <- list(y=y, prior.D=prior.D, prior.M=prior.M, prior.phi.b=prior.phi.b, samples=samples, samples.ml=samples.ml, Z=Z, W=W)
  output <- c(output, results[names(results) != "post.hyperparameters"], ln.ml=ln.ml, ln.BF=ln.ml-ln.ml.null, p.M.given.y=list(post.M.ranked), 
              post.var.exp=list(results$post.phi.sq/(results$post.phi.sq + 1)), post.hap.effects=list(results$post.MCbeta+results$post.delta[,1]),
              p.K.given.y=list(table(results$post.K)/samples))
  output <- output[order(names(output))]
}

#' @keywords internal
stirling.first.unsigned <- function(J){
  s <- diag(1, J)
  
  if (J>1){
    for (n in 2:J){
      for (k in 1:n){
        s[n,k] <- ifelse(k == 1, (n-1)*s[n-1, k], (n-1)*s[n-1, k] + s[n-1, k-1])
      }
    }
  }
  
  s[J,]
}

#' @keywords internal
ln.K.prior.crp.marginalized <- function(K, J, a, b){
  s <- stirling.first.unsigned(J)
  
  density.K.concentration <- Vectorize(function(x){
    exp(K*log(x) + lgamma(x) - lgamma(x + J) + dgamma(x, a, b, log=T))
  })
  
  log(integrate(density.K.concentration, lower=0, upper=Inf)$value) + log(s[K]) - log(sum(s)) + lfactorial(J)
}

#' Calculate hyperparameters for Chinese restaurant process with gamma prior distribution on the concentration parameter
#'
#' Uses optimization to find hyperparameters that approximately yield the specified prior probabilities. May not work well when J is large or target prior probabilities are extreme.
#'
#' @param J number of customers (haplotypes) to be partitioned
#' @param p.1.target prior probability of 1 partition
#' @param p.J.target prior probability of J partitions
#'
#' @return vector of c(shape, rate) hyperparamters for the gamma distribution
#' 
#' @examples
#' calc.concentration.prior(8, 0.05, 0.01)
#'
#' @export
calc.concentration.prior <- function(J, p.1.target, p.J.target){
  distance <- function(c){
    p1 <- exp(ln.K.prior.crp.marginalized(1, J, c[1], c[2]))
    pJ <- exp(ln.K.prior.crp.marginalized(J, J, c[1], c[2]))
    
    (1/2)*log((p1-p.1.target)^2 + (pJ-p.J.target)^2)
  }
  
  optim(c(1,1), distance)$par
}

#' Ewens's sampling formula with optional gamma prior on the concentration parameter
#'
#' Sample allelic series (paritions) from Ewen's sampling formula, optionally informed by user-specified tree(s). Trees must be in coalescent units for appropriate inference.
#'
#' @param samples number of samples
#' @param trees either a user-specified tree(s) of class "phylo" ("multiPhylo"), detailed in the 'ape' package, or an integer with the number of leaves to be partitioned
#' @param prior.alpha prior type c("fixed","gamma") for the concentration parameter, see examples for format. type "beta.prime" supported but sampler appears unstable, use ewenss.calc function instead
#' @param verbose optionally report function progress
#'
#' @return list of allelic series IDs and probabilities, formatted as prior.M object for TIMBR function
#' 
#' @examples
#' #specifying hyperparameters for gamma prior using calc.concentration.prior
#' hyperparam <- calc.concentration.prior(8, 0.05, 0.01)
#' prior.alpha <- list(type="gamma", shape=hyperparam[1], rate=hyperparam[2])
#' 
#' #running the sampler without user-specified trees; compare with target prior probabilities
#' prior.M <- ewenss.sampler(100000, 8, prior.alpha)
#' exp(prior.M$ln.probs[prior.M$M.IDs=="0,0,0,0,0,0,0,0"])
#' exp(prior.M$ln.probs[prior.M$M.IDs=="0,1,2,3,4,5,6,7"])
#' 
#' #running the sampler with a user-specified tree and fixed concentration parameter; compare with tree structure
#' tree <- ape::rcoal(8, LETTERS[1:8])
#' ape::plot.phylo(tree)
#' prior.alpha <- list(type="fixed", alpha=1)
#' prior.M <- ewenss.sampler(100000, tree, prior.alpha)
#' head(prior.M$M.IDs)
#' head(exp(prior.M$ln.probs))
#'
#' @export
ewenss.sampler <- function(samples, trees, prior.alpha, verbose=T){
  sample.M.ID.from.tree <- function(iter, tree=NULL){
    #sample coalescent tree if tree is unspecified
    if (is.null(tree)){
      tree <- ape::rcoal(J, LETTERS[1:J])
    }
    
    #decompose tree into basis V and lengths l
    tree.decomposed <- decompose.tree(tree)
    V <- tree.decomposed$V
    l <- tree.decomposed$l
    L <- sum(l)
    
    #optionally sample the mutation rates for the poisson process
    if (prior.alpha$type!="fixed"){
      if (prior.alpha$type=="beta.prime"){
        prior.alpha.rate <- rgamma(1, prior.alpha.b, prior.alpha.q)
      }
      alpha <- rgamma(iter, prior.alpha.shape, prior.alpha.rate)
    }
    lambda <- alpha*L/2
    
    #calculate null probabilities and sample non-zero numbers of mutations from truncated poisson
    p.null <- dpois(0, lambda)
    theta <- rpois(iter, lambda + log(1 - runif(iter)*(1 - exp(-lambda)))) + 1
    
    #probabilistically assign mutations to branches in proportion to branch lengths
    B <- rbind(sapply(theta, rmultinom, n=1, prob=l[-length(l)])>0, T)
    
    #collect allelic series
    m <- apply(B, 2, function(y){m.rename(apply(V[,y], 1, Position, f=function(x){x==1}))})
    
    #return allelic series and null probabilities
    list(m, p.null)
  }
  
  #optionally disable reporting
  if (verbose){
    print("Iterating Ewens's sampling formula")
  }
  
  #specify alpha or prior hyperparameters
  if (prior.alpha$type=="gamma"){
    prior.alpha.shape <- prior.alpha$shape
    prior.alpha.rate <- prior.alpha$rate
  } else if (prior.alpha$type=="fixed"){
    alpha <- prior.alpha$alpha
  } else if (prior.alpha$type=="beta.prime"){
    prior.alpha.shape <- prior.alpha$a
    prior.alpha.b <- prior.alpha$b
    prior.alpha.q <- ifelse(is.null(prior.alpha$q), 1, prior.alpha$q)
  }
  
  #iterate Ewens's sampling formula using specified trees
  if (is.numeric(trees)){
    #sample allelic series from random coalescent trees
    J <- trees
    results <- replicate(samples, sample.M.ID.from.tree(1))
  } else {
    #set class to multiphylo if single tree is specified
    if (class(trees)=="phylo"){
      trees <- c(trees)
    }
    
    #store number of leaves on the tree
    J <- trees[[1]]$Nnode + 1
    
    #calculate required number of samples for each specified tree
    n.trees <- length(trees)
    iter <- rep(floor(samples/n.trees), n.trees)
    remainder <- samples%%n.trees
    if (remainder != 0){
      iter[1:remainder] <- iter[1:remainder] + 1
    }
    
    #sample allelic series from specified trees
    results <- sapply(1:n.trees, function(x){sample.M.ID.from.tree(iter[x], trees[[x]])})
  }
  
  df <- data.frame(M.IDs=unlist(results[1,]), wt=(1-unlist(results[2,]))/samples, stringsAsFactors=F)
  df <- dplyr::count(df, M.IDs, wt = wt)
  df <- dplyr::add_row(df, M.IDs=paste(rep(0,J),collapse=","), n=1-sum(df$n))
  df <- dplyr::arrange(df, dplyr::desc(n))
  
  list(model.type="list", M.IDs=df$M.IDs, ln.probs=log(df$n), hash.names=T)
}

#' Create an additive genetic design matrix
#'
#' Returns the additive design matrix for the number of specified haplotypes, in the specified format.
#'
#' @param J number of haplotypes
#' @param type format of the diplotype matrix, currently only supports "happy"
#'
#' @return matrix of (half) haplotype dosages corresponding to each diplotype state
#' 
#' @examples
#' additive.design(8, "happy")
#'
#' @export
additive.design <- function(J, type){
  A <- diag(J)
  
  if (type=="happy"){
    for (j in 2:J){
      for (k in 1:(j-1)){
        A <- rbind(A, 0)
        A[nrow(A),k] <- 0.5
        A[nrow(A),j] <- 0.5
      }
    }    
  }
  
  A
}

#' Plot Haplotype Effects from TIMBR Output
#'
#' Plots posterior haplotype effect densities from TIMBR output
#'
#' @param TIMBR.output results object from the TIMBR function
#' @param colors an optional vector of colors for each haplotype density
#' @param file.path an optional file path for saving the plot as a PNG
#' @param plot.width PNG plot width
#' @param plot.height PNG plot height
#' @param TIMBR.output.bkgrd optional second results object for background plot
#' @param colors.bkgrd optional colors for background plot
#' @param transparency optional vector of transparencies for foreground and background plots, ignored if second results object is not specified
#'
#' @return plot of the posterior haplotype effect densities
#' 
#' @examples
#' #example data
#' data(mcv.data)
#' str(mcv.data)
#' 
#' #call TIMBR using CRP
#' results <- TIMBR(mcv.data$y, mcv.data$prior.D, mcv.data$prior.M$crp)
#' 
#' #plot haplotype effects
#' TIMBR.plot.haplotypes(results)
#'
#' @export
TIMBR.plot.haplotypes <- function(TIMBR.output, colors=NULL, file.path=NULL, plot.width=960, plot.height=480, TIMBR.output.bkgrd=NULL,
                                  colors.bkgrd=NULL, transparency=c(0.7,0.4)){
  densities <- apply(TIMBR.output$post.hap.effect, 2, density)
  J <- length(densities)
  
  scale.y <- max(sapply(densities, function(x){max(x$y)})*1.2)
  min.x <- min(sapply(densities, function(x){min(x$x)}))
  max.x <- max(sapply(densities, function(x){max(x$x)}))
  
  if (!is.null(TIMBR.output.bkgrd)){
    densities.bkgrd <- apply(TIMBR.output.bkgrd$post.hap.effect, 2, density)
    if(length(densities.bkgrd)!=J){
      stop("J does not match for TIMBR.output and TIMBR.output.bkgrd")
    } 
    scale.y <- max(sapply(densities.bkgrd, function(x){max(x$y)})*1.2, scale.y)
    min.x <- min(sapply(densities.bkgrd, function(x){min(x$x)}), min.x)
    max.x <- max(sapply(densities.bkgrd, function(x){max(x$x)}), max.x)
  }
  
  par(cex.axis=1.1, cex.lab=1.1, cex.main=1.2, cex.sub=1.1)
  
  if (!is.null(file.path)){
    png(file.path, height=plot.height, width=plot.width)
  }
  
  plot(1, type="n", xlab="Phenotype", ylab="Haplotype", xlim=c(min.x, max.x), ylim=c(0, scale.y*J), axes=FALSE)
  Axis(side=2, at=scale.y*(0:(J-1)+0.5), labels=LETTERS[J:1], las=1, tick=F)
  Axis(side=1)
  
  if (is.null(colors)){
    colors <- ifelse(J==8, c("#9000E0","#F00000","#00A000","#00A0F0","#1010F0","#F08080","#808080","#F0F000"), rep("#4D4D4D", J))
  }
  
  
  if (!is.null(TIMBR.output.bkgrd)){
    if (is.null(colors.bkgrd)){
      colors.bkgrd <- scales::alpha(colors, transparency[2])
      colors <- scales::alpha(rep("#4D4D4D", J), transparency[1])
    }
    
    for (i in 1:J){
      densities.bkgrd[[(J+1)-i]]$y <- densities.bkgrd[[(J+1)-i]]$y + (i-1)*scale.y
      polygon(densities.bkgrd[[(J+1)-i]], col=colors.bkgrd[i], lwd=1)
    }
  }
  
  for (i in 1:J){
    densities[[(J+1)-i]]$y <- densities[[(J+1)-i]]$y + (i-1)*scale.y
    polygon(densities[[(J+1)-i]], col=colors[i], lwd=1)
  }
  
  if (!is.null(file.path)){
    dev.off()
  }
}

#' @keywords internal
consistency.index <- function(J, return.setparts=F){
  partitions.all <- partitions::setparts(J)
  colnames(partitions.all) <- apply(partitions.all, 2, m.rename)
  
  M1 <- apply(partitions.all, 2, function(x){M <- matrix(0, length(x), max(x)); M[cbind(1:8, x)] <- 1; MMt <- tcrossprod(M); MMt[upper.tri(MMt)]})
  M0 <- M1[,apply(partitions.all, 2, max)==2]
  
  index <- apply(M1, 2, function(y){apply(M0, 2, function(x){match(-1, x-y, 0)==0})})
  
  if (return.setparts){
    list(index=index, setparts=partitions.all)
  } else {
    index
  }
}

#' Approximate Bayes Factors from TIMBR Output
#'
#' Approximate Bayes factors for various hypotheses from TIMBR output
#'
#' @param TIMBR.output results object from the TIMBR function, prior.M method must be CRP or uniform
#' @param type "all" - report BFs for all allelic series; "merge" - report BFs for biallelic series (i.e. merge analysis); "consistent" - report BFs for consistent merge analysis; prior.M - report BF for list-type prior.M
#'
#' @return a named vector of approximate Bayes Factors; appoximations that underflow are not reported
#' 
#' @examples
#' #example data
#' data(mcv.data)
#' str(mcv.data)
#' 
#' #call TIMBR using CRP
#' results <- TIMBR(mcv.data$y, mcv.data$prior.D, mcv.data$prior.M$crp)
#' 
#' #calculate biallelic consistency
#' TIMBR.approx(results)
#'
#' @export
TIMBR.approx <- function(TIMBR.output, type="all", ln.ml = F, return.prior=F){
  if (TIMBR.output$prior.M$model.type=="crp"){
    if (TIMBR.output$prior.M$prior.alpha.type=="gamma"){
      prior.alpha <- list(type="gamma", shape=TIMBR.output$prior.M$prior.alpha.shape, rate=TIMBR.output$prior.M$prior.alpha.rate)
    } else if (TIMBR.output$prior.M$prior.alpha.type=="fixed"){
      prior.alpha <- list(type="fixed", alpha=TIMBR.output$prior.M$prior.alpha)
    } else if (TIMBR.output$prior.M$prior.alpha.type=="beta.prime"){
      if (is.null(TIMBR.output$prior.M$prior.alpha.q)){
        TIMBR.output$prior.M$prior.alpha.q <- 1
      }
      prior.alpha <- list(type="beta.prime", a=TIMBR.output$prior.M$prior.alpha.a, b=TIMBR.output$prior.M$prior.alpha.b, q=TIMBR.output$prior.M$prior.alpha.q)
    }
  } else if (TIMBR.output$prior.M$model.type=="uniform"){
    ln.prior.uniform <- -ln.bell(ncol(TIMBR.output$prior.D$A))
  } else if (TIMBR.output$prior.M$model.type!="list" & type=="consistent"){
    stop("prior.M$model.type is not crp/uniform or list (consistent type only)")
  }
  
  if (class(type)!="list"){
    if (type=="all"){
      if (TIMBR.output$prior.M$model.type=="crp"){
        ln.prior <- sapply(names(TIMBR.output$p.M.given.y), function(x){dcrp(m.from.M.ID(x), prior.alpha)})
      } else if (TIMBR.output$prior.M$model.type=="uniform"){
        ln.prior <- rep(ln.prior.uniform, length(TIMBR.output$p.M.given.y))
        names(ln.prior) <- names(TIMBR.output$p.M.given.y)
      }
      
      ln.BF <- rev(sort(TIMBR.output$ln.BF + log(TIMBR.output$p.M.given.y) - ln.prior))
    } else if (type=="merge"){
      m.list <- lapply(names(TIMBR.output$p.M.given.y), m.from.M.ID)
      biallelic <- sapply(m.list, max)==2
      
      if (TIMBR.output$prior.M$model.type=="crp"){
        ln.prior <- sapply(m.list[biallelic], dcrp, prior.alpha=prior.alpha)
      } else if (TIMBR.output$prior.M$model.type=="uniform"){
        ln.prior <- rep(ln.prior.uniform, length(m.list[biallelic]))
        names(ln.prior) <- names(m.list)[biallelic]
      }
      
      ln.BF <- rev(sort(TIMBR.output$ln.BF + log(TIMBR.output$p.M.given.y[biallelic]) - ln.prior))
    } else if (type=="consistent"){
      index <- consistency.index(ncol(TIMBR.output$prior.D$A), T)
      if (TIMBR.output$prior.M$model.type=="crp"){
        ln.prior <- apply(index$setparts, 2, dcrp, prior.alpha=prior.alpha)
      } else if (TIMBR.output$prior.M$model.type=="uniform"){
        ln.prior <- rep(ln.prior.uniform, ncol(index$setparts))
        names(ln.prior) <- colnames(index$setparts)
      } else if (TIMBR.output$prior.M$model.type=="list"){
        ln.prior <- TIMBR.output$prior.M$ln.probs
        names(ln.prior) <- TIMBR.output$prior.M$M.IDs
        ln.prior <- ln.prior[colnames(index$index)]
        ln.prior[is.na(ln.prior)] <- -Inf
        names(ln.prior) <- colnames(index$index)
      }
      
      index.by.ln.prior <- as.matrix(index$index)%*%diag(ln.prior)
      index.by.ln.prior <- apply(index.by.ln.prior, 2, function(x){x[x==0 | is.nan(x)] <- -Inf; x})
      index.by.ln.prior <- index.by.ln.prior - matrixStats::rowLogSumExps(index.by.ln.prior)
      colnames(index.by.ln.prior) <- colnames(index$index)
      index.by.ln.prior <- index.by.ln.prior[,names(TIMBR.output$p.M.given.y)]
      ln.prior <- ln.prior[names(TIMBR.output$p.M.given.y)]
      
      ln.BF <- rev(sort(TIMBR.output$ln.BF + apply(index.by.ln.prior, 1, function(x){matrixStats::logSumExp(x + log(TIMBR.output$p.M.given.y) - ln.prior)})))
    }
  } else {
    if (TIMBR.output$prior.M$model.type=="crp"){
      ln.prior <- apply(sapply(names(TIMBR.output$p.M.given.y), TIMBR:::m.from.M.ID), 2, TIMBR:::dcrp, prior.alpha=prior.alpha)
    } else if (TIMBR.output$prior.M$model.type=="uniform"){
      ln.prior <- rep(ln.prior.uniform, length(TIMBR.output$p.M.given.y))
      names(ln.prior) <- colnames(TIMBR.output$p.M.given.y)
    }
    
    BFs <- rev(sort(TIMBR.output$ln.BF + log(TIMBR.output$p.M.given.y) - ln.prior))
    BFs <- BFs[prior.M$M.IDs]
    BFs[is.na(BFs)] <- -Inf
    names(BFs) <- prior.M$M.IDs
    
    ln.BF <- matrixStats::logSumExp(BFs + prior.M$ln.probs)
  }
  
  if (ln.ml){
    ln.BF <- ln.BF + TIMBR.output$ln.ml.null
  }
  
  if (return.prior){
    output <- cbind(ln.BF, ln.prior[names(ln.BF)])
    colnames(output)[2] <- "ln.prior"
    output
  } else {
    ln.BF
  }
}

#' @keywords internal
decompose.tree <- function(tree){
  J <- length(tree$tip.label)
  
  #store basis matrix V
  V <- matrix(0, J, 2*J-1)
  rownames(V) <- tree$tip.label
  nodepaths <- ape::nodepath(tree)
  V[cbind(unlist(lapply(1:J, function(i) rep(i, length(nodepaths[[i]])))), unlist(nodepaths))] <- 1
  colSums.V <- matrixStats::colSums2(V)
  
  #store branch lengths l
  l <- c(tree$edge.length, 0)
  l <- l[order(c(tree$edge[,2], Position(function(x){x==J}, colSums.V)))]
  
  #reorder V and l
  order.colSums.V <- order(colSums.V)
  l <- l[order.colSums.V]
  V <- V[sort(rownames(V)), order.colSums.V]
  
  #return V and l
  list(V=V, l=l)
}

#' @keywords internal
dcrp <- function(m, prior.alpha, log.p=T){
  J <- length(m)
  J.k <- table(m, dnn=NULL)
  K <- length(J.k)
  
  if (prior.alpha$type=="gamma"){
    shape <- prior.alpha$shape
    rate <- prior.alpha$rate
    
    density.crp.concentration <- Vectorize(function(x){
      exp(lgamma(x) - lgamma(x+J) + (shape+K-1)*log(x) - rate*x)
    })
    
    ln.p <- log(integrate(density.crp.concentration, lower=0, upper=Inf)$value) + sum(lgamma(J.k)) + shape*log(rate) - lgamma(shape)
    
  } else if (prior.alpha$type=="fixed"){
    alpha <- prior.alpha$alpha
    
    ln.p <- lgamma(alpha) - lgamma(alpha+J) + K*log(alpha) + sum(lgamma(J.k))
  } else if (prior.alpha$type=="beta.prime"){
    a <- prior.alpha$a
    b <- prior.alpha$b
    q <- ifelse(is.null(prior.alpha$q), 1, prior.alpha$q)
    
    density.crp.concentration <- Vectorize(function(x){
      exp(lgamma(x) - lgamma(x+J) + (a+K-1)*log(x) - (a+b)*log(1+x/q))
    })
    
    ln.p <- log(integrate(density.crp.concentration, lower=0, upper=Inf)$value) + sum(lgamma(J.k)) - a*log(q) - lbeta(a,b)
  }
  
  ifelse(log.p, ln.p, exp(ln.p))
}

#' Direct calculation of Ewens's sampling formula with optional gamma or beta prime prior on the concentration parameter
#'
#' Directly calculates probabilties for allelic series (paritions) under Ewen's sampling formula, optionally informed by a user-specified tree. Trees must be in coalescent units for appropriate inference.
#'
#' @param trees either a user-specified tree of class "phylo" ("multiPhylo"), detailed in the 'ape' package, or an integer with the number of leaves to be partitioned
#' @param prior.alpha prior type c("fixed","gamma","beta.prime") for the concentration parameter, see examples for format
#'
#' @return list of allelic series IDs and probabilities, formatted as prior.M object for TIMBR function
#' 
#' @examples
#' #specifying hyperparameters for gamma prior using calc.concentration.prior
#' hyperparam <- calc.concentration.prior(8, 0.05, 0.01)
#' prior.alpha <- list(type="gamma", shape=hyperparam[1], rate=hyperparam[2])
#' 
#' #running the sampler without user-specified trees; compare with target prior probabilities
#' prior.M <- ewenss.calc(8, prior.alpha)
#' exp(prior.M$ln.probs[prior.M$M.IDs=="0,0,0,0,0,0,0,0"])
#' exp(prior.M$ln.probs[prior.M$M.IDs=="0,1,2,3,4,5,6,7"])
#' 
#' #running the sampler with a user-specified tree and fixed concentration parameter; compare with tree structure
#' tree <- ape::rcoal(8, LETTERS[1:8])
#' ape::plot.phylo(tree)
#' prior.alpha <- list(type="fixed", alpha=1)
#' prior.M <- ewenss.calc(tree, prior.alpha)
#' head(prior.M$M.IDs)
#' head(exp(prior.M$ln.probs))
#'
#' @export
ewenss.calc <- function(tree, prior.alpha){
  ln.prob.and.M.ID.from.B.ID <- function(B.ID){
    #function to calculate probability of branch mutation configuration
    B <- as.logical(intToBits(B.ID-1)[1:(ncol(V)-1)])
    M.ID <- m.rename(apply(V[, c(B, T), drop=F], 1, Position, f=function(x){x==1}))
    
    if (prior.alpha$type=="fixed"){
      ln.prob <- sum(ln.p[cbind(1:nrow(ln.p), as.integer(B)+1)])
    } else if (prior.alpha$type=="gamma"){
      combinations <- as.matrix(expand.grid(replicate(sum(B), 0:1, simplify=FALSE)))
      sign <- ifelse(B.ID==1, 1, (-1)^(apply(combinations, 1, sum)%%2))
      params <- matrix(TRUE, 2^sum(B), length(l))
      params[, which(B)] <- as.logical(combinations)
      half.l <- 0.5*l
      b.prime <- prior.alpha.rate + apply(params, 1, function(x){sum(half.l[x])})
      
      ln.prob <- prior.alpha.shape*log(prior.alpha.rate) + log(sum(sign*b.prime^(-prior.alpha.shape)))
    } else if (prior.alpha$type=="beta.prime"){
      density.ewens.beta.prime <- Vectorize(function(x){
        ln.p <- -x*l[-length(l)]/2
        ln.p <- cbind(ln.p, log(1-exp(ln.p)))

        exp(sum(ln.p[cbind(1:nrow(ln.p), as.integer(B)+1)]))*x^(prior.alpha.a-1)*(1+x/prior.alpha.q)^(-prior.alpha.a-prior.alpha.b)
      })
      
      ln.prob <- log(integrate(density.ewens.beta.prime, lower=0, upper=Inf, abs.tol=0)$value) - lbeta(prior.alpha.a, prior.alpha.b) - log(prior.alpha.q) - (prior.alpha.a-1)*log(prior.alpha.q)
    }
    
    c(M.ID, ln.prob)
  }
  
  if (is.numeric(tree)){
    #enumerate all set partitions by partition class
    partitions.all <- apply(partitions::parts(tree), 2, partitions::setparts)
    
    #calculate probabilities for each partition class
    ln.probs <- unlist(sapply(1:length(partitions.all), function(x){rep(dcrp(partitions.all[[x]][,1], prior.alpha), ncol(partitions.all[[x]]))}))
    
    #normalize total to correct for approximation
    if (prior.alpha$type!="fixed"){
      ln.probs <- ln.probs - matrixStats::logSumExp(ln.probs)
    }
    
    #generate partition names and prior.M object
    M.IDs <- apply(do.call(cbind, partitions.all), 2, m.rename)
    list(model.type="list", M.IDs=M.IDs, ln.probs=ln.probs, hash.names=T)
  } else {
    #decompose tree into basis V and lengths l
    tree.decomposed <- decompose.tree(tree)
    V <- tree.decomposed$V
    l <- tree.decomposed$l
    
    if (prior.alpha$type=="fixed"){
      #store mutation probabilities for each branch
      ln.p <- -prior.alpha$alpha*l[-length(l)]/2
      ln.p <- cbind(ln.p, log(1-exp(ln.p)))
    } else if (prior.alpha$type=="gamma"){
      prior.alpha.shape <- prior.alpha$shape
      prior.alpha.rate <- prior.alpha$rate
    } else if (prior.alpha$type=="beta.prime"){
      prior.alpha.a <- prior.alpha$a
      prior.alpha.b <- prior.alpha$b
      prior.alpha.q <- ifelse(is.null(prior.alpha$q), 1, prior.alpha$q)
    }
    
    #calculate probabilties for all combinations of branch mutations and collapse by M.ID
    df <- do.call(rbind, lapply(1:(2^(ncol(V)-1)), ln.prob.and.M.ID.from.B.ID))
    df <- data.frame(M.IDs=df[,1], ln.probs=as.numeric(df[,2]), stringsAsFactors=F)
    df <- dplyr::summarize(dplyr::group_by(df, M.IDs), ln.probs = matrixStats::logSumExp(ln.probs))
    df <- dplyr::arrange(df, dplyr::desc(ln.probs))
    
    #normalize total to correct for approximation
    if (prior.alpha$type=="beta.prime"){
      df$ln.probs <- df$ln.probs - matrixStats::logSumExp(df$ln.probs)
    }
    
    list(model.type="list", M.IDs=df$M.IDs, ln.probs=df$ln.probs, hash.names=T)
  }
}

#' Construct consistent prior from existing list-type prior.M
#'
#' Updates prior.M to be consistent with biallelic contrast encoded by M.ID
#'
#' @param prior.M prior.M object for use as TIMBR input, model.type must be list
#' @param M.ID string that denotes the specified biallelic contrast
#'
#' @return updated prior.M object that is consistent with M.ID
#' 
#' @examples
#' #example data
#' data(mcv.data)
#' 
#' #call TIMBR using CRP
#' results <- TIMBR(mcv.data$y, mcv.data$prior.D, mcv.data$prior.M$crp)
#' 
#' #approximate biallelic consistency BFs for all contrasts
#' head(TIMBR.approx(results))
#' 
#' #generate list object for CRP prior
#' prior.alpha <- list(type="gamma", shape=results$prior.M$prior.alpha.shape, rate=results$prior.M$prior.alpha.rate)
#' prior.M <- ewenss.calc(8, prior.alpha)
#' 
#' #specify biallelic contrast for consistency and update prior.M
#' M.ID <- "0,1,0,1,1,0,0,0"
#' prior.M <- TIMBR.consistent(prior.M, M.ID)
#' 
#' #analyze results, compare with approximate BF from earlier
#' results.consistent <- TIMBR(mcv.data$y, mcv.data$prior.D, prior.M)
#' results.consistent$ln.BF
#'
#' @export
TIMBR.consistent <- function(prior.M, M.ID){
  m <- m.from.M.ID(M.ID)
  
  if (max(m)!=2){
    stop("M.ID is not biallelic")
  }
  
  J <- length(m)
  
  index <- consistency.index(J, F)
  index <- index[which(row.names(index)==M.ID),]
  index <- index[prior.M$M.IDs]
  
  prior.M$M.IDs <- prior.M$M.IDs[index]
  
  ln.probs <- prior.M$ln.probs[index]
  prior.M$ln.probs <- ln.probs - matrixStats::logSumExp(ln.probs)
  
  prior.M
}

#' Circos Plot of Pairwise Partition Probabilities from TIMBR Object
#'
#' Plots the probability that each pair of haplotypes is partitioned together on a circos plot; includes heatmap of relative effect sizes when using TIMBR results
#'
#' @param TIMBR.object results object from the TIMBR function; or prior.M object of "list" type 
#' @param file.path an optional file path for saving the plot as a PNG
#' @param plot.width PNG plot width
#' @param plot.height PNG plot height
#' @param post.summary posterior summary for relative effect size; select from c("mean", "median", "mode")
#' @param colors colors for heatmap of relative effect sizes
#' @param color.res number of colors to use for the heatmap color ramp
#'
#' @return circos plot of pairwise partition probabilities
#' 
#' @examples
#' #example data
#' data(mcv.data)
#' str(mcv.data)
#' 
#' #call TIMBR using CRP
#' results <- TIMBR(mcv.data$y, mcv.data$prior.D, mcv.data$prior.M$crp)
#' 
#' #plot partition probabilities
#' TIMBR.plot.circos(results)
#'
#' @export
TIMBR.plot.circos <- function(TIMBR.object, file.path=NULL, plot.width=480, plot.height=480,
                              post.summary="mean", colors=c("blue", "white", "red"), color.res=1000){
  #calculate pairwise probabilities for haplotype groupings
  if (is.null(TIMBR.object$y)){
    E.MMt <- lapply(1:length(TIMBR.object$ln.probs), function(x){M <- TIMBR:::M.matrix.from.ID(TIMBR.object$M.IDs[x]); exp(TIMBR.object$ln.probs[x])*tcrossprod(M)})
  } else {
    E.MMt <- lapply(1:length(TIMBR.object$p.M.given.y), function(x){M <- TIMBR:::M.matrix.from.ID(names(TIMBR.object$p.M.given.y)[x]); TIMBR.object$p.M.given.y[x]*tcrossprod(M)})
  }
  
  E.MMt <- Reduce("+", E.MMt)
  J <- ncol(E.MMt)
  rownames(E.MMt) <- LETTERS[1:J]
  colnames(E.MMt) <- rownames(E.MMt)
  
  #optimize order by minimizing distance of lines drawn on unit circle, weighed by probabilities
  orders <- combinat::permn(2:J, function(x){c(1,x)})
  locations <- cbind(sinpi(2/J*(0:(J-1))), cospi(2/J*(0:(J-1))))
  
  distance <- diag(0,J)
  distance[t(combn(1:J, 2))] <- apply(combn(1:8, 2), 2, function(x){sqrt(sum((locations[x[1],] - locations[x[2],])^2))})
  distance <- distance[upper.tri(distance)]
  
  order.dist <- sapply(orders, function(x){E.MMt.order <- E.MMt[x, x]; sum(E.MMt.order[upper.tri(E.MMt.order)]*distance)})
  best.order <- unlist(orders[which.min(order.dist)])
  E.MMt <- E.MMt[best.order, best.order]
  
  #plot connnections
  if (!is.null(file.path)){
    png(file.path, height=plot.height, width=plot.width)
  }
  
  circlize::circos.initialize(LETTERS[best.order], xlim=cbind(rep(0,J),rep(1,J)))
  circlize::circos.trackPlotRegion(y=rep(0,J), ylim=c(0,1))
  
  for (i in 1:(J-1)){
    for (j in (i+1):J){
      color <- rgb(0,0,0, max=255, alpha=E.MMt[i,j]*255)
      circlize::circos.link(rownames(E.MMt)[i], c(0.45,0.55), rownames(E.MMt)[j], c(0.45,0.55), col=color)
    }
  }
  
  #calculate signed squared effect divided by sum of squared effect and error variance
  if (is.null(TIMBR.object$y)){
    effects <- rep(0, J)
  } else {
    diff <- TIMBR.object$post.hap.effects - apply(TIMBR.object$post.hap.effects, 1, mean)
    sign <- (-1)^(diff < 0)
    sum.sq <- diff^2
    var.exp <- sum.sq / (sum.sq + TIMBR.object$post.sigma.sq)
    var.exp.signed <- var.exp*sign
    
    #posterior summary of transformed effects
    if (post.summary=="mean"){
      effects <- colMeans(var.exp.signed)
    } else if (post.summary=="mode"){
      effects <- apply(var.exp.signed, 2, function(x){dens <- density(x, from=-1, to=1); dens$x[which.max(dens$y)]})
    } else if (post.summary=="median"){
      effects <- apply(var.exp.signed,2,median)
    }
  }
  names(effects) <- LETTERS[1:J]
  
  #add effects and labels
  colors <- colorRampPalette(colors)(color.res+1)
  colors <- colors[round(color.res*(effects+1)/2)+1]
  
  for (i in 1:J){
    circlize::circos.rect(0, 0, 1, 1, LETTERS[i], col=colors[i])
    circlize::circos.text(0.5, 0.5, LETTERS[i], LETTERS[i], cex=1.5)
  }
  
  #clean-up and publish plot
  circlize::circos.clear()
  
  if (!is.null(file.path)){
    dev.off()
  }
}

#' @keywords internal
TIMBR.joint <- function(TIMBR.output.1, TIMBR.output.2){
  if (!identical(TIMBR.output.1$prior.M, TIMBR.output.2$prior.M)){
    stop("prior.M does not match")
  }
  
  ln.ml.h0 <- TIMBR.output.1$ln.ml + TIMBR.output.2$ln.ml
  
  approx.1 <- TIMBR.approx(TIMBR.output.1, ln.ml=T, return.prior=T)
  approx.2 <- TIMBR.approx(TIMBR.output.2, ln.ml=T, return.prior=T)
  
  data <- cbind(approx.1, approx.2[match(rownames(approx.1), rownames(approx.2)),])
  data <- data[!is.na(data[,3]),c(1,3,2)]
  
  ln.ml.hA <- matrixStats::logSumExp(rowSums(data))
  
  ln.BF <- ln.ml.hA - ln.ml.h0
  joint.p.M.given.Y <- rev(sort(exp(rowSums(data) - ln.ml.hA)))
  
  list(ln.BF=ln.BF, joint.p.M.given.Y=joint.p.M.given.Y)
}

#' @keywords internal
var.pop <- function(x){
  var(x)*(length(x) - 1)/length(x)
}

#' @keywords internal
scale.pop <- function(x){
  if (length(x)==1){
    0
  } else {
    x <- x - mean(x)
    x/sqrt(var.pop(x))
  }
}

#' @keywords internal
simulate.population <- function(M, N.J, var.exp, spacing="equal"){
  J <- nrow(M)
  K <- ncol(M)
  N <- J*N.J
  M.list <- apply(M, 1, match, x=1)
  
  #balanced homozygous population
  D.list <- rep(1:J, each=N/J)
  D <- matrix(0, N, J + choose(J,2))
  D[cbind(1:N, D.list)] <- 1
  
  #additive design matrix
  A <- additive.design(J, "happy")
  
  #sample standard normally- or equally-spaced effects, center at zero, set population variance equal to var.exp, sort by size
  if (spacing=="normal"){
    B <- rnorm(K)
  } else if (spacing=="equal"){
    B <- 1:K
  }
  B <- sort(scale.pop(B)*sqrt(var.exp))
  
  #sample standard normal error, center at zero, set population variance equal to 1-var.exp
  e <- scale.pop(rnorm(N))*sqrt(1-var.exp)
  
  #define phenotype
  y <- A[D.list,]%*%B[M.list] + e
  
  #return list of formatted input for TIMBR
  list(y=y, prior.D=list(P=D, A=A, fixed.diplo=T), B=B)
}

#' @keywords internal
log1mexp <- function(ln.p){
  sapply(ln.p, function(x){ifelse(x < -log(2), log1p(-exp(x)), log(-expm1(x)))})
}

#' @keywords internal
TIMBR2 <- function(y, prior.D, prior.M, prior.phi.b=1, samples=10000, samples.ml=10000, Z=NULL, W=NULL, verbose=T){
  
  TIMBR.sampler <- function(iterations, calc.null.ml=T, update.M=T, update.alpha=T){
    
    nglm.hyperparameters.ml <- function(MC, calc.partial.ml=T){
      #compute hyperparameters for the normal-gamma linear model
      if (update.M){
        d <- ncol(MC)
        n.params <- p + d
      }
      
      if (update.M | !fixed.diplo){
        X <- cbind(Z, DA%*%MC)
        XtWy <- crossprod(X,Wy)
        ZtWDAMC <- ZtWDA%*%MC
        CtMtAtDtWDAMC <- crossprod(MC, AtDtWDA)%*%MC
        V.star.inv <- rbind(cbind(ZtWZ, ZtWDAMC),cbind(t(ZtWDAMC),CtMtAtDtWDAMC))
      }
      
      if (d!=0){
        V.star.inv[cbind((p+1):n.params, (p+1):n.params)] <- CtMtAtDtWDAMC[cbind(1:d, 1:d)] + rep(phi.sq^(-1),d)
      }
      
      if (n.params > 0){
        L.Vt.inv <- backsolve(chol(V.star.inv), diag(n.params))
        V.star <- tcrossprod(L.Vt.inv)
      } else {
        L.Vt.inv <- V.star.inv
        V.star <- V.star.inv
      }
      
      m.star <- V.star%*%XtWy
      psi.star <- ytWy - c(crossprod(m.star, V.star.inv)%*%m.star)
      
      return.list <- list("kappa.star"=n, "psi.star"=psi.star, "m.star"=m.star, "V.star"=V.star, "L.Vt.inv"=L.Vt.inv)
      
      if (calc.partial.ml){
        partial.ln.ml <- -0.5*kappa.star*log(psi.star) - 0.5*d*log(phi.sq) + 0.5*determinant(V.star)$modulus[1]
        return.list$partial.ln.ml <- partial.ln.ml
      }
      
      return.list
    }
    
    sample.crp.concentration <- function(){
      #sample latent variable zeta
      ln.zeta <- log(rbeta(1, alpha+1, J))
      
      #sample alpha conditional on zeta
      z <- (prior.alpha.shape+K-1)/(J*(prior.alpha.rate-ln.zeta))
      pi <- z/(1+z)
      
      if (rbinom(1,1,pi)==1){
        alpha <- rgamma(1, prior.alpha.shape+K, prior.alpha.rate-ln.zeta)
      } else {
        alpha <- rgamma(1, prior.alpha.shape+K-1, prior.alpha.rate-ln.zeta)
      }
      
      alpha
    }
    
    #precompute matrix products if model and diplotypes are fixed
    if (!update.M & fixed.diplo){
      DAMC <- DA%*%MC
      X <- cbind(Z, DAMC)
      XtWy <- crossprod(X,Wy)
      ZtWDAMC <- ZtWDA%*%MC
      CtMtAtDtWDAMC <- crossprod(MC, AtDtWDA)%*%MC
      V.star.inv <- rbind(cbind(ZtWZ, ZtWDAMC),cbind(t(ZtWDAMC),CtMtAtDtWDAMC))
    }
    
    #create objects to store results
    post.M <- matrix(NA, iterations, J)
    post.MCbeta <- matrix(NA, iterations, J)
    post.delta <- matrix(NA, iterations, p)
    post.phi.sq <- rep(NA, iterations)
    post.sigma.sq <- rep(NA, iterations)
    post.alpha <- rep(NA, iterations)
    post.hyperparameters <- vector("list", iterations)
    post.K <- rep(NA, iterations)
    post.y.hat <- matrix(NA, iterations, n)
    p.D.given.y <- matrix(0, n, ncol.P)
    
    #iterate sampler
    for (i in 1:iterations){
      if (i%%1000==0 & verbose){
        print(i)
      }
      
      #compute matrix quantitites that depend on D
      if (!fixed.diplo){
        DA <- A[D.list,,drop=F]
        ZtWDA <- ZtW%*%DA
        AtDtWDA <- crossprod(DA*sqrt.W)
        AtDtWy <- crossprod(DA, Wy)
      }
      
      #sample allelic series matrix M
      if (update.M){
        if (model.type=="crp"){
          j.order <- 1:J
        } else {
          #randomize order due to non-exchangeable prior
          j.order <- sample(1:J)
        }
        
        #sample each row of M conditional on the other rows
        for (j in j.order){
          
          if (j!=j.order[1]){
            M.current <- list(M=M, M.list=M.list, M.posteriors=M.posteriors[[M.indicator]], new.index=M.list[j])
          } else {
            M.current <- list()
          }
          
          #set current row to zero and update matrix columns if necessary
          M[j,] <- 0
          
          if (!(1 %in% M[,M.list[j]])){
            M.update.index <- M.list[-j] > M.list[j]
            M.list[-j][M.update.index] <- M.list[-j][M.update.index] - 1
            M <- M[,-M.list[j],drop=F]
            M.current$new.index <- ncol(M)+1
          }
          
          M.list[j] <- NA
          K <- ncol(M)
          C <- contrast.list[[K]]
          
          #enumerate all possible assignments of current row of M
          M.list.space <- lapply(1:(K+1), function(x){M.list[j] <- x; M.list})
          MC.space <- lapply(1:K, function(x){C[M.list.space[[x]],,drop=F]})
          MC.space[[K+1]] <- contrast.list[[K+1]][M.list.space[[K+1]],,drop=F]
          
          #calculate prior for all possible assignments of current row of M
          if (model.type=="crp"){
            #analytic form for exchangeable prior
            colsums.M <- matrixStats::colSums2(M)
            M.ln.prior <- log(c(colsums.M, alpha))
          } else if (model.type=="uniform"){
            #constant non-exchangeable prior
            M.ln.prior <- rep(0, K+1)
          } else if (model.type=="list" | model.type=="mixture"){
            #arbitrary non-exchangeable prior or a mixture of that prior with the CRP
            #requires looking up priors via hash table using a unique naming scheme for M
            M.space.vec <- lapply(1:(K+1), function(x){M.list[j] <- x; M.list})
            
            #compute unique names for possible values of M
            if (hash.names){
              #use a hash table to map string IDs of M to their unique names
              M.space.name <- sapply(M.space.vec, paste, collapse=",")
              M.space.key <- lapply(M.space.name, function(x){prior.M.names[[x]]})
              M.space.key.null <- which(sapply(M.space.key, is.null))
              
              #compute unique names for new string IDs and update hash table
              if (length(M.space.key.null) != 0){
                M.space.key[M.space.key.null] <- sapply(M.space.vec[M.space.key.null], m.rename)
                list2env(setNames(M.space.key[M.space.key.null], M.space.name[M.space.key.null]), envir = prior.M.names)
              }
            } else {
              #compute unique names for possible settings of M
              M.space.key <- lapply(M.space.vec, m.rename)
            }
            
            #use the unique names to look up priors for possible settings of M using hash table
            M.ln.prior <- lapply(M.space.key, function(x){prior.M.hash[[x]]})
            M.ln.prior.null <- which(sapply(M.ln.prior, is.null))
            
            if (length(M.ln.prior.null)!=0){
              if (model.type=="list"){
                #values of M that are not in the hash table have probability zero
                M.ln.prior[M.ln.prior.null] <- -Inf
              } else {
                #compute mixture prior for new M
                missing.weighted.input <- sapply(M.space.key[M.ln.prior.null], function(x){prior.M.input[[x]]})
                missing.weighted.input <- sapply(missing.weighted.input, function(x){ifelse(is.null(x), -Inf, x + prior.M.weight.ln)})
                missing.weighted.crp <- sapply(M.space.vec[M.ln.prior.null], dcrp, prior.alpha=list(type="gamma", shape=prior.alpha.shape, rate=prior.alpha.rate)) + prior.M.weight.ln.1minus
                missing.mixture <- sapply(1:length(M.ln.prior.null), function(x){matrixStats::logSumExp(c(missing.weighted.input[x], missing.weighted.crp[x]))})
                
                #update hash table with mixture prior for new M
                M.ln.prior[M.ln.prior.null] <- missing.mixture
                list2env(setNames(M.ln.prior[M.ln.prior.null], M.space.key[M.ln.prior.null]), envir = prior.M.hash)
              }
            }
            M.ln.prior <- unlist(M.ln.prior)
          }
          
          #calculate t-distributed likelihood for all possible assignments of current row of M
          
          if (j==j.order[1]){
            #M.posteriors <- lapply(MC.space, nglm.hyperparameters.ml)
            M.posteriors <- lapply(1:(K+1), function(x){ifelse(M.ln.prior[x]==-Inf, 0, nglm.hyperparameters.ml(MC.space[[x]]))})
          } else {
            M.posteriors <- vector("list", K+1)
            M.posteriors[[M.current$new.index]] <- M.current$M.posteriors
            #M.posteriors[-M.current$new.index] <- lapply(MC.space[-M.current$new.index], nglm.hyperparameters.ml)
            M.posteriors[-M.current$new.index] <- lapply((1:(K+1))[-M.current$new.index], function(x){ifelse(M.ln.prior[x]==-Inf, 0, nglm.hyperparameters.ml(MC.space[[x]]))})
          }
          
          M.ln.ml <- unlist(lapply(M.posteriors, function(x){x$partial.ln.ml}))
          
          #combine likelihoods with priors and scale by normalizing constant
          M.prob <- M.ln.ml + M.ln.prior
          M.prob <- exp(M.prob - matrixStats::logSumExp(M.prob))
          
          #sample assignment for current row of M from categorical distribution
          M.indicator <- match(rmultinom(1,1,M.prob), x=1)
          
          if (isTRUE(M.indicator==M.current$new.index) & isTRUE(M.current$M.list[j]!=M.current$new.index)){
            M.list <- M.current$M.list
            M <- M.current$M
            K <- K + 1
          } else {
            M.list[j] <- M.indicator
            
            #update M
            if (M.indicator > ncol(M)){
              M <- cbind(M,0)
              M[j, M.indicator] <- 1
              K <- K + 1
            } else {
              M[j, M.indicator] <- 1
            }
          }
        }
        #update quantities that depend on M
        C <- contrast.list[[K]]
        MC <- M%*%C
        AMC <- A%*%MC
        d <- ncol(C)
      }
      
      #sample concentration parameter if using CRP
      if (update.alpha){
        alpha <- sample.crp.concentration()
        if (prior.alpha.type=="beta.prime"){
          prior.alpha.rate <- rgamma(1, prior.alpha.b + prior.alpha.shape, prior.alpha.q + alpha)
        }
      }
      
      #sample error variance and linear coefficients from conjugate normal-gamma distribution 
      #note that coefficients are scaled as eta=beta/lambda
      #store hyperparameters for normal-gamma distribution 
      if (update.M){
        M.posteriors <- M.posteriors[[M.indicator]]
      } else {
        M.posteriors <- nglm.hyperparameters.ml(MC)
      }
      
      #sample error variance from inv-gamma distribution
      sigma.sq.inv <- rgamma(1, 0.5*M.posteriors$kappa.star, 0.5*M.posteriors$psi.star)
      sigma.sq <- sigma.sq.inv^(-1)
      sigma <- sqrt(sigma.sq)
      
      #sample linear coefficients from normal distribution and account for scaling
      theta <- c(M.posteriors$m.star + sigma*(M.posteriors$L.Vt.inv%*%rnorm(p+d)))
      
      if (p > 0){
        theta[-(1:p)] <- theta[-(1:p)]/lambda
        delta <- theta[1:p]
        eta <- theta[-(1:p)]
      } else {
        theta <- theta/lambda
        delta <- rep(0,p)
        eta <- theta
      }
      
      #update tau.sq from conjugate inv-gamma distribution
      tau.sq.shape <- prior.phi.b + 0.5*d
      tau.sq.rate <- 0.5 + 0.5*sigma.sq.inv*sum(eta^(2))
      tau.sq <- rgamma(1, tau.sq.shape, tau.sq.rate)^(-1)
      
      #update lambda from conjugate normal distribution
      Z.delta <- Z%*%delta
      y.prime <- c(y-Z.delta)
      MCeta <- MC%*%eta
      
      if (update.M | !fixed.diplo){
        lambda.var <- (sigma.sq.inv*c(crossprod(MCeta, AtDtWDA)%*%MCeta) + 1)^(-1)
        lambda.mean <- sum(sigma.sq.inv*y.prime*W*lambda.var*c(DA%*%MCeta))
      } else {
        lambda.var <- (sigma.sq.inv*c(crossprod(eta, CtMtAtDtWDAMC)%*%eta) + 1)^(-1)
        lambda.mean <- sum(sigma.sq.inv*y.prime*W*lambda.var*c(DAMC%*%eta))
      }
      
      lambda <- rnorm(1, lambda.mean, sqrt(lambda.var))
      
      #update quantities that depend on tau.sq and lambda
      phi.sq <- lambda^2*tau.sq
      beta <- lambda*eta
      MCbeta <- lambda*(MCeta)
      AMCbeta <- A%*%MCbeta
      
      #update dipltypes jointly from independent categorical distributions
      if (!fixed.diplo){
        #calculate independent normal likelihood for all possible assignments of each row of D
        D.ln.ml <- matrix(dnorm(rep(y.prime, each=ncol.P), rep(AMCbeta, n), rep(sigma*sqrt.W^(-1), each=ncol.P), log=T), n, ncol.P, byrow=T)
        
        #combine likelihood with prior and scale to prevent underflow of all probabilities
        D.prob <- D.ln.ml + ln.P
        D.prob <- exp(D.prob - matrixStats::rowMaxs(D.prob))
        
        #sample assignment for each row of D from independent categorical distributions
        #probabilities are normalized by rmultinom
        D <- t(apply(D.prob, 1, rmultinom, n=1, size=1))
        D.list <- apply(D, 1, match, x=1)
      }
      
      #store posterior samples and hyperparameters
      post.M[i,] <- M.list
      post.MCbeta[i,] <- MCbeta
      post.delta[i,] <- delta
      post.phi.sq[i] <- phi.sq
      post.sigma.sq[i] <- sigma.sq
      post.alpha[i] <- alpha
      post.hyperparameters[[i]] <- M.posteriors[1:4]
      post.K[i] <- K
      post.y.hat[i,] <- AMCbeta[D.list] + Z.delta
      p.D.given.y <- p.D.given.y + D
    }
    
    #report unique names for posterior samples of M
    if (update.M==F){
      post.M <- rep(m.rename(post.M[1,]), iterations)
    } else if (hash.names){
      post.M <- apply(post.M, 1, function(x){prior.M.names[[paste(x, collapse=",")]]})
    } else {
      post.M <- apply(post.M, 1, m.rename)
    }
    
    #calculate marginal posterior diplotype probabilities
    p.D.given.y <- p.D.given.y/iterations
    
    #update variable state for potential reduced run of the sampler
    D <<- D
    D.list <<- D.list
    lambda <<- lambda
    tau.sq <<- tau.sq
    phi.sq <<- phi.sq
    
    #return posterior samples and hyperparameters
    posterior.results <- list("post.M"=post.M, "post.MCbeta"=post.MCbeta, "post.delta"=post.delta, "post.sigma.sq"=post.sigma.sq, "post.phi.sq"=post.phi.sq, 
                              "p.D.given.y"=p.D.given.y, "post.hyperparameters"=post.hyperparameters, "post.K"=post.K, "post.y.hat"=post.y.hat)
    
    if (update.alpha){
      posterior.results$post.alpha <- post.alpha
    }
    
    if (calc.null.ml){
      update.M <- T
      posterior.results$ln.ml.null <- nglm.hyperparameters.ml(matrix(1,J,1)%*%sumtozero.contrast(1), calc.partial.ml=T)$partial.ln.ml
    }
    
    posterior.results
  }
  
  ####################
  #precompute invariant hyperparameters and matrix products
  #specify starting values for sampler
  n <- length(y)
  
  #default covariate matrix includes a mean
  if (is.null(Z)){
    Z <- matrix(1,n,1)
  }
  
  #default covariance matrix assumes no replicates
  if (is.null(W)){
    W <- rep(1,n)
  }
  
  model.type <- prior.M$model.type
  fixed.diplo <- prior.D$fixed.diplo
  
  P <- prior.D$P
  A <- prior.D$A
  
  p <- ncol(Z)
  J <- ncol(A)
  ncol.P <- ncol(P)
  sqrt.W <- sqrt(W)
  ZtWZ <- crossprod(Z*sqrt.W)
  ZtW <- t(Z*W)
  Wy <- W*y
  ytWy <- sum(y*Wy)
  kappa.star <- n
  
  phi.sq <- 0.5/prior.phi.b
  lambda <- 1
  tau.sq <- phi.sq/(lambda^2)
  
  D <- matrix(0, n, ncol.P)
  D.list <- apply(P, 1, which.max)
  D[cbind(1:n, D.list)] <- 1
  
  if (fixed.diplo){
    DA <- A[D.list,,drop=F]
    ZtWDA <- ZtW%*%DA
    AtDtWDA <- crossprod(DA*sqrt.W)
    AtDtWy <- crossprod(DA, Wy)
    P <- D
  } else {
    ln.P <- log(P)
  }
  
  hash.names <- F
  
  if (model.type=="fixed"){
    alpha <- NA
    M <- M.matrix.from.ID(prior.M$M.IDs)
    K <- ncol(M)
    C <- sumtozero.contrast(K)
    MC <- M%*%C
    AMC <- A%*%MC
    d <- ncol(C)
    n.params <- p + d
  } else {
    contrast.list <- lapply(1:J, sumtozero.contrast)
  }
  
  if (model.type=="crp"){
    prior.alpha.type <- prior.M$prior.alpha.type
    if (prior.alpha.type=="gamma"){
      prior.alpha.shape <- prior.M$prior.alpha.shape
      prior.alpha.rate <- prior.M$prior.alpha.rate
      alpha <- prior.alpha.shape/prior.alpha.rate
    } else if (prior.alpha.type=="fixed"){
      alpha <- prior.M$prior.alpha
    } else if (prior.alpha.type=="beta.prime"){
      prior.alpha.shape <- prior.M$prior.alpha.a
      prior.alpha.b <- prior.M$prior.alpha.b
      prior.alpha.q <- ifelse(is.null(prior.M$prior.alpha.q), 1, prior.M$prior.alpha.q)
      prior.alpha.rate <- prior.alpha.b/prior.alpha.q
      alpha <- prior.alpha.shape/prior.alpha.rate
    }
    M <- matrix(1, J, 1)
  } else if (model.type=="uniform"){
    alpha <- NA
    M <- matrix(1, J, 1)
  } else if (model.type=="list" | model.type=="mixture"){
    prior.M.hash <- new.env(hash = T)
    alpha <- NA
    M <- M.matrix.from.ID(prior.M$M.IDs[which.max(prior.M$ln.probs)])
    
    if (model.type=="list"){
      list2env(setNames(as.list(prior.M$ln.probs), prior.M$M.IDs), envir = prior.M.hash)
    } else {
      prior.M.input <- new.env(hash = T)
      list2env(setNames(as.list(prior.M$ln.probs), prior.M$M.IDs), envir = prior.M.input)
      prior.alpha.shape <- prior.M$prior.alpha.shape
      prior.alpha.rate <- prior.M$prior.alpha.rate
      prior.M.weight.ln <- log(prior.M$weight)
      prior.M.weight.ln.1minus <- log(1-prior.M$weight)
    }
    
    hash.names <- prior.M$hash.names
    if (hash.names){
      prior.M.names <- new.env(hash = T)
    }
  }
  
  M.list <- apply(M, 1, match, x=1)
  
  ####################
  #iterate TIMBR sampler
  if (verbose){
    print("Sampling from the full posterior")
  }
  
  if (model.type=="crp"){
    if (prior.alpha.type=="gamma" | prior.alpha.type=="beta.prime"){
      results <- TIMBR.sampler(samples)
    } else if (prior.alpha.type=="fixed"){
      results <- TIMBR.sampler(samples, update.alpha=F)
    }
  } else if (model.type=="fixed"){
    results <- TIMBR.sampler(samples, update.M=F, update.alpha=F)
  } else if (model.type=="uniform" | model.type=="list" | model.type=="mixture"){
    results <- TIMBR.sampler(samples, update.alpha=F)
  }
  
  #collect posterior probability of the null model
  ln.ml.null <- results$ln.ml.null
  
  post.M.ranked <- -sort(-table(results$post.M))/samples
  post.M.null <- post.M.ranked[names(post.M.ranked)==paste(rep(0, J), collapse=",")]
  
  if(length(post.M.null)==0){
    post.M.null <- 0
  }
  
  #if posterior probability of null model is relatively high, use this to calculate marginal likelihood
  if (names(post.M.ranked[1])==paste(rep(0, J), collapse=",") | post.M.null>=0.01){
    if (model.type=="crp"){
      if (prior.alpha.type=="gamma"){
        ln.ml <- ln.ml.null + dcrp(rep(0,J), list(type="gamma", shape=prior.alpha.shape, rate=prior.alpha.rate)) - log(post.M.null)
      } else if (prior.alpha.type=="fixed"){
        ln.ml <- ln.ml.null + dcrp(rep(0,J), list(type="fixed", alpha=prior.M$prior.alpha)) - log(post.M.null)
      } else if (prior.alpha.type=="beta.prime"){
        ln.ml <- ln.ml.null + dcrp(rep(0,J), list(type="beta.prime", a=prior.alpha.shape, b=prior.alpha.b, q=prior.alpha.q)) - log(post.M.null)
      }
    } else if (model.type=="fixed"){
      ln.ml <- ln.ml.null
    } else if (model.type=="uniform"){
      ln.ml <- ln.ml.null - ln.bell(J) - log(post.M.null)
    } else if (model.type=="list" | model.type=="mixture"){
      ln.ml <- ln.ml.null + prior.M.hash[[paste(rep(0,J), collapse=",")]] - log(post.M.null)
    }
  } else {
    #obtain samples of the normal-gamma hyperparameters conditional on a single model matrix M
    if (model.type=="fixed"){
      #these are given when the model.type is fixed
      nglm.hyperparameters <- results$post.hyperparameters
      post.phi.sq <- results$post.phi.sq
      p.D.given.y <- results$p.D.given.y
      samples.ml <- samples
    } else {
      #set M to the MAP and update dependent quantities
      M <- M.matrix.from.ID(names(post.M.ranked)[1])
      M.list <- apply(M, 1, match, x=1)
      K <- ncol(M)
      C <- contrast.list[[K]]
      MC <- M%*%C
      AMC <- A%*%MC
      d <- ncol(C)
      n.params <- p + d
      
      #iterate reduced sampler run with M fixed at MAP
      if (verbose){
        print("Sampling from the conditional posterior")
      }
      
      reduced.results <- TIMBR.sampler(samples.ml, update.M=F, calc.null.ml=F, update.alpha=F)
      nglm.hyperparameters <- reduced.results$post.hyperparameters
      post.phi.sq <- reduced.results$post.phi.sq
      p.D.given.y <- reduced.results$p.D.given.y
    }
    
    #set remaining variables to values with high posterior probability given M fixed at MAP
    kappa.star <- n
    psi.star <- sapply(nglm.hyperparameters, function(x){x$psi.star})
    m.star <- lapply(nglm.hyperparameters, function(x){x$m.star})
    V.star <- lapply(nglm.hyperparameters, function(x){x$V.star})
    
    theta <- matrixStats::rowMeans2(sapply(1:length(m.star), function(x){sapply(unlist(m.star[x]), function(x){x})}))
    
    if (p > 0){
      delta <- theta[1:p]
      beta <- theta[-(1:p)]
    } else {
      delta <- rep(0,p)
      beta <- theta
    }
    
    sigma.sq <- mean((kappa.star/psi.star)^(-1))
    
    if (!fixed.diplo){
      D <- matrix(0, n, ncol.P)
      D.list <- apply(p.D.given.y, 1, which.max)
      D[cbind(1:n, D.list)] <- 1
    }
    
    #calculate partial marginal likelihood at point of high posterior probability
    Zdelta <- Z%*%delta
    AMCbeta <- A%*%MC%*%beta
    
    p1 <- sum(dnorm(y, AMCbeta[D.list] + Zdelta, sqrt(sigma.sq*W^(-1)), log=T))
    p4 <- log(sigma.sq)
    p5 <- ln.beta.prior.marginalized(beta, sigma.sq, prior.phi.b)
    p7 <- matrixStats::logSumExp(dgamma(sigma.sq^(-1), 0.5*kappa.star, 0.5*psi.star, log=T))-log(samples.ml)
    p8 <- matrixStats::logSumExp(sapply(1:samples.ml, function(x){mvtnorm::dmvnorm(theta, m.star[[x]], V.star[[x]]*sigma.sq, log=T)}))-log(samples.ml)
    c <- -0.5*n*log(pi) + lgamma(0.5*kappa.star) - 0.5*sum(-log(W))
    
    if (!fixed.diplo){
      p2 <- sum(ln.P[cbind(1:n, D.list)])
      
      D.ln.ml <- matrix(dnorm(rep((y - Zdelta), each=ncol.P), rep(AMCbeta, n), rep(sqrt(sigma.sq)*sqrt.W^(-1), each=ncol.P), log=T), n, ncol.P, byrow=T)
      D.ln.prob <- D.ln.ml + ln.P
      D.ln.prob <-  D.ln.prob - matrixStats::rowLogSumExps(D.ln.prob)
      p9 <- sum(D.ln.prob[cbind(1:n, D.list)])
    } else {
      p2 <- 0
      p9 <- 0
    }
    
    if (model.type=="crp"){
      if (prior.alpha.type=="gamma"){
        p3 <- dcrp(apply(M, 1, match, x=1), list(type="gamma", shape=prior.alpha.shape, rate=prior.alpha.rate))
      } else if (prior.alpha.type=="fixed"){
        p3 <- dcrp(apply(M, 1, match, x=1), list(type="fixed", alpha=prior.M$prior.alpha))
      } else if (prior.alpha.type=="beta.prime"){
        p3 <- dcrp(apply(M, 1, match, x=1), list(type="beta.prime", a=prior.alpha.shape, b=prior.alpha.b, q=prior.alpha.q))
      }
      p6 <- log(post.M.ranked[1])
    } else if (model.type=="fixed"){
      p3 <- 0
      p6 <- 0
    } else if (model.type=="uniform"){
      p3 <- -ln.bell(J)
      p6 <- log(post.M.ranked[1])
    } else if (model.type=="list" | model.type=="mixture"){
      p3 <- prior.M.hash[[names(post.M.ranked[1])]]
      p6 <- log(post.M.ranked[1])
    }
    
    names(p6) <- NULL
    ln.ml <- p1 + p2 + p3 + p4 + p5 - p6 - p7 - p8 - p9 - c
  }
  
  #return posterior samples, marginal densities, and marginal likelihood
  output <- list(y=y, prior.D=prior.D, prior.M=prior.M, prior.phi.b=prior.phi.b, samples=samples, samples.ml=samples.ml, Z=Z, W=W)
  output <- c(output, results[names(results) != "post.hyperparameters"], ln.ml=ln.ml, ln.BF=ln.ml-ln.ml.null, p.M.given.y=list(post.M.ranked), 
              post.var.exp=list(results$post.phi.sq/(results$post.phi.sq + 1)), post.hap.effects=list(results$post.MCbeta+results$post.delta[,1]),
              p.K.given.y=list(table(results$post.K)/samples))
  output <- output[order(names(output))]
}
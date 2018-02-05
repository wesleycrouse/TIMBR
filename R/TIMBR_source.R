#' @keywords internal
sumtozero.contrast <- function(k){
  u <- 1/((k-1)^(1/2))
  v <- (-1+(k^(1/2)))*((k-1)^(-3/2))
  w <- (k-2)*v+u
  
  C <- matrix(-v,k-1,k-1)
  diag(C) <- rep(w,k-1)
  rbind(C, rep(-u,k-1))
}

#' @keywords internal
model.rename <- function(m){
  unique.m <- unique(m)
  sapply(m, function(y){match(x=y, unique.m)-1})
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
model.matrix.from.ID <- function(M.ID){
  m <- as.numeric(unlist(strsplit(M.ID, ","))) + 1
  J <- length(m)
  M <- matrix(0, J, max(m))
  M[cbind(1:J, m)] <- 1
  M
}

#' @keywords internal
beta.prior.marginalized <- function(beta, sigma.sq, prior.v.b, prior.v.a=0.5){
  d <- length(beta)
  
  #parameters for hypergeometric U
  z.U <- 0.5*sum(beta^2)/sigma.sq
  a.U <- prior.v.b + 0.5*d
  b.U <- prior.v.a + 0.5*d
  
  (-0.5*d)*log(2) - 0.5*d*log(pi) - 0.5*d*log(sigma.sq) - lbeta(0.5, prior.v.b) + lgamma(a.U) + log(gsl::hyperg_U(a.U, b.U, z.U))
}

#' @keywords internal
model.prior.crp.marginalized <- function(m, prior.alpha.a, prior.alpha.b, scale.by.max=F){
  J.k <- table(m, dnn=NULL)
  K <- length(J.k)
  J <- length(m)
  sum.lgamma.J.k <- sum(lgamma(J.k))
  
  C <- sum.lgamma.J.k + prior.alpha.a*log(prior.alpha.b) - lgamma(prior.alpha.a)
  
  density.crp.concentration <- Vectorize(function(x, log=F, scale.by = 0, trans.domain = F){
    if (trans.domain){x <- x/(1-x)}
    density <- lgamma(x) - lgamma(x+J) + (prior.alpha.a+K-1)*log(x) - prior.alpha.b*x - scale.by
    ifelse(log, density, exp(density))
  }, vectorize.args = "x")
  
  if (scale.by.max){
    max.ln.density <- optimize(density.crp.concentration, c(0,1), log=T, trans.domain=T, maximum=T)$objective
    integral <- integrate(density.crp.concentration, lower=0, upper=Inf, scale.by=max.ln.density)
    log(integral$value) + C + max.ln.density
  } else {
    integral <- integrate(density.crp.concentration, lower=0, upper=Inf)
    log(integral$value) + C
  }
}

#' Tree-based Inference of Multiallelism via Bayesian Regression
#'
#' Posterior samples and Bayes Factors using the TIMBR model
#'
#' @param y vector of phenotype values for each strain
#' @param prior.D list of inputs for the prior distribution of strain diplotype states; see data(mcv.data) for an example
#' @param prior.M list of inputs for the prior distribution of the allelic series model; see data(mcv.data) for examples
#' @param prior.v.b shape parameter for the beta prime prior distribution on the variance component
#' @param sample number of samples to draw from the full posterior
#' @param samples.ml number of samples to draw from the condtiional posterior (if necessary)
#' @param Z design matrix for covariates; first column must be a vector of ones, which is the default
#' @param W vector of replicates for each strain; one replicate per strain by default
#' @param verbose optionally report function progress
#'
#' @return returns a list of input parameters, posterior samples, and the marginal likelihood
#' 
#' @examples
#' #example data
#' data(mcv.data)
#' str(mcv.data)
#' 
#' #call TIMBR using CRP
#' results <- TIMBR(mcv.data$y, mcv.data$prior.D, mcv.data$prior.M$crp)
#' 
#' #calculate the Bayes Factor
#' results$ln.ml - results$ln.ml.null
#' 
#' #report posterior probabilities for the top allelic series models
#' head(rev(sort(table(results$post.M))))/results$samples
#'
#' @export
TIMBR <- function(y, prior.D, prior.M, prior.v.b=1, samples=1000, samples.ml=1000, Z=NULL, W=NULL, verbose=T){
  
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
      #sample latent variable eta
      zeta <- rbeta(1, alpha+1, J)
      ln.zeta <- log(zeta)
      
      #sample alpha conditional on eta
      z <- (prior.alpha.a+k-1)/(J*(prior.alpha.b-ln.zeta))
      pi <- z/(1+z)
      
      if (rbinom(1,1,pi)==1){
        alpha <- rgamma(1, prior.alpha.a+k, (prior.alpha.b-ln.zeta))
      } else {
        alpha <- rgamma(1, prior.alpha.a+k-1, (prior.alpha.b-ln.zeta))
      }
      
      alpha
    }
    
    #precompute matrix products if model and diplotypes are fixed
    if (!update.M & fixed.diplo){
      X <- cbind(Z, DA%*%MC)
      XtWy <- crossprod(X,Wy)
      ZtWDAMC <- ZtWDA%*%MC
      CtMtAtDtWDAMC <- crossprod(MC, AtDtWDA)%*%MC
      V.star.inv <- rbind(cbind(ZtWZ, ZtWDAMC),cbind(t(ZtWDAMC),CtMtAtDtWDAMC))
      DAMC <<- DA%*%MC
    }
    
    #create objects to store results
    post.M <- matrix(NA, iterations, J)
    post.D <- matrix(0, n, ncol.P)
    post.MCbeta <- matrix(NA, iterations, J)
    post.delta <- matrix(NA, iterations, p)
    post.phi.sq <- rep(NA, iterations)
    post.sigma.sq <- rep(NA, iterations)
    post.alpha <- rep(NA, iterations)
    post.hyperparameters <- vector("list", iterations) 
    
    #iterate sampler
    for (i in 1:iterations){
      if (i%%1000==0 & verbose){
        print(i)
      }
      
      #compute matrix quantitites that depend on D
      if (!fixed.diplo){
        DA <- D%*%A
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
          #set current row to zero and update matrix columns if necessary
          M[j,] <- 0
          
          if (!(1 %in% M[,M.list[j]])){
            M.update.index <- M.list[-j] > M.list[j]
            M.list[-j][M.update.index] <- M.list[-j][M.update.index] - 1
            M <- M[,-M.list[j],drop=F]
          }
          
          M.list[j] <- NA
          k <- ncol(M)
          C <- contrast.list[[k]]
          
          #calculate t-distributed likelihood for all possible assignments of current row of M
          MC.space <- lapply(1:k,function(x){M[j,x]<-1; M%*%C})
          MC.space[[k+1]] <- cbind(M,c(rep(0,j-1),1,rep(0,J-j)))%*%contrast.list[[k+1]]
          M.posteriors <- lapply(MC.space, nglm.hyperparameters.ml)
          M.ln.ml <- unlist(lapply(M.posteriors, function(x){x$partial.ln.ml}))
          
          #calculate prior for all possible assignments of current row of M
          if (model.type=="crp"){
            #analytic form for exchangeable prior
            colsums.M <- matrixStats::colSums2(M)
            M.ln.prior <- log(c(colsums.M, alpha)) - log(sum(colsums.M + alpha))
          } else if (model.type=="uniform"){
            #constant non-exchangeable prior
            M.ln.prior <- rep(0, k+1)
          } else if (model.type=="list" | model.type=="mixture"){
            #arbitrary non-exchangeable prior or a mixture of that prior with the CRP
            #requires looking up priors via hash table using a unique naming scheme for M
            M.space.vec <- lapply(1:(k+1), function(x){M.list[j] <- x; M.list})
            
            #compute unique names for possible values of M
            if (hash.names){
              #use a hash table to map string IDs of M to their unique names
              M.space.name <- sapply(M.space.vec, paste, collapse=",")
              M.space.key <- lapply(M.space.name, function(x){prior.M.names[[x]]})
              M.space.key.null <- which(sapply(M.space.key, is.null))
              
              #compute unique names for new string IDs and update hash table
              if (length(M.space.key.null) != 0){
                M.space.key[M.space.key.null] <- sapply(M.space.vec[M.space.key.null], function(x){paste(model.rename(x), collapse=",")})
                list2env(setNames(M.space.key[M.space.key.null], M.space.name[M.space.key.null]), envir = prior.M.names)
              }
            } else {
              #compute unique names for possible settings of M
              M.space.key <- lapply(M.space.vec, function(x){paste(model.rename(x), collapse=",")})
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
                missing.weighted.crp <- sapply(M.space.vec[M.ln.prior.null], model.prior.crp.marginalized, prior.alpha.a=prior.alpha.a, prior.alpha.b=prior.alpha.b) + prior.M.weight.ln.1minus
                missing.mixture <- sapply(1:length(M.ln.prior.null), function(x){matrixStats::logSumExp(c(missing.weighted.input[x], missing.weighted.crp[x]))})
                
                #update hash table with mixture prior for new M
                M.ln.prior[M.ln.prior.null] <- missing.mixture
                list2env(setNames(M.ln.prior[M.ln.prior.null], M.space.key[M.ln.prior.null]), envir = prior.M.hash)
              }
            }
            M.ln.prior <- unlist(M.ln.prior)
          }
          
          #combine likelihoods with priors and scale by normalizing constant
          M.prob <- M.ln.ml + M.ln.prior
          M.prob <- exp(M.prob - matrixStats::logSumExp(M.prob))
          
          #sample assignment for current row of M from categorical distribution
          M.indicator <- match(rmultinom(1,1,M.prob), x=1)
          M.list[j] <- M.indicator
          
          #update M
          if (M.indicator > ncol(M)){
            M <- cbind(M,0)
            M[j, M.indicator] <- 1
            k <- k + 1
          } else {
            M[j, M.indicator] <- 1
          }
        }
        #update quantities that depend on M
        C <- contrast.list[[k]]
        MC <- M%*%C
        AMC <- A%*%MC
        d <- ncol(C)
      }
      
      #sample concentration parameter if using CRP
      if (update.alpha){
        alpha <- sample.crp.concentration()
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
      tau.sq.shape <- prior.v.b + 0.5*d
      tau.sq.rate <- 0.5 + 0.5*sigma.sq.inv*sum(eta^(2))
      tau.sq <- rgamma(1, tau.sq.shape, tau.sq.rate)^(-1)
      
      #update lambda from conjugate normal distribution
      y.prime <- c(y-Z%*%delta)
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
      
      #update dipltypes jointly from independent categorical distributions
      if (!fixed.diplo){
        AMCbeta <- AMC%*%beta
        
        #calculate independent normal likelihood for all possible assignments of each row of D
        D.ln.ml <- matrix(dnorm(rep(y.prime, each=ncol.P), rep(AMCbeta, n), rep(sigma*sqrt.W^(-1), each=ncol.P), log=T), n, ncol.P, byrow=T)
        
        #combine likelihood with prior and scale to prevent underflow of all probabilities
        D.prob <- D.ln.ml + ln.P
        D.prob <- exp(D.prob - matrixStats::rowMaxs(D.prob))
        
        #sample assignment for each row of D from independent categorical distributions
        #probabilities are normalized by rmultinom
        D <- t(apply(D.prob, 1, rmultinom, n=1, size=1))
      }
      
      #store posterior samples and hyperparameters
      post.M[i,] <- M.list
      post.MCbeta[i,] <- MCbeta
      post.delta[i,] <- delta
      post.phi.sq[i] <- phi.sq
      post.sigma.sq[i] <- sigma.sq
      post.alpha[i] <- alpha
      post.hyperparameters[[i]] <- M.posteriors[1:4]
      
      if (!fixed.diplo){
        post.D <- post.D + D.prob
      } else {
        post.D <- post.D + D
      }
    }
    
    #report unique names for posterior samples of M
    if (update.M==F){
      post.M <- rep(paste(model.rename(post.M[1,]), collapse=","), iterations)
    } else if (hash.names){
      post.M <- apply(post.M, 1, function(x){prior.M.names[[paste(x, collapse=",")]]})
    } else {
      post.M <- apply(post.M, 1, function(x){paste(model.rename(x), collapse=",")})
    }
    
    #calculate marginal posterior diplotype probabilities
    post.D <- post.D/iterations
    
    #update variable state for potential reduced run of the sampler
    D <<- D
    lambda <<- lambda
    tau.sq <<- tau.sq
    phi.sq <<- phi.sq
    
    #return posterior samples and hyperparameters
    posterior.results <- list("post.D"=post.D, "post.M"=post.M, "post.MCbeta"=post.MCbeta, "post.delta"=post.delta, "post.sigma.sq"=post.sigma.sq, "post.phi.sq"=post.phi.sq, "post.hyperparameters"=post.hyperparameters)
    
    if (update.alpha){
      posterior.results$post.alpha <-  post.alpha
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
  
  phi.sq <- 0.5/prior.v.b
  lambda <- 1
  tau.sq <- phi.sq/(lambda^2)
  
  D <- matrix(0, n, ncol.P)
  D[cbind(1:n, apply(P, 1, which.max))] <- 1
  
  if (fixed.diplo){
    DA <- D%*%A
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
    M <- model.matrix.from.ID(prior.M$M.IDs)
    C <- sumtozero.contrast(ncol(M))
    MC <- M%*%C
    AMC <- A%*%MC
    d <- ncol(C)
    n.params <- p + d
  } else {
    contrast.list <- lapply(1:J, sumtozero.contrast)
  }
  
  if (model.type=="crp"){
    prior.alpha.a <- prior.M$prior.alpha.a
    prior.alpha.b <- prior.M$prior.alpha.b
    alpha <- prior.alpha.a/prior.alpha.b
    M <- matrix(1, J, 1)
  } else if (model.type=="uniform"){
    alpha <- NA
    M <- matrix(1, J, 1)
  } else if (model.type=="list" | model.type=="mixture"){
    prior.M.hash <- new.env(hash = T)
    alpha <- NA
    M <- model.matrix.from.ID(prior.M$M.IDs[which.max(prior.M$probs)])
    
    if (model.type=="list"){
      list2env(setNames(as.list(log(prior.M$probs)), prior.M$M.IDs), envir = prior.M.hash)
    } else {
      prior.M.input <- new.env(hash = T)
      list2env(setNames(as.list(log(prior.M$probs)), prior.M$M.IDs), envir = prior.M.input)
      prior.alpha.a <- prior.M$prior.alpha.a
      prior.alpha.b <- prior.M$prior.alpha.b
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
    results <- TIMBR.sampler(samples)
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
  if (post.M.null>=0.01){
    if (model.type=="crp"){
      ln.ml <- ln.ml.null + model.prior.crp.marginalized(rep(0,J), prior.alpha.a, prior.alpha.b) - log(post.M.null)
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
      post.D <- results$post.D
      samples.ml <- samples
    } else {
      #set M to the MAP and update dependent quantities
      M <- model.matrix.from.ID(names(post.M.ranked)[1])
      M.list <- apply(M, 1, match, x=1)
      C <- contrast.list[[ncol(M)]]
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
      post.D <- reduced.results$post.D
    }
    
    #set remaining variables to values with high posterior probability given M fixed at MAP
    kappa.star <- n
    psi.star <- sapply(nglm.hyperparameters, function(x){x$psi.star})
    m.star <- lapply(nglm.hyperparameters, function(x){x$m.star})
    V.star <- lapply(nglm.hyperparameters, function(x){x$V.star})
    
    theta <- rowMeans(sapply(1:length(m.star), function(x){sapply(unlist(m.star[x]), function(x){x})}))
    
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
      D[cbind(1:n, apply(post.D, 1, which.max))] <- 1
    }
    
    #calculate partial marginal likelihood at point of high posterior probability
    Zdelta <- Z%*%delta
    AMCbeta <- A%*%MC%*%beta
    
    p1 <- sum(dnorm(y, Zdelta + D%*%AMCbeta, sqrt(sigma.sq*W^(-1)), log=T))
    p4 <- log(sigma.sq)
    p5 <- beta.prior.marginalized(beta, sigma.sq, prior.v.b)
    p7 <- matrixStats::logSumExp(dgamma(sigma.sq^(-1), 0.5*kappa.star, 0.5*psi.star, log=T))-log(samples.ml)
    p8 <- matrixStats::logSumExp(sapply(1:samples.ml, function(x){mvtnorm::dmvnorm(theta, m.star[[x]], V.star[[x]]*sigma.sq, log=T)}))-log(samples.ml)
    c <- -0.5*n*log(pi) + lgamma(0.5*kappa.star) - 0.5*sum(-log(W))
    
    if (!fixed.diplo){
      D.states <- apply(D,1,match,x=1)
      p2 <- sum(ln.P[cbind(1:n, D.states)])
      
      D.ln.ml <- matrix(dnorm(rep((y - Zdelta), each=ncol.P), rep(AMCbeta, n), rep(sqrt(sigma.sq)*sqrt.W^(-1), each=ncol.P), log=T), n, ncol.P, byrow=T)
      D.ln.prob <- D.ln.ml + ln.P
      D.ln.prob <-  D.ln.prob - matrixStats::rowLogSumExps(D.ln.prob)
      p9 <- sum(D.ln.prob[cbind(1:n, D.states)])
    } else {
      p2 <- 0
      p9 <- 0
    }
    
    if (model.type=="crp"){
      p3 <- model.prior.crp.marginalized(apply(M, 1, match, x=1), prior.alpha.a, prior.alpha.b)
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
  
  #return posterior samples and marginal likelihood
  output <- list("y"=y, "prior.D"=prior.D, "prior.M"=prior.M, "prior.v.b"=prior.v.b, "samples"=samples, "samples.ml"=samples.ml, "Z"=Z, "W"=W)
  output <- c(output, results[names(results) != "post.hyperparameters"])
  output$ln.ml <- ln.ml
  output
}
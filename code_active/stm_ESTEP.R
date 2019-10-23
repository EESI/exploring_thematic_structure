LOCAL <- function (model, documents, nsims){
  
  cat('Running *corrected* local theta routine.\n')
  
  sigma <- model$sigma
  siginv <- model$invsigma
  mu <- model$mu$mu
  lambda <- model$eta
  beta <- lapply(model$beta$logbeta, exp)
  betaindex <- model$settings$covariates$betaindex
  out <- vector(mode = "list", length = nrow(lambda))
  for (i in 1:nrow(lambda)) {
    doc.ct <- documents[[i]][2, ]
    eta <- lambda[i, ]
    theta <- model$theta[i, ]
    doc.beta <- beta[[betaindex[i]]][, documents[[i]][1, 
                                                      ]]
    hess <- ln.hess(eta, theta, doc.beta, doc.ct, siginv)
    nu <- try(chol2inv(chol.default(hess)), silent = TRUE)
    if (class(nu) == "try-error") {
      
      if (NCOL(mu) == 1) mu.i <- mu else mu.i <- mu[,i]
      optim.out <- optim(par = eta, fn = lhoodcpp, gr = gradcpp, 
                         method = "BFGS", control = list(maxit = 500), 
                         doc_ct = doc.ct, mu = mu.i, siginv = siginv, 
                         beta = doc.beta)
      
      eta <- optim.out$par
      theta <- softmax(c(eta, 0))
      hess <- ln.hess(eta, theta, doc.beta, doc.ct, siginv)
      nu <- try(chol2inv(chol.default(hess)), silent = TRUE)
      if (class(nu) == "try-error") {
        newton.out <- newton(eta, doc.ct, mu.i, siginv, 
                             doc.beta, hess, max.its = 1000)
        nu <- newton.out$nu
        eta <- newton.out$eta
      }
    }
    mat <- rmvnorm(nsims, eta, nu)
    mat <- cbind(mat, 0)
    out[[i]] <- exp(mat - row.lse(mat))
  }
  return(out)
}

NEWTON <- function (eta, doc.ct, mu, siginv, beta, hess, max.its = 1000){
  
  cat('Running *corrected* Newton method.\n')
  
  its <- 0
  search <- function(x, dir, eta, ...) {
    lhoodcpp(eta = (eta + x * dir), ...)
  }
  while (its < max.its) {
    dir <- -solve(hess) %*% gradcpp(eta, beta, doc.ct, mu, siginv)
    maxint <- 2
    opt <- optimize(f = search, interval = c(-2, maxint), 
                    dir = dir, eta = eta, doc_ct = doc.ct, mu = mu, siginv = siginv, 
                    beta = beta, maximum = FALSE)
    while (opt$objective > lhoodcpp(eta, beta, doc.ct, mu, siginv)) {
      maxint <- min(maxint/2, opt$minimum - 0.00001 * opt$minimum)
      opt <- optimize(f = search, interval = c(0, maxint), 
                      dir = dir, eta = eta, doc_ct = doc.ct, mu = mu, 
                      siginv = siginv, beta = beta, maximum = FALSE)
    }

    eta <- eta + opt$minimum * dir
    expeta <- c(exp(eta), 1)
    theta <- expeta/sum(expeta)
    hess <- ln.hess(eta, theta, beta, doc.ct, siginv)
    nu <- try(chol2inv(chol.default(hess)), silent = TRUE)
    if (class(nu) != "try-error") {
      break
    }
    else {
      its <- its + 1
    }
  }
  if (its == max.its) {
    warning(sprintf("Document covariance matrix didn't converge after %i Newton iterations. Using nearest positive definite matrix to inverse hessian.", 
                    max.its))
    nu <- nearPD(solve(hess))
  }
  return(list(eta = eta, nu = nu))
}

ESTEP <- function(documents, beta.index, update.mu, #null allows for intercept only model  
                  beta, lambda.old, mu, sigma, 
                  verbose) {
  
  cat('Running *adaptable* estep.\n')
  
  #quickly define useful constants
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  ctevery <- ifelse(N>100, floor(N/100), 1)
  if(!update.mu) mu.i <- as.numeric(mu)
  
  # 1) Initialize Sufficient Statistics 
  sigma.ss <- diag(0, nrow=(K-1))
  beta.ss <- vector(mode="list", length=A)
  for(i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
  
  # 2) Precalculate common components
  sigobj <- try(chol.default(sigma), silent=TRUE)
  if(class(sigobj)=="try-error") {
    if (verbose) cat('\nTry error on chol()\nCalculating inv(sigma) via solve()\n')
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  }else{
    siginv <- chol2inv(sigobj)
    if(!any(is.nan(siginv))){
      if(verbose) cat('\nCalculating inv(sigma) via chol2inv(chol())\n')
      sigmaentropy <- sum(log(diag(sigobj)))
    }else{
      if (verbose) cat('\nNaN values on chol2inv(chol())\nCalculating inv(sigma) via solve()\n')
      siginv <- solve(sigma)
      sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    }
  }
  
  # 3) Document Scheduling
  # For right now we are just doing everything in serial.
  # the challenge with multicore is efficient scheduling while
  # maintaining a small dimension for the sufficient statistics.
  for(i in 1:N) {
    #update components
    doc <- documents[[i]]
    words <- doc[1,]
    aspect <- beta.index[i]
    init <- lambda.old[i,]
    if(update.mu) mu.i <- mu[,i]
    beta.i <- beta[[aspect]][,words,drop=FALSE]
    
    #infer the document
    doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, 
                                     doc=doc, sigmaentropy=sigmaentropy)
    
    # update sufficient statistics 
    sigma.ss <- sigma.ss + doc.results$eta$nu
    beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
    bound[i] <- doc.results$bound
    lambda[[i]] <- c(doc.results$eta$lambda)
    if(verbose && i%%ctevery==0) cat(".")
  }
  if(verbose) cat("\n") #add a line break for the next message.
  
  #4) Combine and Return Sufficient Statistics
  lambda <- do.call(rbind, lambda)
  return(list(sigma=sigma.ss, beta=beta.ss, bound=bound, lambda=lambda))
}

tmpfun <- get("estep", envir = asNamespace("stm"))
environment(ESTEP) <- environment(tmpfun)
attributes(ESTEP) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace("estep", ESTEP, ns="stm")

tmpfun <- get("thetapost.local", envir = asNamespace("stm"))
environment(LOCAL) <- environment(tmpfun)
attributes(LOCAL) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace("thetapost.local", LOCAL, ns="stm")

tmpfun <- get("newton", envir = asNamespace("stm"))
environment(NEWTON) <- environment(tmpfun)
attributes(NEWTON) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace("newton", NEWTON, ns="stm")
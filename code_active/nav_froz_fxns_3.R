doc.to.ijv <- function(documents, fixzeroindex=TRUE) {
  #Turns our format into triplet format (can be zero indexed)
  indices <- unlist(lapply(documents, '[',1,)) #grab the first row
  if((0 %in% indices) & fixzeroindex) indices <- indices + 1 #if zero-indexed, fix it.
  counts <- lapply(documents, '[',2,)  #grab the second row but preserve the list structure for a moment
  VsubD <- unlist(lapply(counts,length)) #grab the number of unique words per document
  rowsums <- unlist(lapply(counts,sum)) #grab the number of tokens per documents
  docids <- rep(1:length(documents), times=VsubD) #add row numbers
  counts <- unlist(counts) #unlist the count structure
  #now we return using the convention for the simple triplet matrix,
  #plus the row sums which we use in DMR.
  return(list(i=as.integer(docids), j=as.integer(indices), v=as.integer(counts), rowsums=as.integer(rowsums)))
}

gram <- function(mat) {
  nd <- rowSums(mat)
  mat <- mat[nd>=2,] #its undefined if we don't have docs of length 2
  nd <- nd[nd>=2]
  divisor <- nd*(nd-1)
  #clearer code is removed below for memory efficiency
  #Htilde <- mat/sqrt(divisor)
  #Hhat <- diag(colSums(mat/divisor))
  #Q <- crossprod(Htilde) - Hhat
  Q <- crossprod(mat/sqrt(divisor)) - diag(colSums(mat/divisor))
  #if(min(Q)<0) Q@x[Q@x < 0 ] <- 0
  return(as.matrix(Q))
}

fastAnchor <- function(Qbar, K, verbose=TRUE) {
  basis <- c()
  rowSquaredSums <- rowSums(Qbar^2) #StabilizedGS
  
  for(i in 1:K) {
    basis[i] <- which.max(rowSquaredSums) #83-94
    
    maxval <- rowSquaredSums[basis[i]]
    normalizer <- 1/sqrt(maxval) #100
    
    #101-103
    Qbar[basis[i],] <- Qbar[basis[i],]*normalizer 
    
    #For each row
    innerproducts <- Qbar%*%Qbar[basis[i],] #109-113
    
    #Each row gets multiplied out through the Basis
    project <- as.numeric(innerproducts)%o%Qbar[basis[i],] #118
    
    #Now we want to subtract off the projection but
    #first we should zero out the components for the basis
    #vectors which we weren't intended to calculate
    project[basis,] <- 0 #106 (if statement excluding basis vectors)
    Qbar <- Qbar - project #119
    rowSquaredSums <- rowSums(Qbar^2)
    rowSquaredSums[basis] <- 0 #here we cancel out the components we weren't calculating.
    if(verbose) cat(".")
  }
  return(basis)
}

recoverL2 <- function(Qbar, anchor, wprob, verbose=TRUE, ...) {
  #NB: I've edited the script to remove some of the calculations by commenting them
  #out.  This allows us to store only one copy of Q which is more memory efficient.
  #documentation for other pieces is below.
  #' \item{R}{a matrix of dimensions K by K that contains the topic covariances.}
  #' \item{condprob}{a list of exponentiated gradient results.  useful for checking convergence.}
  
  #Qbar <- Q/rowSums(Q)
  X <- Qbar[anchor,]
  XtX <- tcrossprod(X)
  
  #Word by Word Solve For the Convex Combination
  condprob <- vector(mode="list", length=nrow(Qbar))
  for(i in 1:nrow(Qbar)) {
    if(i %in% anchor) { 
      #if its an anchor we create a dummy entry that is 1 by definition
      vec <- rep(0, nrow(XtX))
      vec[match(i,anchor)] <- 1
      condprob[[i]] <- list(par=vec)
    } else {
      y <- Qbar[i,]
      condprob[[i]] <- expgrad(X,y,XtX, ...)
    }
  }
  #Recover Beta (A in this notation)
  #  Now we have p(z|w) but we want the inverse
  weights <- lapply(condprob, function(x) x$par)
  weights <- do.call(rbind, weights)
  A <- weights*wprob
  A <- t(A)/colSums(A)
  
  #Recover The Topic-Topic Covariance Matrix
  #Adag <- mpinv(A)  
  #R <- t(Adag)%*%Q%*%Adag
  return(list(A=A))
  #return(list(A=A, R=R, condprob=condprob))
}

#' Exponentiated Gradient with L2 Loss
#'
#' Find the optimal convex combination of features X which approximate the vector y under
#' and L2 loss.
#' 
#' An implementation of RecoverL2 based on David Mimno's code.  In this setting the
#' objective and gradient can both be kernalized leading to faster computations than
#' possible under the KL loss.  Specifically the computations are no longer dependent
#' on the size of the vocabulary making the operation essentially linear on the total
#' vocab size.
#' The matrix X'X can be passed as an argument in settings
#' like spectral algorithms where it is constant across multiple optimizations.  
#'
#' @param X the transposed feature matrix.
#' @param y the target vector
#' @param XtX optionally a precalculated crossproduct.
#' @param alpha an optional initialization of the parameter.  Otherwise it starts at 1/nrow(X).
#' @param tol convergence tolerance
#' @param max.its maximum iterations to run irrespective of tolerance
#' @return 
#' \item{par}{optimal weights}
#' \item{its}{number of iterations run}
#' \item{converged}{logical indicating if it converged}
#' \item{entropy}{entropy of the resulting weights}
#' \item{log.sse}{log of the sum of squared error}
#' @export
expgrad <- function(X, y, XtX=NULL, alpha=NULL, tol=1e-7, max.its=500) {
  if(is.null(alpha)) alpha <- 1/nrow(X) 
  alpha <- matrix(alpha, nrow=1, ncol=nrow(X))
  if(is.null(XtX)) XtX <- tcrossprod(X)
  
  ytX <- y%*%t(X)
  converged <- FALSE
  eta <- 50
  sse.old <- Inf
  its <- 1
  while(!converged) {
    #find the gradient (y'X - alpha'X'X)
    grad <- (ytX - alpha%*%XtX) #101-105
    sse <- sum(grad^2) #106, sumSquaredError
    grad <- 2*eta*grad
    maxderiv <- max(grad)    
    
    #update parameter 
    alpha <- alpha*exp(grad-maxderiv)
    #project parameter back to space
    alpha <- alpha/sum(alpha)
    
    converged <- abs(sqrt(sse.old)-sqrt(sse)) < tol
    if(its==max.its) break
    sse.old <- sse
    its <- its + 1
  } 
  entropy <- -1*sum(alpha*log(alpha))
  return(list(par=as.numeric(alpha), its=its, converged=converged,
              entropy=entropy, log.sse=log(sse)))
}

#calculate the Moore Penrose Pseudo-Inverse
# adapted from the mASS function ginv
mpinv <- function (X) {
  
  if (!is.matrix(X)) X <- as.matrix(X)
  
  tol <- sqrt(.Machine$double.eps)
  Xsvd <- svd(X)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                                               t(Xsvd$u[, Positive, drop = FALSE]))
}

tsneAnchor <- function(Qbar) {
  if(!(requireNamespace("Rtsne",quietly=TRUE) & requireNamespace("geometry", quietly=TRUE))){
    stop("Please install the Rtsne and geometry packages to use this setting.")
  } 
  #project to 3-D
  proj <- try(Rtsne::Rtsne(Qbar, dims=3) , silent=TRUE)
  if(class(proj)=="try-error") {
    #if this failed it is probably duplicates which Rtsne cannot handle
    dup <- duplicated(Qbar)
    if(!any(dup)) stop("an unknown error has occured in Rtsne")
    
    dup <- which(dup)
    for(r in dup) {
      row <- Qbar[r,]
      row[row>0] <- runif(sum(row>0),0,1e-5) # add a bit of noise to non-zero duplicates
      row <- row/sum(row) #renormalize
      Qbar[r,] <- row
    }
    #and now do it again
    proj <- Rtsne::Rtsne(Qbar, dims=3) 
  }
  
  hull <- geometry::convhulln(proj$Y)
  anchor <- sort(unique(c(hull)))
  return(anchor)
}

kappa.init <- function(documents, K, V, A, interactions) {
  kappa.out <- list()
  #Calculate the baseline log-probability (m)
  freq <- matrix(unlist(documents),nrow=2) #break it into a matrix
  freq <- split(freq[2,], freq[1,]) #shift into list by word type
  m <- unlist(lapply(freq, sum)) #sum over the word types
  m <- m/sum(m)
  #m <- log(m)
  m <- log(m) - log(mean(m)) #logit of m
  kappa.out$m <- m
  
  #Defining parameters
  aspectmod <- A > 1
  if(aspectmod) {
    interact <- interactions 
  } else {
    interact <- FALSE
  }
  
  #Create the parameters object
  parLength <- K + A*aspectmod + (K*A)*interact
  kappa.out$params <- vector(mode="list",length=parLength)
  for(i in 1:length(kappa.out$params)) {
    kappa.out$params[[i]] <- rep(0, V)
  }
  
  #Create a running sum of the kappa parameters starting with m
  kappa.out$kappasum <- vector(mode="list", length=A)
  for (a in 1:A) {
    kappa.out$kappasum[[a]] <- matrix(m, nrow=K, ncol=V, byrow=TRUE)
  }
  
  #create covariates. one element per item in parameter list.
  #generation by type because its conceptually simpler
  if(!aspectmod & !interact) {
    kappa.out$covar <- list(k=1:K, a=rep(NA, parLength), type=rep(1,K))
  }
  if(aspectmod & !interact) {
    kappa.out$covar <- list(k=c(1:K,rep(NA,A)), a=c(rep(NA, K), 1:A), type=c(rep(1,K), rep(2,A)))      
  }
  if(interact) {
    kappa.out$covar <- list(k=c(1:K,rep(NA,A), rep(1:K,A)), 
                            a=c(rep(NA, K), 1:A, rep(1:A,each=K)), 
                            type=c(rep(1,K), rep(2,A), rep(3,K*A)))            
  }
  return(kappa.out)
}

stm.init <- function(documents, settings) {
  
  K <- settings$dim$K
  V <- settings$dim$V
  A <- settings$dim$A
  N <- settings$dim$N
  mode <- settings$init$mode
  nits <- settings$init$nits 
  alpha <- settings$init$alpha 
  eta <- settings$init$eta 
  burnin <- settings$init$burnin 
  
  #mode=spectral
  verbose <- settings$verbose
  # (1) Prep the Gram matrix
  docs <- doc.to.ijv(documents)
  mat <- sparseMatrix(docs$i,docs$j, x=docs$v)
  rm(docs)
  wprob <- colSums(mat)
  wprob <- wprob/sum(wprob)
  
  Q <- gram(mat)
  #verify that there are no zeroes
  keep <- NULL
  Qsums <- rowSums(Q)
  
  if(any(Qsums==0)) {
    #if there are zeroes, we want to remove them for just the anchor word procedure.
    temp.remove <- which(Qsums==0)
    keep <- which(Qsums!=0)
    Q <- Q[-temp.remove,-temp.remove]
    Qsums <- Qsums[-temp.remove]
    wprob <- wprob[-temp.remove]
  }
  
  Q <- Q/Qsums
  
  # (2) anchor words
  anchor <- fastAnchor(Q, K=K, verbose=verbose)
  # (3) recoverKL
#   beta <- recoverL2(Q, anchor, wprob, verbose=verbose)$A
  
#   if(!is.null(keep)) {
#     #if there were zeroes, reintroduce them
#     #assign missing compoinents the machine double epsilon
#     #and renormalize just in case.
#     beta.new <- matrix(0, nrow=K, ncol=V)
#     beta.new[,keep] <- beta
#     beta.new[,temp.remove] <- .Machine$double.eps 
#     beta <- beta.new/rowSums(beta.new)  
#     rm(beta.new)
#   }
  
  # (4) generate other parameters
  mu <- matrix(0, nrow=(K-1),ncol=1)
  sigma <- diag(20, nrow=(K-1))
  lambda <- matrix(0, nrow=N, ncol=(K-1))
  
  #turn beta into a list and assign it for each aspect
#   beta <- rep(list(beta),A)
#   model <- list(mu=mu, sigma=sigma, beta=beta, lambda=lambda)
  model <- list(mu=mu, sigma=sigma, lambda=lambda)
  #initialize the kappa vectors
  if(!settings$kappa$LDAbeta) {
    model$kappa <- kappa.init(documents, K, V, A, interactions=settings$kappa$interactions)
  }
  return(model)
}

estep <- function(documents, beta.index, update.mu, #null allows for intercept only model  
                  beta, lambda.old, mu, sigma, 
                  verbose) {
  
  #quickly define useful constants
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  ctevery <- ifelse(N>100, floor(N/100), 1)
  if(!update.mu) mu.i <- as.numeric(mu)
  
  # 1) Initialize Sufficient Statistics 
  sigma.ss <- diag(0, nrow=(K-1))
#   beta.ss <- vector(mode="list", length=A)
#   for(i in 1:A) {
#     beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
#   }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
  
  # 2) Precalculate common components
  sigobj <- try(chol.default(sigma), silent=TRUE)
  if(class(sigobj)=="try-error") {
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
    sigmaentropy <- sum(log(diag(sigobj)))
    siginv <- chol2inv(sigobj)
  }
  # 3) Document Scheduling
  # For right now we are just doing everything in serial.
  # the challenge with multicore is efficient scheduling while
  # maintaining a small dimension for the sufficient statistics.
  Z <- list()
  doc_sums <- matrix(0,N,K)
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
    
    phi <- t(t(doc.results$phis)/colSums(doc.results$phis)) # K dim multinom par, therefore topics on the simplex?
    z <- sapply(1:ncol(phi), function(w) which(rmultinom(1,1:K,phi[,w])==1))
    Z[[i]] <- z
    doc_sums[i,z] <- doc_sums[i,z] + 1
    
    # update sufficient statistics 
    sigma.ss <- sigma.ss + doc.results$eta$nu
#     beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
    bound[i] <- doc.results$bound
    lambda[[i]] <- c(doc.results$eta$lambda)
    if(verbose && i%%ctevery==0) cat(".")
  }
  if(verbose) cat("\n") #add a line break for the next message.
  
  #4) Combine and Return Sufficient Statistics
  lambda <- do.call(rbind, lambda)
  return(list(sigma=sigma.ss, bound=bound, lambda=lambda, doc_sums=doc_sums, Z=Z))
}

logisticnormalcpp <- function(eta, mu, siginv, beta, doc, sigmaentropy) {
  doc.ct <- doc[2,]
  Ndoc <- sum(doc.ct)
  #even at K=100, BFGS is faster than L-BFGS
  optim.out <- optim(par=eta, fn=lhoodcpp, gr=gradcpp,
                     method="BFGS", control=list(maxit=500),
                     doc_ct=doc.ct, mu=mu,
                     siginv=siginv, beta=beta)
  
  #Solve for Hessian/Phi/Bound returning the result
  hpbcpp(optim.out$par, doc_ct=doc.ct, mu=mu,
         siginv=siginv, beta=beta,
         sigmaentropy=sigmaentropy)
}

lhoodcpp <- function(eta, beta, doc_ct, mu, siginv) {
  .Call('stm_lhoodcpp', PACKAGE = 'stm', eta, beta, doc_ct, mu, siginv)
}

gradcpp <- function(eta, beta, doc_ct, mu, siginv) {
  .Call('stm_gradcpp', PACKAGE = 'stm', eta, beta, doc_ct, mu, siginv)
}

hpbcpp <- function(eta, beta, doc_ct, mu, siginv, sigmaentropy) {
  .Call('stm_hpbcpp', PACKAGE = 'stm', eta, beta, doc_ct, mu, siginv, sigmaentropy)
}

opt.mu <- function(lambda, mode=c("CTM","Pooled", "L1"), covar=NULL, enet=NULL) {
  
  #When there are no covariates we use the CTM method
  if(mode=="CTM") {
    mu <- matrix(colMeans(lambda), ncol=1)
    return(list(mu=mu, gamma=NULL))
  }
}

opt.sigma <- function(nu, lambda, mu, sigprior) {  
  
  #find the covariance
  if(ncol(mu)==1) {
    covariance <- crossprod(sweep(lambda, 2, STATS=as.numeric(mu), FUN="-"))
  } else {
    covariance <- crossprod(lambda-t(mu)) 
  }
  sigma <- (covariance + nu)/nrow(lambda) #add to estimation variance
  sigma <- diag(diag(sigma),nrow=nrow(nu))*sigprior + (1-sigprior)*sigma #weight by the prior
  return(sigma)
}

opt.beta <- function(beta.ss, kappa, settings) {
  #if its standard lda just row normalize
  if(is.null(kappa)) return(list(beta=list(beta.ss[[1]]/rowSums(beta.ss[[1]]))))
}

convergence.check <- function(bound.ss, convergence, settings) {
  
  #unpack the relevant pieces from the settings
  verbose <- settings$verbose
  emtol <- settings$convergence$em.converge.thresh
  maxits <- settings$convergence$max.em.its
  
  #initialize the convergence object if empty
  if(is.null(convergence)) convergence <- list(bound=c(), its=1, converged=FALSE, stopits=FALSE)
  
  #fill in the current bound
  convergence$bound[convergence$its] <- sum(bound.ss)
  
  #if not the first iteration
  if(convergence$its > 1) {
    old <- convergence$bound[convergence$its-1]
    new <- convergence$bound[convergence$its]
    convergence.check <- (new-old)/abs(old)
    #if(convergence.check < emtol & convergence.check > 0) {
    if(convergence.check < emtol) {
      convergence$converged <- TRUE
      convergence$stopits <- TRUE
      if(verbose) cat("Model Converged \n")
      return(convergence)
    }
  }
  
  if(convergence$its==maxits) {
    if(verbose) cat("Model Terminated Before Convergence Reached \n")
    convergence$stopits <- TRUE
    return(convergence)
  }
  convergence$its <- convergence$its + 1
  return(convergence)
}

row.lse <- function(mat) {
  matrixStats::rowLogSumExps(mat)
}

ctm_frozen <- function(fit,documents,vocabulary,seed=NULL,max.em.its=500,emtol=1e-5,
                       avg_iters=15,verbose=TRUE,reportevery=5,true_doc_content=TRUE,
                       data=data,covariate=NULL,
                       parallel=TRUE,nc=60){
   
   #beta_frozen <- list(exp(fit$beta$logbeta[[1]]))
   #beta_frozen <- lapply(fit$beta$logbeta,exp) # i moved this below after stm.init
   
   K <- fit$settings$dim$K
   N <- length(documents)
   
   
   if (is.null(covariate)){
      if (!is.null(fit$settings$call$content)) content <- TRUE else content <- FALSE
   }else if (covariate != 0){
      content <- FALSE
      fit$beta$logbeta <- list(fit$beta$logbeta[[covariate]])
   }else{
      content <- FALSE
      fit$beta$logbeta <- list(do.call('rbind',fit$beta$logbeta))
      K <- nrow(list(do.call('rbind',fit$beta$logbeta))[[1]])
   }
   
   
   prevalence <- NULL
   init.type <- 'Spectral' # if I also set K=0, then estimates topics via Lee and Mimno 2014
     
   interactions <- TRUE # only relevant for covariates
   ngroups <- 1
   model <- NULL # allows one to restart an existing model
   sigma.prior <- 0 # strength of regularization for diagnolized covariance matrix
   gamma.prior <- c("Pooled", "L1")
   kappa.prior <- c("L1", "Jeffreys")
   control <- list()
   
   wcountvec <- unlist(lapply(documents, function(x) rep(x[1,], times=x[2,])),use.names=FALSE)
   wcounts <- list(Group.1=sort(unique(wcountvec)))
   V <- length(wcounts$Group.1)
   wcounts$x <- tabulate(wcountvec)
   rm(wcountvec)
   
   xmat <- NULL
   
   
   if (content){
      termobj <- terms(eval(fit$settings$call$content),data=data)
      if(is.null(data)) stop('Make sure you set the data used for fit.')
      char <- rownames(attr(termobj, "factors"))[1]
      yvar <- as.factor(data[[char]])
      yvarlevels <- levels(yvar)
      betaindex <- as.numeric(yvar)
   }else{
      yvarlevels <- NULL
      betaindex <- rep(1, length(documents))
   }
   
   
   
   A <- length(unique(betaindex)) #define the number of aspects
   
   settings <- list(dim=list(K=K, A=A, 
                             V=V, N=N, wcounts=wcounts),
                    verbose=verbose,
                    topicreportevery=reportevery,
                    convergence=list(max.em.its=max.em.its, em.converge.thresh=emtol),
                    covariates=list(X=xmat, betaindex=betaindex, yvarlevels=yvarlevels),
                    gamma=list(mode="CTM", prior=NULL, enet=1), # gamma is for prev so ignore match.arg(gamma.prior)
                    sigma=list(prior=sigma.prior),
                    kappa=list(LDAbeta=TRUE, interactions=FALSE, 
                               fixedintercept=TRUE, mstep=list(tol=.001, maxit=3),
                               contrast=FALSE),
                    tau=list(mode=NULL, tol=1e-5,
                             enet=1,nlambda=250, lambda.min.ratio=.001, ic.k=2,
                             maxit=1e4),
                    init=list(mode=init.type, nits=50, burnin=25, alpha=(50/K), eta=.01,
                              s=.05, p=3000, d.group.size=2000), 
                    seed=seed,
                    ngroups=ngroups)
   
   if (content){
      settings$kappa$LDAbeta <- FALSE #can't do LDA topics with a covariate 
   }else{
      settings$kappa$interactions <- FALSE #can't have interactions without a covariate.
   }
   
   
   legalargs <-  c("tau.maxit", "tau.tol", 
                   "fixedintercept","kappa.mstepmaxit", "kappa.msteptol", 
                   "kappa.enet", "nlambda", "lambda.min.ratio", "ic.k", "gamma.enet",
                   "nits", "burnin", "alpha", "eta", "contrast",
                   "rp.s", "rp.p", "rp.d.group.size", "SpectralRP")
   
   # Process the Seed
   if(is.null(settings$seed)) {
      #if there is no seed, choose one and set it, recording for later
      seed <- floor(runif(1)*1e7) 
      set.seed(seed)
      settings$seed <- seed
   } else {
      #otherwise just use the provided seed.
      set.seed(settings$seed)
   }
   
   
   globaltime <- proc.time()  # maybe remove?
   #verbose <- settings$verbose
   #ngroups <- settings$ngroups # probably remove
   
   #initialize
   model <- stm.init(documents, settings)
   beta <- list(beta=lapply(fit$beta$logbeta,exp))
   beta$kappa <- fit$beta$kappa
   
   #unpack
   mu <- list(mu=model$mu)
   sigma <- model$sigma
   
   #beta <- list(beta=model$beta)
   #if(!is.null(model$kappa)) beta$kappa <- model$kappa
   
   lambda <- model$lambda
   convergence <- NULL 
   #discard the old object
   rm(model)
   
   #Pull out some book keeping elements
   ntokens <- sum(settings$dim$wcounts$x)
   betaindex <- settings$covariates$betaindex # since no covariates, intercept? therefore 1?
   stopits <- FALSE
   suffstats <- vector(mode="list", length=ngroups)
   
   ############
   #Step 2: Run EM
   ############
   while(!stopits) {
      
      #####
      # Non-Blocked Updates
      #####
      
      t1 <- proc.time()
      
      #run the model
      suffstats <- estep(documents=documents, beta.index=betaindex, 
                         update.mu=(!is.null(mu$gamma)),  
                         beta$beta, lambda, mu$mu, sigma, 
                         verbose)
      
      t1 <- proc.time()
      
      Z <- suffstats$Z
      sigma.ss <- suffstats$sigma
      lambda <- suffstats$lambda
      #   beta.ss <- suffstats$beta
      bound.ss <- suffstats$bound
      
      #do the m-step
      mu <- opt.mu(lambda=lambda, mode=settings$gamma$mode, 
                   covar=settings$covariates$X, settings$gamma$enet)
      sigma <- opt.sigma(nu=sigma.ss, lambda=lambda, 
                         mu=mu$mu, sigprior=settings$sigma$prior)
      #   beta <- opt.beta(beta.ss, beta$kappa, settings)
      
      
      #Convergence
      convergence <- convergence.check(bound.ss, convergence, settings)
      stopits <- convergence$stopits
      
      cat('converged:',convergence$converged,'\n')
      
   }
   
   ############
   #Step 3: Sample Z from post
   ############
   
   samp_posterior <- samp_z_parallel(documents=documents, beta.index=betaindex, 
                                     update.mu=(!is.null(mu$gamma)),  
                                     beta$beta, lambda, mu$mu, sigma, 
                                     avg_iters,
                                     verbose,
                                     true_doc_content=true_doc_content,
                                     parallel=parallel,nc=nc)
   Z_bar <- samp_posterior$Z_bar
   doc_sums <- samp_posterior$doc_sums
   #topic_sums <- samp_posterior$topic_sums
   #doc_sums_replicates <- samp_posterior$doc_sums_replicates
   #beta_sums_replicates <- samp_posterior$beta_sums_replicates
   
   #######
   #Step 3: Construct Output
   #######
   time <- (proc.time() - globaltime)[3]
   #convert the beta back to log-space
   beta$logbeta <- beta$beta
   for(i in 1:length(beta$logbeta)) {
      beta$logbeta[[i]] <- log(beta$logbeta[[i]])
   }
   beta$beta <- NULL
   lambda <- cbind(lambda,0)
   model <- list(mu=mu, sigma=sigma, beta=beta, 
                 Z=Z,Z_bar=Z_bar,
                 doc_sums=doc_sums,
                 #topic_sums=topic_sums,beta_sums_replicates=beta_sums_replicates,doc_sums_replicates=doc_sums_replicates,
                 settings=settings,
                 vocabulary=vocabulary, convergence=convergence, 
                 theta=exp(lambda - row.lse(lambda)), 
                 eta=lambda[,-ncol(lambda), drop=FALSE],
                 invsigma=solve(sigma), time=time, version=utils::packageDescription("stm")$Version)
   class(model) <- "STM"  
   
   return(model)
}


samp_z <- function(documents, beta.index, update.mu, #null allows for intercept only model  
                   beta, lambda.old, mu, sigma, 
                   avg_iters,
                   verbose,
                   true_doc_content=TRUE) {
   
   
   V <- ncol(beta[[1]])
   K <- nrow(beta[[1]])
   N <- length(documents)
   A <- length(beta)
   if(!update.mu) mu.i <- as.numeric(mu)
   
   sigobj <- try(chol.default(sigma), silent=TRUE)
   if(class(sigobj)=="try-error") {
      sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
      siginv <- solve(sigma)
   } else {
      sigmaentropy <- sum(log(diag(sigobj)))
      siginv <- chol2inv(sigobj)
   }
   
   doc_sums_total <- matrix(0,N,K)
   topic_sums_total <- array(0,c(K,N,V))
   doc_sums_replicates <- array(0,c(avg_iters,N,K))
   beta_sums_replicates <- array(0,c(avg_iters,K,V))
   for (avg in 1:avg_iters){
      
      cat('\nexploring posterior of z for',avg_iters,'iterations:',avg,'\n')
      
      for(i in 1:N) {
         #update components
         doc <- documents[[i]]
         words <- doc[1,]
         words_counts <- doc[2,]
         aspect <- beta.index[i]
         init <- lambda.old[i,]
         if(update.mu) mu.i <- mu[,i]
         beta.i <- beta[[aspect]][,words,drop=FALSE]
         
         #infer the document
         doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, 
                                          doc=doc, sigmaentropy=sigmaentropy)
         
         phi <- t(t(doc.results$phis)/colSums(doc.results$phis)) # K dim multinom par, therefore topics on the simplex?
         
         for (j in 1:length(words)){
            
            if (true_doc_content){
               z <- rowSums(rmultinom(words_counts[j],1:K,phi[,j]))
            }else{
               z <- rmultinom(1,1:K,phi[,j]) # check if orientation is correct
            }
            
            ww <- words[j]
            
            beta_sums_replicates[avg,,ww] <- beta_sums_replicates[avg,,ww]+ z
            doc_sums_replicates[avg,i,] <- doc_sums_replicates[avg,i,] + z
            doc_sums_total[i,] <- doc_sums_total[i,] + z
            topic_sums_total[,i,ww] <- topic_sums_total[,i,ww] + z     
         }
         
         if(verbose) cat(".")
      }
      
   }
   
   cat('\nexploration of posterior of z completed','\n')
   return(list(Z_bar=doc_sums_total/rowSums(doc_sums_total),
               doc_sums_replicates=doc_sums_replicates,
               beta_sums_replicates=beta_sums_replicates,
               topic_sums=topic_sums_total))
}

samp_z_parallel <- function(documents, beta.index, update.mu, #null allows for intercept only model  
                            beta, lambda.old, mu, sigma, 
                            avg_iters,
                            verbose,
                            true_doc_content=TRUE,
                            parallel, nc) {
   
   
   V <- ncol(beta[[1]])
   K <- nrow(beta[[1]])
   N <- length(documents)
   A <- length(beta)
   if(!update.mu) mu.i <- as.numeric(mu)
   
   sigobj <- try(chol.default(sigma), silent=TRUE)
   if(class(sigobj)=="try-error") {
      sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
      siginv <- solve(sigma)
   } else {
      sigmaentropy <- sum(log(diag(sigobj)))
      siginv <- chol2inv(sigobj)
   }
   
   
   if (parallel){
      
      require(snow)
      require(doSNOW)
   
      cl <- makeCluster(nc,type='SOCK')
      clusterExport(cl,c('logisticnormalcpp','gradcpp','lhoodcpp','hpbcpp'))
      registerDoSNOW(cl)
      
      for (avg in 1:avg_iters){
               
         doc_sums_total <- matrix(0,N,K)
         

         cat('\nexploring posterior of z for',avg_iters,'iterations:',avg,'\n')
         explore_z <- foreach(i=1:N,.verbose=FALSE,.errorhandling="stop",.packages=c("Matrix","stm"),.combine='rbind') %dopar% {
            #update components
            doc_sum_vec <- numeric(K)
            
            doc <- documents[[i]]
            words <- doc[1,]
            words_counts <- doc[2,]
            aspect <- beta.index[i]
            init <- lambda.old[i,]
            if(update.mu) mu.i <- mu[,i]
            beta.i <- beta[[aspect]][,words,drop=FALSE]
            
            #infer the document
            doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, 
                                             doc=doc, sigmaentropy=sigmaentropy)
            
            phi <- t(t(doc.results$phis)/colSums(doc.results$phis)) # K dim multinom par, therefore topics on the simplex?
            
            for (j in 1:length(words)){
               
               if (true_doc_content){
                  z <- rowSums(rmultinom(words_counts[j],1:K,phi[,j]))
               }else{
                  z <- rmultinom(1,1:K,phi[,j]) # check if orientation is correct
               }
               
               doc_sum_vec <- doc_sum_vec + z    
            }
            
            doc_sum_vec
         }
         
         doc_sums_total <- doc_sums_total + explore_z
         
      }
      
      if(exists('cl')) stopCluster(cl)
      
   }else{
      
      doc_sums_total <- matrix(0,N,K)
      for (avg in 1:avg_iters){
         
         cat('\nexploring posterior of z for',avg_iters,'iterations:',avg,'\n')
         
         for(i in 1:N) {
            #update components
            doc <- documents[[i]]
            words <- doc[1,]
            words_counts <- doc[2,]
            aspect <- beta.index[i]
            init <- lambda.old[i,]
            if(update.mu) mu.i <- mu[,i]
            beta.i <- beta[[aspect]][,words,drop=FALSE]
            
            #infer the document
            doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, 
                                             doc=doc, sigmaentropy=sigmaentropy)
            
            phi <- t(t(doc.results$phis)/colSums(doc.results$phis)) # K dim multinom par, therefore topics on the simplex?
            
            for (j in 1:length(words)){
               
               if (true_doc_content){
                  z <- rowSums(rmultinom(words_counts[j],1:K,phi[,j]))
               }else{
                  z <- rmultinom(1,1:K,phi[,j]) # check if orientation is correct
               }
               
               ww <- words[j]
               
               doc_sums_total[i,] <- doc_sums_total[i,] + z 
            }
            
            if(verbose) cat(".")
         }
         
      }
      
   }
   
   cat('\nexploration of posterior of z completed','\n')
   return(list(Z_bar=doc_sums_total/rowSums(doc_sums_total),
               doc_sums=doc_sums_total))
}


balance_data <- function(metadata,freq,label,prop,R=FALSE){
  set_meta_table <- floor(table(metadata)*prop)
  set_meta_floor <- which.min(set_meta_table)
  set_samp_1 <- sample(which(metadata == names(set_meta_floor)),
                       set_meta_table[set_meta_floor],
                       replace=R)
  set_samp_2 <- sample(which(metadata != names(set_meta_floor)),
                       set_meta_table[set_meta_floor],
                       replace=R)
  set_samp_idx <- c(set_samp_1,set_samp_2)
  set_samp <- freq[set_samp_idx,]
  set_labels <- ifelse(metadata[set_samp_idx]==label,1,0)
  
  return(list(freq=set_samp,labels=set_labels))
}

collapse_beta <- function(otu_taxa,data_mat,metadata,topic_coefs=NULL,taxon=2,top=20){
  otu_taxa <- otu_taxa[rownames(data_mat),]
  target_taxon <- data.frame(long=apply(otu_taxa[,1:taxon],1,function(i) paste0(i,collapse='|')),
                             lowest_name=otu_taxa[,taxon],
                             id=rownames(otu_taxa))
  
  target_length <- nrow(unique(target_taxon))
  taxon_table <- data.frame(data_mat,taxon=target_taxon$long) %>%
    group_by(taxon) %>%
    summarise_each(funs(sum)) 
  
  tax_short <- substring(str_to_lower(colnames(otu_taxa)[taxon]),1,3)
  taxon_ids <- data.frame(long=taxon_table$taxon,
                          ids=paste0(tax_short,1:nrow(taxon_table)),
                          row.names=paste0(tax_short,1:nrow(taxon_table)))
  taxon_table <- dplyr::select(taxon_table,-taxon) 
  
  rownames(taxon_table) <- rownames(taxon_ids)
  
  taxon_table <- t(t(taxon_table)/colSums(taxon_table))
  taxon_table <- taxon_table[rowSums(taxon_table) != 0,]
  
  if (nrow(taxon_table) > top){
     if (!is.null(topic_coefs)){
        top_split <- 0
        top_start <- top
        coef_sig_pos_idx <- as.integer(names(topic_coefs)[topic_coefs>0])
        coef_sig_neg_idx <- as.integer(names(topic_coefs)[topic_coefs<0])
        
        pos_mean <- if (length(coef_sig_pos_idx) > 1) rowMeans(taxon_table[,coef_sig_pos_idx]) else taxon_table[,coef_sig_pos_idx]
        neg_mean <- if (length(coef_sig_neg_idx) > 1) rowMeans(taxon_table[,coef_sig_neg_idx]) else taxon_table[,coef_sig_neg_idx]
        
        while (length(top_split) < top_start){
           
           if ((length(pos_mean) > 0) & (length(pos_mean) > 0)){
              top_split <- unique(c(order(pos_mean,decreasing=TRUE)[1:(top/2)],
                                    order(neg_mean,decreasing=TRUE)[1:(top/2)]))
           }else if(length(pos_mean) == 0){
              top_split <- order(neg_mean,decreasing=TRUE)[1:(top/2)] 
           }else{
              top_split <- order(pos_mean,decreasing=TRUE)[1:(top/2)]    
           }
           top <- top + 2
        }
        taxon_table <- taxon_table[top_split,]
     }else{
        taxon_table <- taxon_table[order(rowMeans(taxon_table),decreasing=TRUE)[1:top],]
     }
  }

  rownames(taxon_table) <- taxon_ids[rownames(taxon_table),]$long
  
  return(taxon_table)
}

collapse_beta_kegg <- function(kegg_meta, kegg_mat, topic_coefs=NULL, top=12, method='mean', var_len=1, var_min=10^-5, filter_unknown=TRUE){
  
  kegg_mat <- kegg_mat[rowSums(kegg_mat) != 0,]
  kegg_mat_names <- rownames(kegg_mat)
  
  pw_unique <- unique(unlist(kegg_meta))
  pw_unique <- pw_unique[pw_unique != ""]
  
  collapse_mat <- matrix(0,ncol(kegg_mat),length(pw_unique),
                         dimnames=list(paste0('T',1:ncol(kegg_mat)),pw_unique))
  for (i in 1:nrow(kegg_mat)){
    id <- kegg_mat_names[i]
    pw <- kegg_meta[[id]]
    pw <- pw[pw != ""]
    collapse_mat[,pw] <- collapse_mat[,pw] + kegg_mat[id,]
  }
  collapse_mat <- t(collapse_mat/rowSums(collapse_mat))
  
  if (filter_unknown){
    remove_words <- c('Function unknown','None')
    remove_idx <- which(stringr::str_detect(rownames(collapse_mat),
                                      paste(remove_words,collapse='|')))
    collapse_mat <- collapse_mat[-remove_idx,]
  }
  
  if (method == 'mean'){
    if (nrow(collapse_mat) > top){
      if (!is.null(topic_coefs)){
        top_split <- 0
        top_start <- top
        coef_sig_pos_idx <- as.integer(names(topic_coefs)[topic_coefs>0])
        coef_sig_neg_idx <- as.integer(names(topic_coefs)[topic_coefs<0])
        
        coll_pos <- collapse_mat[,coef_sig_pos_idx]
        coll_neg <- collapse_mat[,coef_sig_neg_idx]
        if (length(coef_sig_pos_idx) == 1) coll_pos <- matrix(coll_pos,ncol=1) ### testing
        if (length(coef_sig_neg_idx) == 1) coll_neg <- matrix(coll_neg,ncol=1) ### testing
        
        while(length(top_split) < top_start){
          top_split <- unique(c(order(rowMeans(coll_pos),decreasing=TRUE)[1:(top/2)],
                                order(rowMeans(coll_neg),decreasing=TRUE)[1:(top/2)]))
          top <- top + 2
        }
        collapse_mat <- collapse_mat[top_split,]
      }else{
        collapse_mat <- collapse_mat[order(rowMeans(collapse_mat),decreasing=TRUE)[1:top],]
      }
    }
  }else if(method == 'sd'){
    if (nrow(collapse_mat) > top){
      if (!is.null(topic_coefs)){
        top_split <- 0
        top_start <- top
        coef_sig_pos_idx <- as.integer(names(topic_coefs)[topic_coefs>0])
        coef_sig_neg_idx <- as.integer(names(topic_coefs)[topic_coefs<0])
         
        coll_pos <- collapse_mat[,coef_sig_pos_idx]
        coll_neg <- collapse_mat[,coef_sig_neg_idx]
        if (length(coef_sig_pos_idx) == 1) coll_pos <- matrix(coll_pos,ncol=1) ### testing
        if (length(coef_sig_neg_idx) == 1) coll_neg <- matrix(coll_neg,ncol=1) ### testing
        
        while(length(top_split) < top_start){
          top_split <- unique(c(order(apply(coll_pos,1,
                                            function(x) ifelse(sum(x > var_min) > var_len,sd(x),0)),
                                      decreasing=TRUE)[1:(top/2)],
                                order(apply(coll_neg,1,
                                            function(x) ifelse(sum(x > var_min) > var_len,sd(x),0)),
                                      decreasing=TRUE)[1:(top/2)]))
          top <- top + 2
        }
        collapse_mat <- collapse_mat[top_split,]
      }else{
        collapse_mat <- collapse_mat[order(apply(collapse_mat,1,
                                                 function(x) ifelse(sum(x > var_min) > var_len,sd(x),0)),
                                           decreasing=TRUE)[1:top],]
      }
    }
  }
  
  return(collapse_mat)
}

pw_counter <- function(kegg_metadata,kegg_counts,coef_sig){
  pw_names <- unique(unlist(kegg_metadata))
  pw_count <- matrix(0,length(pw_names),2,dimnames=list(pw_names,c('CDp','CDn')))
  coef_p_idx <- as.integer(names(coef_sig[coef_sig > 0]))
  coef_n_idx <- as.integer(names(coef_sig[coef_sig < 0]))
  kegg_risk_rownames <- rownames(kegg_counts)
  for (i in 1:nrow(kegg_counts)){
    pw_update <- kegg_metadata[[kegg_risk_rownames[i]]]
    pw_count[pw_update,1] <- pw_count[pw_update,1] + sum(kegg_counts[i,coef_p_idx])
    pw_count[pw_update,2] <- pw_count[pw_update,2] + sum(kegg_counts[i,coef_n_idx])
  }
  pw_count <- t(t(pw_count)/colSums(pw_count))
  return(pw_count)
}

plot_pw <- function(pw_name,coefp,coefn,kegg_metadata_risk,kegg_risk_ra,other=FALSE){
  K <- ncol(kegg_risk_ra)
  coefp_idx <- as.integer(names(coefp))
  coefn_idx <- as.integer(names(coefn))
  other_idx <- (1:K)[!(1:K %in% c(coefp_idx,coefn_idx))]
  kegg_risk_pw_subset <- names(kegg_metadata_risk)[sapply(kegg_metadata_risk,function(x) pw_name %in% x)]
  kegg_risk_pw_subset <- kegg_risk_pw_subset[kegg_risk_pw_subset %in% rownames(kegg_risk_ra)]
  plot_df <- data.frame()
  for (i in coefp_idx){
    df_temp <- data.frame(Abundance=kegg_risk_ra[kegg_risk_pw_subset,i],
                          Background=rowMeans(kegg_risk_ra[kegg_risk_pw_subset,-i]),
                          Topic=paste0('T',i),
                          Diagnosis='CDp')
    plot_df <- rbind(plot_df,df_temp)
  }
  for (i in coefn_idx){
    df_temp <- data.frame(Abundance=kegg_risk_ra[kegg_risk_pw_subset,i],
                          Background=rowMeans(kegg_risk_ra[kegg_risk_pw_subset,-i]),
                          Topic=paste0('T',i),
                          Diagnosis='CDn')
    plot_df <- rbind(plot_df,df_temp)
  }
  if (other){
    for (i in other_idx){
      df_temp <- data.frame(Abundance=kegg_risk_ra[kegg_risk_pw_subset,i],
                            Background=rowMeans(kegg_risk_ra[kegg_risk_pw_subset,-i]),
                            Topic='Other',
                            Diagnosis='Unassociated')
      plot_df <- rbind(plot_df,df_temp)
    }
  }
  pp <- ggplot(plot_df,aes(x=Background,y=Abundance,colour=Topic)) + 
    geom_point(size=4,alpha=.5) + geom_abline(slope=1) +
    facet_wrap(~Diagnosis) +
    ggtitle(pw_name)
  print(pp)
}

load_fits <- function(load_file_name){
  loaded <- readRDS(load_file_name)
  fit <<- loaded$fit
  fit_frozen <<- loaded$fit_frozen
}

plot_sig_pw <- function(pw,metadata,topic_max){
  metadata <- metadata[names(topic_max)]
  pw_idx <- which(sapply(metadata, function(x) pw %in% x))
  
  y1 <- asin(sqrt(topic_max))
  x1 <- 1:length(y1)
  
  x2 <- which(names(topic_max) %in% names(pw_idx))
  y2 <- y1[x2]
  
  print(qplot(x1,y1,colour=I('gray')) + geom_point(aes(x2,y2,colour=I('red'))))
}

pw_test <- function(pw_names,topic_max,metadata,alt='greater',min_thres,model='wilcox'){
  
  cat('\n')
  record <- data.frame()
  metadata <- metadata[names(topic_max)]
  
  if (!(model %in% c('t.test','wilcox'))) stop('Must choose t.test or wilcox.')
  
  for (i in 1:length(pw_names)){ 
    
    pw <- pw_names[i]
    pw_idx <- which(sapply(metadata , function(x) pw %in% x))
    
    if (length(pw_idx) < min_thres) next # decide which number to use here
    
    mu1 <- log10(topic_max[pw_idx] + 10-4)
    mu2 <- log10(topic_max[-pw_idx] + 10-4)    
    
    if(model == 't.test'){
      mod <- t.test(mu1,mu2,alternative=alt)
    }else{
      mod <- wilcox.test(mu1,mu2,alternative=alt)
    }
    
    record <- rbind(record,data.frame(pw=pw,stat=mod$statistic,p=mod$p.value))
    
    cat(i,' ',sep='')
  }
  
  return(record)
  
}


pw_analysis <- function(kegg_risk,kegg_metadata_risk,coef_sig,max_cores=60,test='wilcox',min_thres=5,alpha=.05){
   
   require(snow)
   require(doSNOW)
   
   pw_names <- unique(unlist(kegg_metadata_risk))
   
   ncores <- ifelse(K <= max_cores, K, max_cores)
   cl <- makeCluster(ncores,type='SOCK')
   clusterExport(cl,c('pw_test'))
   registerDoSNOW(cl)
   record_matrix <- foreach(i=1:K,
                            .verbose=FALSE,.errorhandling="pass",
                            .packages=c("tidyr","dplyr"),
                            .combine='rbind') %dopar% {
      topic_max <- kegg_risk[,i]
      pw_out <- pw_test(pw_names,topic_max,kegg_metadata_risk,'greater',min_thres,test)
      pw_out$topic <- i
      
      pw_out
   }
   if(exists('cl')) stopCluster(cl)

   record_matrix$padj <- p.adjust(record_matrix$p,'BH')
   record_matrix$sig <- as.integer(record_matrix$padj < alpha)

   sig_matrix_topic <- record_matrix %>%
      mutate(sig=ifelse(sig==1,'s1','s0')) %>%
      group_by(topic,sig) %>%
      summarise(n_sig=n()) %>%
      ungroup %>%
      spread(sig,n_sig) %>%
      mutate(prop=s1/(s1+s0))
   
   sig_matrix_pw <- record_matrix %>%
      mutate(sig=ifelse(sig==1,'s1','s0')) %>%
      group_by(pw,sig) %>%
      summarise(n_sig=n()) %>%
      ungroup %>%
      spread(sig,n_sig,fill=0) %>%
      mutate(prop=s1/(s1+s0))
   
   sig_matrix_topic_pw <- record_matrix %>%
      dplyr::select(-p,-padj,-stat) %>%
      spread(pw,sig) %>%
      dplyr::select(-topic)

   p1 <- sig_matrix_topic_pw %>%
            mutate(topic = row_number(),
                   coef = coef_k,
                   corr = ifelse(coef_k>0,'CD+','CD-')) %>%
            gather(pathway,sig,-topic,-coef,-corr) %>%
            group_by(pathway) %>%
            mutate(total = sum(sig),
                   sig = factor(sig,level=c('0','1')),
                   topic = factor(topic,levels=topic[order(abs(coef_k))])) %>%
            ungroup() %>%
            filter(total > 0) %>%
            ggplot(aes(y=pathway,x=topic,fill=sig)) + geom_tile(colour='white',size=.1) +
            facet_wrap(~corr,scales='free_x') +
            scale_fill_manual(values=c('black','orange')) +
            theme(axis.text.x=element_text(angle=0,hjust=1))
   
   print(p1)
   
   return(list(topic=sig_matrix_topic,
               pw=sig_matrix_pw,
               topic_pw=sig_matrix_topic_pw,
               figure=p1))
}


remove_duplicates <- function(scores,nc){
   pairs <- as.matrix(scores[,1:2])
   pairs <- mclapply(1:NROW(pairs),function(i) as.vector(sort(pairs[i,])),mc.cores=nc)
   return(scores[!duplicated(pairs),])
}

correct_combined_score <- function(scores,prior=.063,omit=NULL){
   
   if (is.null(omit)) return(scores)
   
   score <- scores[,3:9]
   cat('1.')
   score[,omit] <- 0
   cat('2.')
   score <- score*.001
   cat('3.')
   score <- score - prior
   cat('4.')
   score[score < 0] <- 0
   cat('5.')
   score <- score/(1-prior)
   cat('6.')
   score[score > 1] <- 1
   cat('7.')
   score <- (1-apply(1-score,1,prod))*(1-prior) + prior
   cat('8.')
   score[score < 0] <- 0
   cat('9.')
   score[score > 1] <- 1
   cat('10.')
   score <- round(score*1000,0)
   cat('11.')

   scores[,10] <- score
   cat('12\n')
   
   return(scores)
}

create_ref_set <- function(prior=.063,cutoff=700,omit=NULL){
  
  scores <- readr::read_delim('~/StringDB/COG.links.detailed.v10.txt.gz',delim=' ')
  #prior <- local_prior_pp <- 0.063 ## note: when Michael recalculated this prior, it went to 0.052. This will, however, not be changed for STRING 8.1
  cat('\nCOG scores loaded (',nrow(scores), ' rows).\n',sep='')
  
  t1 <- Sys.time()
  scores <- remove_duplicates(scores,60)
  t2 <- Sys.time()
  cat('Removed duplicate rows due to order of pairs (',nrow(scores),' rows remaining) (time elapsed = ',round(t2-t1,1),'.\n',sep='')
  
  cat('Correcting scores (12 steps).\n')
  t1 <- Sys.time()
  scores <- correct_combined_score(scores,prior,omit=omit)
  t2 <- Sys.time()
  cat('Combined scores corrected (time elapsed = ',round(t2-t1,1),' m)\n',sep='')
  
  scores <- scores[scores$combined_score >= 700,]
  cat('Removed pairs with combined score less than',cutoff,'\n')
  
  return(scores)
}

create_dict <- function(R){
   cogs <- R[,1:2]
   cogs_group2_idx <- which(!(cogs$group2 %in% cogs$group1))
   cogs_group2 <- cogs[cogs_group2_idx,2:1]
   cogs <- rbind(cogs,cogs_group2)
   cogs <- as.matrix(dplyr::arrange(R[,1:2],group1))
   
   cat('Creating dictionary.\n')
   dict12 <- list()
   for (r in 1:NROW(R)){
      cog1 <- cogs[r,1]
      cog2 <- cogs[r,2]
      dict12[[cog1]] <- c(dict12[[cog1]],cog2)
   }
   cat('Dictionary completed.\n')
   
   return(dict12)
}

prep_string_data <- function(biom,ref_set=NULL){
  
  cogs_ra <- as.matrix(biom_data(biom_file))
  cogs_ra <- t(t(cogs_ra)/colSums(cogs_ra))
  cogs <- rownames(cogs_ra)
  
  return(list(cogs_ra=cogs_ra,cogs=cogs))  
  
}

find_matching_pairs <- function(dict,topic_matrix,topic,min_cutoff=.001){
   
   topic_vec <- topic_matrix[,topic]
   topic_vec <- topic_vec[topic_vec > min_cutoff]
   
   topic_comb <- t(apply(combn(names(topic_vec),2),2,sort)) ### needn't filter for duplicates because combm, not permutations
   
   cat('Starting topic ',topic,' (',NROW(topic_comb),' pairs)\n',sep='')
   
   count_direct <- 0
   count_indirect <- 0
   for (r in 1:NROW(topic_comb)){
      cog1 <- topic_comb[r,1]
      cog2 <- topic_comb[r,2]
      cog1_cxns <- dict[[cog1]]
      if (cog2 %in% cog1_cxns){
         count_direct <- count_direct + 1
      }else{
         cog2_cxns <- dict[[cog2]]
         if (length(intersect(cog1_cxns,cog2_cxns)) > 0){
            count_indirect <- count_indirect + 1
         }
      }
   }
   
   counts <- c(count_direct + count_indirect,count_direct,count_indirect)
   names(counts) <- c('Total','Direct','Indirect')
   
   return(counts)
   
}

plot_taxa_heatmap <- function(beta_frozen_ra,beta_otu,beta_meta,coef_k,taxon=6,dist='jaccard',clust='ward.D2'){
   
   require(vegan)
   
   taxon_labels <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
   beta_collapse_1 <- log10(floor(collapse_beta(beta_otu,beta_frozen_ra,beta_meta,NULL,taxon,Inf)*10^6) + 1)
   
   taxa_names <- stringr::str_extract(rownames(beta_collapse_1),'([^__]*)$')
   rownames(beta_collapse_1) <- taxa_names
   
   dismat1 <- vegdist(beta_collapse_1,method=dist)
   dismat2 <- vegdist(t(beta_collapse_1),method=dist)
   
   rc <- hclust(dismat1,clust)
   cc <- hclust(dismat2,clust)
   colors <- colorRampPalette(c("white","lightblue", "darkblue"))(n=200)
   colors_rs <- colorRampPalette(c("darkgreen","lightgreen","white","lightblue", "darkblue"))(n=200)
   
   
   coef_clus <- coef_k
   names(coef_clus) <- paste0('T',1:K)
   coef_clus_order <- order(coef_clus)
   
   
   col_left <- rev(colorRampPalette(c('lightgray','blue'))(sum(coef_clus < 0)))
   col_mid <- colorRampPalette(c('gray'))(sum(coef_clus == 0))
   col_right <- rev(colorRampPalette(c('red','lightgray'))(sum(coef_clus > 0)))
   colcolors <- c(col_left,col_mid,col_right)[order(coef_clus_order)]
   
   gplots::heatmap.2(as.matrix(beta_collapse_1),
                     scale="none",
                     density.info="histogram",
                     denscol='black',
                     labCol=cc$labels,
                     labRow=rc$labels,
                     dendrogram = 'both',
                     Colv=as.dendrogram(cc),
                     Rowv=as.dendrogram(rc),
                     col=colors,
                     ColSideColors=colcolors,
                     trace="none",
                     tracecol='black',
                     key=TRUE,keysize=2,
                     xlab='',ylab='',
                     main=taxon_labels[taxon],
                     #margins=c(10,10),cexRow=1.8,cexCol=1.8,
                     margins=c(4,20)
   )
   
}

topic_select <- function(coef_k,beta_collapse,n=48){
   
   sel <- apply(beta_collapse,1,which.max)
   if (length(sel) >= n) return(sel[1:n])
   
   names(coef_k) <- 1:K
   sel <- unique(c(sel,as.integer(names(sort(abs(coef_k[coef_k != 0]),decreasing=TRUE)))))
   if (length(sel) >= n) return(sel[1:n])
   
   sel <- c(sel,sample(which(!(1:K %in% sel)),n-length(sel),replace=FALSE))
   return(sel)
}


plot_pw_heatmap <- function(kegg_risk_ra,kegg_metadata_risk,coef_k,dist='jaccard',clust='ward.D2',var_thres=.5){
   
   require(vegan)
   
   pws <- table(unlist(sapply(kegg_metadata_risk, identity)))
   
   ra_log <- floor(kegg_risk_ra * 10^6) + 1
   ra_log_pws <- sapply(kegg_metadata_risk[rownames(ra_log)], function(x) names(which.max(pws[x])))
   ra_log <- data.frame(ra_log,pw=ra_log_pws) %>%
      group_by(pw) %>%
      summarise_each(funs(sum)) 
   rnames <- as.character(ra_log$pw)
   ra_log <- ra_log %>%
      dplyr::select(-pw) 
   ra_log <- as.matrix(log10(ra_log))
   rownames(ra_log) <- rnames
   
   rowvar <- apply(ra_log,1,var)
   ra_log <- ra_log[rowvar > quantile(rowvar,var_thres),]
   
   
   dismat1 <- vegdist(ra_log,method=dist)
   dismat2 <- vegdist(t(ra_log),method=dist)
   
   rc <- hclust(dismat1,clust)
   cc <- hclust(dismat2,clust)
   colors <- colorRampPalette(c("white","lightblue", "darkblue"))(n=200)
   colors_rs <- colorRampPalette(c("darkgreen","lightgreen","white","lightblue", "darkblue"))(n=200)
   
   
   coef_clus <- coef_k
   names(coef_clus) <- paste0('T',1:K)
   coef_clus_order <- order(coef_clus)
   
   
   col_left <- rev(colorRampPalette(c('lightgray','blue'))(sum(coef_clus < 0)))
   col_mid <- colorRampPalette(c('gray'))(sum(coef_clus == 0))
   col_right <- rev(colorRampPalette(c('red','lightgray'))(sum(coef_clus > 0)))
   colcolors <- c(col_left,col_mid,col_right)[order(coef_clus_order)]
   
   gplots::heatmap.2(as.matrix(ra_log),
                    main="",
                    scale="none",
                    density.info="histogram",
                    denscol='black',
                    labCol=cc$labels,
                    labRow=rc$labels,
                    dendrogram = 'both',
                    Colv=as.dendrogram(cc),
                    Rowv=as.dendrogram(rc),
                    col=colors,
                    ColSideColors=colcolors,
                    trace="none",
                    tracecol='black',
                    key=TRUE,keysize=2,
                    xlab='',ylab='',
                    #margins=c(10,10),cexRow=1.8,cexCol=1.8,
                    margins=c(4,20)
   )
   
   gplots::heatmap.2(as.matrix(ra_log),
                    main="Row Scaled",
                    scale="row",
                    density.info="histogram",
                    denscol='black',
                    labCol=cc$labels,
                    labRow=rc$labels,
                    dendrogram = 'both',
                    Colv=as.dendrogram(cc),
                    Rowv=as.dendrogram(rc),
                    col=colors_rs,
                    ColSideColors=colcolors,
                    trace="none",
                    tracecol='black',
                    key=TRUE,keysize=2,
                    xlab='',ylab='',
                    #margins=c(10,10),cexRow=1.8,cexCol=1.8,
                    margins=c(4,20)
   )
   
}

pw_counter <- function(kegg_risk_ra,kegg_metadata_risk){
   
   kos_metadata <- kegg_metadata_risk[rownames(kegg_risk_ra)]
   pws_complete <- unique(unlist(sapply(kos_metadata, identity)))
   pws_list <- lapply(1:length(pws_complete),function(x) NULL)
   names(pws_list) <- pws_complete
   
   for (ko in rownames(kegg_risk_ra)){
      pws <- kos_metadata[[ko]]
      for (pw in pws){
         pws_list[[pw]] <- c(pws_list[[pw]],ko)
      }
   }
   
   pws_list
   
}


test_topic <- function(pseudo_kos,k,exact=TRUE){
   target_topic <- pseudo_kos[,k]
   others_topic <- rowMeans(pseudo_kos[,-k])
   #others_topic <- apply(pseudo_kos[,-k],1,function(x) sqrt(sum(x^2)))
   #others_topic <- apply(pseudo_kos[,-k],1,median)
   
   wilcox.test(target_topic,others_topic,paired=TRUE,alternative='greater',exact=exact)$p.value
}

pw_test <- function(pw_list,kegg_risk_ra,K,pseudo_count=10^-8,fdr=TRUE,exact=TRUE,alpha=.05){
   
   pvals_list <- list()
   
   for (pw in names(pw_list)){
      kos <- pw_list[[pw]] 
      if (length(kos) < 5) next      
      topic_kos <- kegg_risk_ra[kos,]
      
      pseudo_kos <- topic_kos + pseudo_count
      pvals <- sapply(1:K,function(k) test_topic(pseudo_kos,k,exact=exact))
      
      pvals_list[[pw]] <- pvals
   }
   
   if (fdr){
      n_tests <- length(unlist(pvals_list))
      pvals_list <- lapply(pvals_list,function(x) p.adjust(x,'BH',n=n_tests))
   }
   
   pvals_list <- lapply(pvals_list,function(x) which(x < alpha))
   pvals_list <- pvals_list[sapply(pvals_list,length) > 0]
   
   return(pvals_list)
   
}

plot_corr_network <- function(fit,coef_k,thres=.1,layout=igraph::layout.auto,...){
   require(igraph)
   
   simple <- topicCorr(fit,'simple',thres)
   huge <- topicCorr(fit,'huge')
   
   cols <- ifelse(coef_k > 0, 'indianred1',
                  ifelse(coef_k < 0, 'lightblue',
                         'gray'))
   
   sig_topics <- unique(unlist(apply(simple$poscor,1,function(x) which(x > thres & x != 1))))
   
   if (length(sig_topics) > 0) print(plot(simple,vertex.color=cols[sig_topics],topics=sig_topics,vlabels=sig_topics,main='',layout=layout,...))
   
   print(plot(huge,vertex.color=cols,vlabels=1:nrow(huge$posadj),main='',layout=layout,...))
   
}

beta_prep <- function(train_fit,test_fit,counts,otu_taxa_xavier,vocab,save_fit,save_dir,save_fit_foldername,dupnum,seed_permuted,beta_type=''){
   
   beta_frozen <- t(exp(train_fit$beta$logbeta[[1]]))
   K <- ncol(beta_frozen)
   colnames(beta_frozen) <- paste0('T',1:K)
   rownames(beta_frozen) <- counts$ids[vocab,'long']
   beta_frozen_ra <<- beta_frozen
   beta_frozen <<- floor(beta_frozen_ra*max(counts$table_clean))
   beta_meta <<- data.frame(Topic=colnames(beta_frozen))
   beta_otu <<- otu_taxa_xavier[rownames(beta_frozen),]
   
   if (save_fit){
      beta_biom <- make_biom(beta_frozen,beta_meta,beta_otu)
      beta_filename_short <- paste0('beta_table_',beta_type)
      beta_filename <- paste0(beta_filename_short,'.biom')
      if (file.exists(file.path(save_dir,save_fit_foldername,beta_filename))){
         beta_filename <- paste0('dup',dupnum,'_',beta_filename)
      }
      write_biom(beta_biom,file.path(save_dir,save_fit_foldername,beta_filename))
      
      #### permute for null model ####
      set.seed(seed_permuted)
      beta_frozen_permuted <- t(apply(beta_frozen,1,function(x) x[sample(1:K,K,replace=FALSE)]))
      colnames(beta_frozen_permuted) <- colnames(beta_frozen)
      beta_biom <- make_biom(beta_frozen_permuted,beta_meta,beta_otu)
      beta_filename <- paste0(beta_filename_short,'_permuted.biom')
      if (file.exists(file.path(save_dir,save_fit_foldername,beta_filename))){
         beta_filename <- paste0('dup',dupnum,'_',beta_filename)
      }
      write_biom(beta_biom,file.path(save_dir,save_fit_foldername,beta_filename))
   }
   
}

eval_labels <- function(save_fit,load_fit,train_Zbar,test_Zbar,train_meta,test_meta,
                        save_dir,save_fit_foldername,save_coef_filename,beta_type='',
                        nc=60){
   
   save_coef_filename <- paste0(save_fit_foldername,'_coefs_',beta_type,'.rds')
   
   if (load_fit){
      
      out <- readRDS(file.path(save_dir,save_fit_foldername,save_coef_filename))
      
      return(out)
   
   }
   
   dx <- eval_diagnosis_prediction(train_Zbar,test_Zbar,train_meta,test_meta,nc) 
   is <- eval_isolation_source_prediction(train_Zbar,test_Zbar,train_meta,test_meta,nc)
   out <- list(dx=dx,is=is)
   
   if (save_fit){
      
      if (file.exists(file.path(save_dir,save_fit_foldername,save_coef_filename))){
         save_coef_filename <- paste0('dup',dupnum,'_',save_coef_filename)
      }
      saveRDS(out,file.path(save_dir,save_fit_foldername,save_coef_filename))
      
   }
   
   return(out)
   
}

ctm_frozen2 <- function(fit,documents,vocabulary,seed=NULL,max.em.its=500,emtol=1e-5,
                        avg_iters=15,verbose=TRUE,reportevery=5,
                        data=data,covariate=NULL,
                        parallel=TRUE,nc=60){
  
  #beta_frozen <- list(exp(fit$beta$logbeta[[1]]))
  #beta_frozen <- lapply(fit$beta$logbeta,exp) # i moved this below after stm.init
  
  K <- fit$settings$dim$K
  N <- length(documents)
  
  
  if (is.null(covariate)){
    if (!is.null(fit$settings$call$content)) content <- TRUE else content <- FALSE
  }else if (covariate != 0){
    content <- FALSE
    fit$beta$logbeta <- list(fit$beta$logbeta[[covariate]])
  }else{
    content <- FALSE
    fit$beta$logbeta <- list(do.call('rbind',fit$beta$logbeta))
    K <- nrow(list(do.call('rbind',fit$beta$logbeta))[[1]])
  }
  
  
  prevalence <- NULL
  init.type <- 'Spectral' # if I also set K=0, then estimates topics via Lee and Mimno 2014
  
  interactions <- TRUE # only relevant for covariates
  ngroups <- 1
  model <- NULL # allows one to restart an existing model
  sigma.prior <- 0 # strength of regularization for diagnolized covariance matrix
  gamma.prior <- c("Pooled", "L1")
  kappa.prior <- c("L1", "Jeffreys")
  control <- list()
  
  wcountvec <- unlist(lapply(documents, function(x) rep(x[1,], times=x[2,])),use.names=FALSE)
  wcounts <- list(Group.1=sort(unique(wcountvec)))
  V <- length(wcounts$Group.1)
  wcounts$x <- tabulate(wcountvec)
  rm(wcountvec)
  
  xmat <- NULL
  
  
  if (content){
    termobj <- terms(eval(fit$settings$call$content),data=data)
    if(is.null(data)) stop('Make sure you set the data used for fit.')
    char <- rownames(attr(termobj, "factors"))[1]
    yvar <- as.factor(data[[char]])
    yvarlevels <- levels(yvar)
    betaindex <- as.numeric(yvar)
  }else{
    yvarlevels <- NULL
    betaindex <- rep(1, length(documents))
  }
  
  
  
  A <- length(unique(betaindex)) #define the number of aspects
  
  settings <- list(dim=list(K=K, A=A, 
                            V=V, N=N, wcounts=wcounts),
                   verbose=verbose,
                   topicreportevery=reportevery,
                   convergence=list(max.em.its=max.em.its, em.converge.thresh=emtol),
                   covariates=list(X=xmat, betaindex=betaindex, yvarlevels=yvarlevels),
                   gamma=list(mode="CTM", prior=NULL, enet=1), # gamma is for prev so ignore match.arg(gamma.prior)
                   sigma=list(prior=sigma.prior),
                   kappa=list(LDAbeta=TRUE, interactions=FALSE, 
                              fixedintercept=TRUE, mstep=list(tol=.001, maxit=3),
                              contrast=FALSE),
                   tau=list(mode=NULL, tol=1e-5,
                            enet=1,nlambda=250, lambda.min.ratio=.001, ic.k=2,
                            maxit=1e4),
                   init=list(mode=init.type, nits=50, burnin=25, alpha=(50/K), eta=.01,
                             s=.05, p=3000, d.group.size=2000), 
                   seed=seed,
                   ngroups=ngroups)
  
  if (content){
    settings$kappa$LDAbeta <- FALSE #can't do LDA topics with a covariate 
  }else{
    settings$kappa$interactions <- FALSE #can't have interactions without a covariate.
  }
  
  
  legalargs <-  c("tau.maxit", "tau.tol", 
                  "fixedintercept","kappa.mstepmaxit", "kappa.msteptol", 
                  "kappa.enet", "nlambda", "lambda.min.ratio", "ic.k", "gamma.enet",
                  "nits", "burnin", "alpha", "eta", "contrast",
                  "rp.s", "rp.p", "rp.d.group.size", "SpectralRP")
  
  # Process the Seed
  if(is.null(settings$seed)) {
    #if there is no seed, choose one and set it, recording for later
    seed <- floor(runif(1)*1e7) 
    set.seed(seed)
    settings$seed <- seed
  } else {
    #otherwise just use the provided seed.
    set.seed(settings$seed)
  }
  
  
  globaltime <- proc.time()  # maybe remove?
  
  settings$init$mode <- 'LDA'
  
  #initialize
  model <- stm.init(documents, settings)
  beta <- list(beta=lapply(fit$beta$logbeta,exp))
  beta$kappa <- fit$beta$kappa
  
  #unpack
  mu <- list(mu=model$mu)
  sigma <- model$sigma
  
  
  lambda <- model$lambda
  convergence <- NULL 
  #discard the old object
  rm(model)
  
  #Pull out some book keeping elements
  ntokens <- sum(settings$dim$wcounts$x)
  betaindex <- settings$covariates$betaindex # since no covariates, intercept? therefore 1?
  stopits <- FALSE
  suffstats <- vector(mode="list", length=ngroups)
  
  ############
  #Step 2: Run EM
  ############
  while(!stopits) {
    
    #####
    # Non-Blocked Updates
    #####
    
    t1 <- proc.time()
    
    #run the model
    suffstats <- stm:::estep(documents=documents, beta.index=betaindex, 
                             update.mu=(!is.null(mu$gamma)),  
                             beta$beta, lambda, mu$mu, sigma, 
                             verbose)
    
    t1 <- proc.time()
    
    Z <- suffstats$Z
    sigma.ss <- suffstats$sigma
    lambda <- suffstats$lambda
    bound.ss <- suffstats$bound
    
    #do the m-step
    mu <- opt.mu(lambda=lambda, mode=settings$gamma$mode, 
                 covar=settings$covariates$X, settings$gamma$enet)
    sigma <- opt.sigma(nu=sigma.ss, lambda=lambda, 
                       mu=mu$mu, sigprior=settings$sigma$prior)
    
    
    #Convergence
    convergence <- convergence.check(bound.ss, convergence, settings)
    stopits <- convergence$stopits
    
    cat('converged:',convergence$converged,'\n')
    
  }
  
  #######
  #Step 3: Construct Output
  #######
  time <- (proc.time() - globaltime)[3]
  #convert the beta back to log-space
  beta$logbeta <- beta$beta
  for(i in 1:length(beta$logbeta)) {
    beta$logbeta[[i]] <- log(beta$logbeta[[i]])
  }
  beta$beta <- NULL
  lambda <- cbind(lambda,0)
  model <- list(mu=mu, sigma=sigma, beta=beta, 
                settings=settings,
                vocabulary=vocabulary, convergence=convergence, 
                theta=exp(lambda - row.lse(lambda)), 
                eta=lambda[,-ncol(lambda), drop=FALSE],
                invsigma=solve(sigma), time=time, version=utils::packageDescription("stm")$Version)
  class(model) <- "STM"  
  
  return(model)
}
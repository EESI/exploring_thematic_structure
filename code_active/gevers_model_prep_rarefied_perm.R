#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) stop("At least one argument must be supplied (input file).n", call.=FALSE)

library(stm)
library(biom)
library(readr)
library(tidyr)
library(dplyr)
library(fastICA)
library(randomForest)
library(stringr)
library(kernlab)
library(Rcpp)
library(parallel)
library(foreach)
library(ape)
library(phyloseq)
library(doParallel)
library(stm)
library(LDAvis)
library(caret)
library(glmnet)
library(ggplot2)
library(knitr)
library(gridExtra)


source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

permTest <- function (formula, stmobj, treatment, nruns = 100, documents, 
                      vocab, data, seed = NULL, stmverbose = TRUE, uncertainty = "Global") 
{
  if (!requireNamespace("clue", quietly = TRUE)) 
    stop("Install the clue package to use this function")
  settings <- stmobj$settings
  if (!is.data.frame(data)) 
    stop("data object must be a data.frame containing treatment variable.")
  if (!(treatment %in% colnames(data))) 
    stop("treatment variable must be in the data.frame")
  if (!all(data[, treatment] %in% c(0, 1))) 
    stop("treatment variable must be binary and coded 0 and 1.")
  prob <- mean(data[, treatment])
  if (is.null(seed)) {
    seed <- floor(runif(1) * 10000000)
    set.seed(seed)
  }
  else {
    set.seed(seed)
  }
  settings$verbose <- stmverbose
  settings$keepHistory <- FALSE
  betaref <- exp(stmobj$beta$logbeta[[1]])
  qeffects <- function(formula, mod, data, uncertainty) {
    prep <- stm::estimateEffect(formula, stmobj = mod, metadata = data, 
                                uncertainty = uncertainty)
    betas <- lapply(prep$parameters, function(x) lapply(x, 
                                                        function(y) stm:::rmvnorm(100, y$est, y$vcov)))
    simbetas <- lapply(betas, do.call, what = rbind)
    parcol <- which(colnames(settings$covariates$X) == treatment)
    par <- lapply(simbetas, function(x) quantile(x[, parcol], 
                                                 c(0.025, 0.5, 0.975)))
    par
  }
  cat("Calculating effects for reference model \n")
  ref.effect <- qeffects(formula, stmobj, data, uncertainty)
  tosave <- list()
  for (i in 1:(nruns - 1)) {
    settings$seed <- floor(runif(1) * 10000000)
    data[, treatment] <- rbinom(n = nrow(data), size = 1, 
                                prob = prob)
    termobj <- terms(formula, data = data)
    if (attr(termobj, "response") == 1) 
      stop("Response variables should not be included in prevalence formula.")
    settings$covariates$X <- model.matrix(termobj, data = data)
    cat(sprintf("Running model %i of %i \n", (i + 1), (nruns)))
    mod <- stm:::stm.control(documents, vocab, settings = settings, model=NULL)
    par <- qeffects(formula, mod, data, uncertainty)
    betamod <- exp(mod$beta$logbeta[[1]])
    align <- clue::solve_LSAP(betaref %*% t(betamod), maximum = TRUE)
    tosave[[i]] <- par[align]
  }
  out <- list(ref = ref.effect, permute = tosave, variable = treatment, 
              seed = seed)
  class(out) <- "STMpermute"
  return(out)
}


K <- as.integer(args[1])
f <- as.integer(args[2]) ### which fit?

#K <- 25
m <- 5
#f <- 4
rare_min <- 1000

dir_name <- sprintf('~/Dropbox/stm_microbiome/data_active/Gevers_ti/stm_s97_rarefied_%s_supervised/models/',rare_min)
fit_filename <- sprintf('model_K_%s_nfits_%s.rds',K,m)

dat <- readRDS(file.path(dir_name,fit_filename))



fit <- dat$fits[[paste0('fit',f)]]
form <- as.formula(fit$settings$call$prevalence)
docs <- dat$docs
vocab <- dat$vocab
meta <- dat$meta


cat('\n',
    'm = ',m,'\n',
    'K = ',K,'\n',
    'formula = ',paste0(form,collapse=''),'\n\n',sep='')

seed_perm <- 89

perm <- permTest(formula=~DIAGNOSIS, fit, 'DIAGNOSIS',
                        nruns=25,
                        documents=docs, vocab=vocab, data=meta, seed=seed_perm,
                        uncertainty="Global",
                        stmverbose=TRUE)


dir_perm_name <- gsub('models','perm',dir_name)
dir.create(dir_perm_name,showWarnings=FALSE,recursive=TRUE)
perm_filename <- sprintf('perm_K_%s_fit%s.rds',K,f)
saveRDS(list(dat=dat,
             K=K,
             m=m,
             perm=perm,
             seed_perm=seed_perm),
        file.path(dir_perm_name,perm_filename))

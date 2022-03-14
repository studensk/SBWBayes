library(tmbstan)
library(TMB)
library(tidyverse)
library(lme4)
library(ellipse)
library(MASS)
library(Matrix)
library(splines)
library(parallel)
all.days.df <- read.csv('data/all_days_df.csv')
source('code/functions.R')
source('code/test_stan.R')
nList <- lme4:::namedList
on.data.grouped <- read.csv('data/grouped.csv')

## compile & load TMB

ptm.start <- proc.time()
cl <- makeCluster(10)
clusterExport(cl, varlist = c('on.data.grouped', 'nList'))
clusterEvalQ(cl,{
  library(tmbstan)
  library(TMB)
  library(tidyverse)
  library(lme4)
  library(ellipse)
  library(MASS)
  library(Matrix)
  library(splines)
  library(parallel)
})

group.lst <- parLapply(cl, 1:10, function(grp) {
  all.days.df <- read.csv('data/all_days_df.csv')
  source('code/functions.R')
  source('code/test_stan.R')
  nList <- lme4:::namedList
  on.data.grouped <- read.csv('data/grouped.csv')
  
  basename <- "stan_mle_structured_parametric"
  cc <- compile(paste0('code/', basename, '.cpp'))
  setwd("code")
  try(dyn.unload(dynlib(basename)),silent=TRUE)
  dyn.load(dynlib(basename))
  setwd('..')
  options(mc.cores = parallel::detectCores())
  chains <- 1
  
  prior.samp <- function(chains) {
    l <- lapply(1:chains, function(x) {
      lst <- list(
        'phi_rho' = rbeta(1, 4, 4),
        'psi_rho' = rbeta(1, 4, 4),
        'y0_rho' = rgamma(1, 10, scale = 0.1),
        'HA' = rgamma(5, 5, scale = 0.2),
        'TL' = rnorm(5, 284, 2),
        'HL' = rgamma(5, 6, scale = 2),
        'TH' = rnorm(5, 304, 2),
        'HH' = rgamma(5, 10, scale = 3),
        's_eps' = exp(rnorm(5, -1.5, 0.1)),
        's_alpha' = exp(rnorm(1, -1.5, 0.3)),
        'alpha' = 0.05
      )
      return(lst)
    })
    return(l)
  }
  
  parms <- list(
    'phi_rho' = 0.5,
    'psi_rho' = 0.5,
    'y0_rho' = 0.5,
    'HL' = rep(8.4, 5),
    'HH' = rep(9.3, 5),
    'HA' = rep(1.4, 5), 
    'TL' = rep(285.90432, 5),
    'TH' = rep(306.47530, 5),
    'TA' = rep(298, 5),
    's_eps' = rep(0.3, 5),
    's_alpha' = 0.22)
  parms1 <- parms
  diagnostics <- list()
  
  set.seed(123)
  
  p <- 'ON'
  stages <- paste0('L', 2:6)
  
  ## Obtain MLEs and covariance estimates for each stage
  
  ## Sample 200 draws from prior
  ps <- prior.samp(100)
  mle.lst <- lapply(1:100, function(ind) {
    parms <- parms1
    dd <- subset(on.data.grouped, group_index != grp)
    dd$stagename <- dd$stage
    dd$stage <- as.numeric(factor(dd$stage)) - 1
    dd2 <- with(dd, nList(temp1,time1,temp2,time2,block,time2d,stage,
                          nobs=as.integer(nobs)))
    dd2$block <- as.numeric(factor(dd2$block)) - 1
    dd2$use_prior <- 0 
    
    nblock <- length(unique(dd2$block))
    nstage <- length(unique(dd2$stage))
    
    pnames <- c("phi_rho", "psi_rho", "y0_rho",
                "HA","TL","HL","TH","HH",
                "s_eps", "s_alpha", 'alpha')
    x <- ps[[ind]]
    x$alpha <- rep(0, nstage)
    
    ff <- MakeADFun(data=dd2,
                    parameters=as.list(x[pnames]),
                    DLL=basename,
                    random = c('alpha'),
                    silent=TRUE)
    
    ## Optimize over objective with aid of gradient function
    opt <- tryCatch({optim(par = ff$par, fn = ff$fn, gr = ff$gr,
                           method = 'BFGS', hessian = TRUE,
                           control = list("maxit" = 1000))}, 
                    error = function(e) {return("Initial value not finite")})
    if (length(opt) == 1) {return('Initial value not finite')}
    sdr <- sdreport(ff)
    opt$par <- ff$env$last.par.best
    hess <- opt$hessian
    vcov <- solve(hess)
    opt$vcov <- vcov
    
    alpha.var <- sdr$diag.cov.random
    names(alpha.var) <- paste0('alpha.', 1:5, '.')
    alpha.var.diag <- diag(alpha.var)
    opt$alpha.var <- alpha.var.diag
    
    vcov.full <- as.matrix(bdiag(vcov, alpha.var.diag))
    all.names <- c(rownames(vcov), names(alpha.var))
    colnames(vcov.full) <- all.names
    rownames(vcov.full) <- all.names
    opt$vcov.full <- vcov.full
    
    if (min(eigen(vcov)$values) < -1) {opt$value <- NA}
    
    if (length(opt) > 1) {
      if (opt$convergence != 0) {
        return("Did not converge")
      }
    }
    return(opt)
  })
  
  ## Save optimal solution with minimal objective value
  w <- which.min(sapply(mle.lst, function(x) {ifelse(length(x) > 1, x$value, NA)}))
  mle <- mle.lst[[w]]  
  hls <- mle$par[which(names(mle$par) == 'HL')]
  mle$par[which(names(mle$par) == 'HL')] <- -hls
  
  vcov <- mle$vcov.full
  eig <- eigen(vcov)
  vals <- pmax(eig$values, 0)
  vecs <- eig$vectors
  vcov.new <- vecs %*% diag(vals) %*% solve(vecs)
  mv <- as.data.frame(mvrnorm(14000, mle$par, vcov.new))
  names(mv) <- rownames(mle$vcov.full)
  
  #mv <- subset(mv, !is.na(log(-HL/HH)))
  c.names <- unlist(sapply(names(mv), function(x) {strsplit(x, split = '[.]')[[1]][1]}))
  t.names <- table(c.names)
  
  col.lst <- lapply(1:length(t.names), function(ind) {
    name <- names(t.names)[ind]
    count <- t.names[ind]
    value.df <- mv[,which(c.names == name)]
    if (count == 5) {values <- unlist(value.df)}
    else {values <- rep(unlist(value.df), 5)}
    new <- data.frame(values)
    names(new) <- name
    return(new)
  })
  col.df <- bind_cols(col.lst)
  col.df$stage <- factor(rep(stages, each = nrow(mv)))
  col.df$rho25 <- mapply(function(stage, phi, psi, y0, alpha, sa) {
    q <- quadratic(stage, phi, psi, y0)
    sa2 <- -0.5*sa*sa
    mult <- exp(alpha*sa + sa2)
    rho <- q*mult
    return(rho)
  }, as.numeric(col.df$stage)-3, col.df$phi_rho, 
  col.df$psi_rho, col.df$y0_rho, col.df$alpha, col.df$s_alpha)
  col.df$group <- grp
  return(col.df)
})
stopCluster(cl)
ptm.end <- proc.time()
ptm.tot <- ptm.end - ptm.start
group.df <- bind_rows(group.lst)
write.csv(group.df, 'data/group_df.csv', row.names = FALSE)


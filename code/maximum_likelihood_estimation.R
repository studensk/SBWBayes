library(tmbstan)
library(TMB)
library(tidyverse)
library(lme4)
library(ellipse)
library(MASS)
library(Matrix)
library(splines)
all.days.df <- read.csv('data/all_days_df.csv')
source('code/functions.R')
source('code/test_stan.R')
load("code/regniere_JIP2015_basic.RData")
nList <- lme4:::namedList

##### Independent Runs #####

## compile & load TMB
basename <- "stan_mle"
cc <- compile(paste0('code/', basename, '.cpp'))
setwd("code")
try(dyn.unload(dynlib(basename)),silent=TRUE)
dyn.load(dynlib(basename))
setwd('..')
options(mc.cores = parallel::detectCores())
chains <- 4

prior.samp <- function(chains) {
  l <- lapply(1:chains, function(x) {
    lst <- list(
      'rho25' = rgamma(1, 2, scale = 0.25),
      'HA' = rgamma(1, 5, scale = 0.2),
      'TL' = rnorm(1, 284, 2),
      'HL' = log(rgamma(1, 6, scale = 2)),
      'TH' = rnorm(1, 304, 2),
      'HH' = log(rgamma(1, 10, scale = 3)),
      's_eps' = exp(rnorm(1, -1.5, 0.1)),
      's_upsilon' = exp(rnorm(1, -2.5, 0.05))
    )
    return(lst)
  })
  return(l)
}

parms <- list('rho25' = 0.22647, 
              'HL' = 8.4,
              'HH' = 9.3,
              'HA' = 1.4,
              'TL' = 285.90432,
              'TH' = 306.47530,
              's_eps' = 0.3,
              's_upsilon' = 0.1)
parms1 <- parms
diagnostics <- list()

set.seed(123)

p <- 'ON'
stages <- paste0('L', 2:6)

## Obtain MLEs and covariance estimates for each stage
stg.lst <- list()
for (s in stages) {
  ## Sample 100 draws from prior
  ps <- prior.samp(100)
  
  ## Find MLEs using each draw as a different set of starting values
  mle.lst <- lapply(1:100, function(ind) {
    parms <- parms1
    dd <- subset(all.days.df, province == 'ON' & stage == s)
    dd <- subset(dd, generation == 'F1')
    dd2 <- with(dd, nList(temp1,time1,temp2,time2,block,time2d,
                          nobs=as.integer(nobs)))
    dd2$block <- dd2$temp1/5 - 1
    dd2$use_prior <- 0
    if (min(dd2$block) == 1) {dd2$block <- dd2$block - 1}
    parms$upsilon <- rep(0, length(unique(dd2$block)))
    pnames <- c("rho25","HA","TL","HL","TH","HH","s_eps",
                "s_upsilon",
                'upsilon')
    print(ind)
    x <- ps[[ind]]
    x <- append(x, list('upsilon' = rep(0, 7)))
    x$rho25 <- rgamma(1, shape = 31, scale = 0.007)
    ff <- MakeADFun(data=dd2,
                    parameters=as.list(x[pnames]),
                    DLL=basename,
                    random='upsilon',
                    silent=TRUE)
    
    ## Optimize over objective with aid of gradient function
    opt <- tryCatch({optim(par = ff$par, fn = ff$fn, gr = ff$gr,
                           method = 'BFGS', hessian = TRUE,
                           control = list("maxit" = 1000))}, 
                    error = function(e) {return("Initial value not finite")})
    sdr <- sdreport(ff)
    opt$par <- ff$env$last.par.best
    hess <- opt$hessian
    vcov <- solve(hess)
    opt$vcov <- vcov
    
    ups.var <- sdr$diag.cov.random
    names(ups.var) <- paste0('upsilon.', seq(5, 35, by = 5), '.')
    ups.var.diag <- diag(ups.var)
    opt$ups.var <- ups.var.diag
    
    vcov.full <- as.matrix(bdiag(vcov, ups.var.diag))
    all.names <- c(rownames(vcov), names(ups.var))
    colnames(vcov.full) <- all.names
    rownames(vcov.full) <- all.names
    opt$vcov.full <- vcov.full
    
    if (min(eigen(vcov.full)$values) < -1) {opt$value <- NA}
    
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
  stg.lst <- append(stg.lst, list(mle))
}
names(stg.lst) <- stages


agm.lst <- lapply(stages, function(st) {
  print(st)
  mle <- stg.lst[[st]]
  vcov <- mle$vcov.full
  eig <- eigen(vcov)
  vals <- pmax(eig$values, 0)
  vecs <- eig$vectors
  vcov.new <- vecs %*% diag(vals) %*% solve(vecs)
  mv <- as.data.frame(mvrnorm(100000, mle$par, vcov.new))
  names(mv) <- rownames(mle$vcov.full)
  mv$HL <- -mv$HL
  
  mv <- subset(mv, !is.na(log(-HL/HH)))
  mv$s_upsilon <- pmax(mv$s_upsilon, 0)
  
  gs.mv <- get_curves(mv, spline = TRUE)
  ag.gs.mv <- aggregate(data = gs.mv, 
                        rate ~ temp,
                        function(x) {
                          quantile(x, probs = c(0.025, 0.5, 0.975))})
  rm(gs.mv)
  rate.df <- as.data.frame(ag.gs.mv$rate)
  names(rate.df) <- c('q025', 'rate',  'q975')
  s.agm <- cbind(ag.gs.mv$temp, rate.df)
  names(s.agm)[1] <- 'temp'
  
  
  s.agm$stage <- st
  return(s.agm)
})
agm.df <- bind_rows(agm.lst)
agm.df$fit <- 'MLE'
names(agm.lst) <- stages

post.df.all.lst <- lapply(stages, function(st) {
  post <- read.csv(paste0('code/output/simulations/posterior_ON_', st, '.csv'))
  post$HL <- -post$HL
  gc <- get_curves(post, spline = TRUE)
  ag <- aggregate(data = gc, rate ~ temp,
                  function(x) {quantile(x, probs = c(0.025, 0.5, 0.975))})
  rate.df <- as.data.frame(ag$rate)
  names(rate.df) <- c('q025', 'rate',  'q975')
  agp <- cbind(ag$temp, rate.df)
  names(agp)[1] <- 'temp'
  agp$stage <- st
  return(agp)
})
post.df.all <- bind_rows(post.df.all.lst)
post.df.all$fit <- 'Bayes'
all.fits <- bind_rows(agm.df, post.df.all)
ontario <- subset(all.days.df, province == 'ON')

ontario$r.est2 <- pmin(ontario$r.est2, 2)

ggplot(data = all.fits) +
  geom_ribbon(aes(x = temp, ymin = q025, ymax = q975, fill = fit), 
              alpha = 0.4) +
  geom_line(aes(x = temp, y = rate, col = fit)) +
  geom_segment(data = ontario,
               aes(xend = temp1, x = temp1, y = r.est1, yend = r.est2, size = nobs),
               alpha = 0.5) +
  geom_point(data = subset(ontario, r.est1 == r.est2),
             aes(x = temp1, y = r.est1, size = nobs), alpha = 0.5) +
  facet_wrap(vars(stage), scales = 'free') + 
  theme_minimal() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    size = 14),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 15, l = 0),
                                    size = 14),
        strip.text = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 12)) +
  labs(x = 'Rearing Temperature',
         y = 'Development Rate',
         fill = 'Fitting Method', 
         col = 'Fitting Method',
         size = '# of Observations')
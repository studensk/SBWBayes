library(tmbstan)
library(TMB)
library(tidyverse)
library(lme4)
all.days.df <- read.csv('data/all_days_df.csv')
source('code/functions.R')
source('code/test_stan.R')
nList <- lme4:::namedList

##### Stage Structured #####
## The following code is analogous to the previous block, but with the implementation
##  of quadratic stage structure
load("code/regniere_JIP2015_basic.RData")
nList <- lme4:::namedList

basename <- "regniere_structured_parametric"
basename_loc <- paste0('code/', basename)
cc <- compile(paste0(basename_loc, '.cpp'))
try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
dyn.load(dynlib(basename_loc))
options(mc.cores = parallel::detectCores())

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
chains <- 4
parm.lst <- prior.samp(chains)

p <- 'ON'
parms <- parms1
dd <- subset(all.days.df, province == p)
dd <- subset(dd, generation == 'F1')
dd$stagename <- dd$stage
dd$stage <- as.numeric(factor(dd$stage)) - 1
dd2 <- with(dd, nList(temp1,time1,temp2,time2,block,time2d,stage,
                      nobs=as.integer(nobs)))
dd2$block <- as.numeric(factor(dd2$block)) - 1
dd2$use_prior <- 1

nblock <- length(unique(dd2$block))
nstage <- length(unique(dd2$stage))

parms$upsilon <- rep(0, nblock)
parms$alpha <- rep(0, nstage)

pnames <- c(#"rho25",
  "phi_rho", "psi_rho", "y0_rho",
  "HA","TL","HL","TH","HH","TA",
  "s_eps", "s_alpha",'alpha')
ff <- MakeADFun(data=dd2,
                parameters=as.list(parms[pnames]),
                DLL=basename,
                random=c('alpha'),
                silent=TRUE)
parm.lst <- lapply(parm.lst, function(x) {
  x$alpha <- rep(0.05, nstage)
  dat <- subset(dd, temp1 == 25 & stagename == 'L4')
  x$y0_rho <- weighted.mean(dat$r.est1, dat$nobs)
  x$HH <- pmin(x$HH, 30)
  x$TH <- pmax(x$TH, 304)
  x$HL <- pmax(x$HL, -15)
  return(x)
})
stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                 chains = chains, iter = 5000, warmup = 2000,
                 control = list(max_treedepth = 15, adapt_delta = 0.999))

stan.array <- as.array(stan1)
Rhat <- max(apply(stan.array, 3, Rhat))
ESS.bulk <- min(apply(stan.array, 3, ess_bulk))
ESS.tail <- min(apply(stan.array, 3, ess_tail))

g <- get_sampler_params(stan1, inc_warmup = FALSE)
g.df <- do.call('rbind', g)
div <- sum(g.df[,'divergent__'])
diagnostics <- append(diagnostics, list(c('Rhat' = Rhat, 
                                          'ESS.bulk' = ESS.bulk,
                                          'ESS.tail' = ESS.tail,
                                          'divergences' = div)))
post.df <- as.data.frame(stan1)

dd2.df <- as.data.frame(dd2)
stages <- c('L2', 'L3', 'L4', 'L5', 'L6')

## Clean model output for plotting and simulations
p.lst <- lapply(1:5, function(ind) {
  stg.ind <- (ind - 1) - 2
  
  alpha <- post.df[, paste0('alpha[', ind, ']')]
  sa2 <- -0.5*(post.df$s_alpha)^2
  alpha.mult <- exp(alpha*post.df$s_alpha + sa2)
  
  a <- (post.df$y0_rho/8)*(post.df$phi_rho + post.df$psi_rho - 2)
  b <- (post.df$y0_rho/4)*(post.df$psi_rho - post.df$phi_rho)
  c <- post.df$y0_rho
  rho <- a*stg.ind^2 + b*stg.ind + c
  rho25 <- rho*alpha.mult 
  q.pars <- post.df[,1:3]
  
  par.cols <- post.df[, grep(paste0('[', ind, ']'), names(post.df))[1:6]]
  pars <- as.data.frame(cbind(q.pars,
                              rho25, 
                              par.cols))
  names(pars) <- c('phi_rho', 'psi_rho', 'y0_rho', 'rho25', 
                   'HA', 'TL', 'HL', 'TH', 'HH', 's_eps')
  pars$HL <- -pars$HL
  pars$TA <- sapply(1:nrow(pars), function(r) {
    row <- pars[r,]
    ta <- with(row, ta.fun(HL, HH, TL, TH))
    return(ta)
  })
  
  stg.df <- subset(dd2.df, stage == ind - 1)
  stg.df <- stg.df[order(stg.df$temp1),]
  blocks <- unique(stg.df$block) + 1
  
  temps <- sort(unique(dd2$temp1))
  
  stage.df <- as.data.frame(cbind(pars, 'alpha' = alpha.mult, 
                                  's_alpha' = post.df$s_alpha))
  stage.df$stage <- stages[ind]
  return(stage.df)
})
p.df <- bind_rows(p.lst)

write.csv(p.df, 'data/model_results.csv', row.names = FALSE)

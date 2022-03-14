library(tmbstan)
library(TMB)
library(tidyverse)
library(lme4)
library(dclone)
all.days.df <- read.csv('data/all_days_df_ind.csv')
source('code/functions.R')
source('code/test_stan.R')
nList <- lme4:::namedList

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

## Ontario data
on.data <- subset(all.days.df, province == 'ON')
on.data <- on.data[order(on.data$cup),]
on.data <- on.data[order(on.data$stage),]

folds <- 10
set.seed(100)
group.ord <- order(runif(nrow(on.data)))
temps <- sort(unique(all.days.df$temp1))
temp.lst <- lapply(temps, function(x) {
  sub <- subset(on.data, temp1 == x)
  sub <- sub[order(sub$stage),]
  n.tot <- length(unique(sub$cup))
  ind.vec <- rep(1:folds, floor(n.tot/folds))
  leftover <- n.tot %% folds
  if (leftover > 0) {extra <- sample(1:10, leftover, replace = FALSE)}
  else {extra <- vector()}
  ind.vec <- c(ind.vec, extra)
  ind.vec.random <- sample(ind.vec, length(ind.vec), replace = FALSE)
  sub$group_index <- rep(ind.vec.random, length(unique(sub$stage)))
  return(sub)
})
on.data.grouped <- bind_rows(temp.lst)
on.data.grouped$nobs <- 1
on.data.grouped <- subset(on.data.grouped, select = -c(cup, X))
#write.csv(on.data.grouped, 'data/grouped.csv', row.names = FALSE)

ptm.start <- proc.time()
cl <- makeCluster(10)
clusterExport(cl, varlist = c('on.data.grouped', 'parms1', 'parm.lst',
                              'prior.samp', 'nList'))
clusterEvalQ(cl,{
  library('TMB')
  library('tmbstan')
  library('rstan')
})
results1 <- parLapply(cl, 1:10, function(x) {
  basename <- "regniere_structured_parametric"
  basename_loc <- paste0('code/', basename)
  cc <- compile(paste0(basename_loc, '.cpp'))
  try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
  dyn.load(dynlib(basename_loc))
  options(mc.cores = parallel::detectCores())
  dd.orig <- subset(on.data.grouped, group_index != x)
  dd <- aggregate(data = dd.orig, nobs ~ ., sum)
  parms <- parms1
  dd$stagename <- dd$stage
  dd$stage <- as.numeric(factor(dd$stage)) - 1
  dd2 <- with(dd, nList(temp1,time1,temp2,time2,block,time2d,stage,
                        nobs=as.integer(nobs)))
  dd2$block <- as.numeric(factor(dd2$block)) - 1
  dd2$use_prior <- 1
  
  nblock <- length(unique(dd2$block))
  nstage <- length(unique(dd2$stage))
  
  parms$alpha <- rep(0, nstage)
  
  pnames <- c("phi_rho", "psi_rho", "y0_rho",
              "HA","TL","HL","TH","HH","TA",
              "s_eps", "s_alpha", 'alpha')
  ff <- MakeADFun(data=dd2,
                  parameters=as.list(parms[pnames]),
                  DLL=basename,
                  random=c('alpha'),
                  silent=TRUE)
  
  parm.lst <- lapply(parm.lst, function(y) {
    y$alpha <- rep(0.05, nstage)
    dat <- subset(dd, temp1 == 25 & stagename == 'L4')
    y$y0_rho <- weighted.mean(1/(dat$time2 - 0.5), dat$nobs)
    y$HH <- pmin(y$HH, 30)
    y$TH <- pmax(y$TH, 304)
    y$HL <- pmax(y$HL, -15)
    return(y)
  })
  stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                   chains = 4, iter = 5000, warmup = 1500,
                   control = list(max_treedepth = 15, adapt_delta = 0.999))
  
  post.df <- as.data.frame(stan1)
  write.csv(post.df, paste0('data/parametric_cv', x, '.csv'))
  
  stan.array <- as.array(stan1)
  Rhat <- max(apply(stan.array, 3, Rhat))
  ESS.bulk <- min(apply(stan.array, 3, ess_bulk))
  ESS.tail <- min(apply(stan.array, 3, ess_tail))
  
  g <- get_sampler_params(stan1, inc_warmup = FALSE)
  g.df <- do.call('rbind', g)
  div <- sum(g.df[,'divergent__'])
  diagnostics <- list(c('Rhat' = Rhat, 
                        'ESS.bulk' = ESS.bulk,
                        'ESS.tail' = ESS.tail,
                        'divergences' = div))
  return(list(stan1, diagnostics))
})
stopCluster(cl)
ptm.end <- proc.time()
ptm.tot <- ptm.end - ptm.start

diag.lst <- lapply(results1, function(x) {x[[2]]})
diag.df <- bind_rows(diag.lst)
diag.df$group.index <- 1:10
write.csv(diag.df, 'data/parametric_diagnostics.csv', row.names = FALSE)


library(tidyverse)
library(tmbstan)
library(TMB)
library(lme4)
library(parallel)

source('code/data_simulation.R')
source('code/test_stan.R')

##### Generate and Clean Data #####
set.seed(123)
ptm <- proc.time()
data.lst <- gen.pops(20)
proc.time() - ptm

priors.lst <- lapply(1:length(data.lst), function(i) {
  x <- data.lst[[i]]
  dat <- x[[1]]
  dat$prior.samp <- i
  return(dat)
})
priors.df <- bind_rows(lapply(priors.lst, as.data.frame))
priors.df$stage <- rep(stages, length(priors.lst))

all.data <- bind_rows(lapply(data.lst, function(x) {x[[2]]}))
names(all.data) <- c('temp1', 'temp2', 'stage', 'time1', 'time2', 'nobs', 'prior.samp')
#all.data$block <- all.data$temp1/5
all.data$block <- paste(all.data$stage, all.data$temp1, sep = '_')
lev.ord <- paste(rep(stages, each = 7), rep(seq(5, 35, by = 5), 5), sep = '_')
all.data$block <- factor(all.data$block)
all.data$block <- factor(all.data$block, levels = lev.ord)
all.data$time2d <- pmax(0, all.data$time2 - 1)

## Remove "skippers"
# all.data <- subset(all.data, time1 > 0 | time2 > 0)

## Correction for early moulters
early <- which(all.data$time2 == 0)
early.sub <- all.data[early,]
early.sub$time2 <- early.sub$time1
early.sub$time1 <- 0
early.sub$time2d <- pmax(early.sub$time2 - 1, 0)
early.sub$temp2 <- early.sub$temp1

all.data[early,] <- early.sub
all.data$time2 <- pmax(all.data$time2, 1)

write.csv(all.data, 'data/sim_data.csv')

##### Run Model #####

## Params to initialize MakeADFun
parms <- list(
  'phi_rho' = 0.5,
  'psi_rho' = 0.5,
  'y0_rho' = 0.38,
  'HL' = 8.1,
  'HH' = 20,
  'HA' = 0.72, 
  'TL' = 285.90432,
  'TH' = 306.47530,
  's_eps' = rep(0.3, 5),
  's_upsilon' = rep(0.6, 5),
  's_alpha' = 0.22)


set.seed(124)

basename <- "regniere_structured"
basename_loc <- paste0('code/', basename)
cc <- compile(paste0(basename_loc, '.cpp'))
try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
dyn.load(dynlib(basename_loc))

nList <- lme4:::namedList
stages <- paste0('L', 2:6)

samps <- unique(all.data$prior.samp)
sv.lst <- lapply(samps, function(x) {
  print(x)
  dd <- subset(all.data, prior.samp == x)
  dd$stagename <- dd$stage
  dd$stage <- as.numeric(factor(dd$stage)) - 1
  
  test.lst <- lapply(1:8, function(chain) {
    ind <- (x-1)*8 + chain
    anyna <- TRUE
    while(anyna) {
      ps <- prior.samp(1)[[1]]
      dd2 <- with(dd, nList(temp1,time1,temp2,time2,block,time2d,stage,
                            nobs=as.integer(nobs)))
      dd2$block <- as.numeric(factor(dd2$block)) - 1
      dd2$use_prior <- 1
      
      nblock <- length(unique(dd2$block))
      nstage <- length(unique(dd2$stage))
      
      parms$upsilon <- rep(0, nblock)
      parms$alpha <- rep(0.05, nstage)
      
      pnames <- c(#"rho25",
        "phi_rho", "psi_rho", "y0_rho",
        "HA","TL","HL","TH","HH",
        "s_eps", 's_upsilon', "s_alpha",
        'upsilon', 'alpha')
      ff <- MakeADFun(data=dd2,
                      parameters=as.list(parms[pnames]),
                      DLL=basename,
                      random=c('upsilon', 'alpha'),
                      silent=TRUE)
      
      mod <- tmbstan(ff, init = unlist(parms[1:(length(parms)-2)]), silent = TRUE, 
                     chains = 0)
      
      ps$upsilon <- rep(0, nblock)
      ps$alpha <- rep(0.05, nstage)
      ps$HL <- abs(ps$HL)
      
      gr <- grad_log_prob(mod, unlist(ps))
      anyna <- any(is.na(gr))
    }
    return(ps)
  })
  
  return(test.lst)
})

## Set up data for model input
cl <- makeCluster(40)
clusterExport(cl, c('all.data', 'parms', 'prior.samp', 'sv.lst'))
clusterEvalQ(cl,{
  library(tidyverse)
  library(tmbstan)
  library(TMB)
  library(lme4)
  library(rstan)
  
  basename <- "regniere_structured"
  basename_loc <- paste0('code/', basename)
  cc <- compile(paste0(basename_loc, '.cpp'))
  try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
  dyn.load(dynlib(basename_loc))
  
  nList <- lme4:::namedList
  stages <- paste0('L', 2:6)
})

gr.lst <- parLapply(cl, 1:20, function(i) {
  chains <- 8
  
  dd <- subset(all.data, prior.samp == i)
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
  
  pnames <- c("phi_rho", "psi_rho", "y0_rho",
              "HA","TL","HL","TH","HH",
              "s_eps", 's_upsilon', "s_alpha",
              'upsilon', 'alpha')
  ff <- MakeADFun(data=dd2,
                  parameters=as.list(parms[pnames]),
                  DLL=basename,
                  random=c('upsilon', 'alpha'),
                  silent=TRUE)
  
  parm.lst <- sv.lst[[i]]
  parm.lst <- lapply(parm.lst, function(x) {
    x$HL <- abs(x$HL)
    x$upsilon <- rep(0, 35)
    x$alpha <- rep(0.05, 5)
    return(x)
  })
  ##### Evaluate Gradient #####
  stan1 <- capture.output({tmbstan(ff, init = parm.lst, silent = TRUE, 
                   chains = 4, iter = 1250, warmup = 750,
                   control= list('max_treedepth' = 15, 'adapt_delta' = 0.9))}, 
                   file = paste0('code/output/iter', i, '.txt'), append = TRUE)
  return(stan1)
  # stan.array <- as.array(stan1)
  # post.df <- as.data.frame(stan1)
  # 
  # post.df$Rhat <- max(apply(stan.array, 3, Rhat))
  # post.df$ESS.bulk <- min(apply(stan.array, 3, ess_bulk))
  # post.df$ESS.tail <- min(apply(stan.array, 3, ess_tail))
  # 
  # g.df <- get_sampler_params(stan1, inc_warmup = FALSE)
  # g.df <- do.call('rbind', g)
  # post.df$divergent <- g.df[,'divergent__']
  # 
  # post.df$iter <- i
  
  #return(post.df)
})
#gr.df <- bind_rows(gr.lst)

stopCluster(cl)

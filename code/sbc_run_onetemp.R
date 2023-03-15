library(tidyverse)
library(tmbstan)
library(TMB)
library(lme4)
library(parallel)

source('code/data_simulation_onetemp.R')

##### Generate and Clean Data #####
set.seed(123)
ptm <- proc.time()
data.lst <- gen.pops(20, sz = 250)
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
names(all.data) <- c('temp1', 'stage', 'time1.orig', 'nobs', 'cup', 'prior.samp')
#all.data$block <- all.data$temp1/5
all.data$block <- paste(all.data$stage, all.data$temp1, sep = '_')
lev.ord <- paste(rep(stages, each = 7), rep(seq(5, 35, by = 5), 5), sep = '_')
all.data$block <- factor(all.data$block)
all.data$block <- factor(all.data$block, levels = lev.ord)
all.data$ind <- paste(all.data$cup, all.data$temp1, all.data$prior.samp, sep = '_')
w2 <- which(all.data$temp1 %in% c(5, 10, 30, 35) & all.data$time1 == 0)
all.data$time1[w2] <- 1

all.data$l2stage <- sapply(all.data$stage, function(x) {ifelse(x == 'L2', 1, 0)})
all.data$time1d <- floor(all.data$time1.orig)
all.data$time1 <- ceiling(all.data$time1.orig) + (1-all.data$l2stage)

write.csv(all.data, 'data/sim_data_onetemp.csv', row.names = FALSE)
write.csv(priors.df, 'data/sim_priors_onetemp.csv', row.names = FALSE)


##### Run Model #####

## Params to initialize MakeADFun
parms <- list(
  'phi_rho' = 0.5,
  'psi_rho' = 0.5,
  'y0_rho' = 0.5,
  'HL' = 8.4,
  'HH' = 9.3,
  'HA' = 1.4, 
  'TL' = 285.90432,
  'TH' = 306.47530,
  's_eps' = rep(0.3, 5),
  's_upsilon' = rep(0.08, 5),
  's_alpha' = 0.22)

samps <- unique(all.data$prior.samp)

basename <- "regniere_structured_onetemp"
basename_loc <- paste0('code/', basename)
cc <- compile(paste0(basename_loc, '.cpp'))
try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
dyn.load(dynlib(basename_loc))

sv.lst <- lapply(samps, function(x) {
  print(x)
  dd <- subset(all.data, prior.samp == x)
  dd$stagename <- dd$stage
  dd$stage <- as.numeric(factor(dd$stage)) - 1
  
  test.lst <- lapply(1:4, function(chain) {
    ind <- (x-1)*4 + chain
    anyna <- TRUE
    while(anyna) {
      ps <- prior.samp(1)[[1]]
      dd2 <- with(dd, nList(temp1,time1,time1d,stage,
                            nobs=as.integer(nobs)))
      dd2$t_block1 <- dd2$stage*7 + dd2$temp1/5 - 1
      dd2$use_prior <- 1
      
      nblock <- length(unique(dd2$stage))*length(unique(dd2$temp1))
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
      
      gr <- grad_log_prob(mod, unlist(ps[1:(length(ps)-1)]))
      #print(gr)
      anyna <- any(is.na(gr))
    }
    return(ps)
  })
  
  return(test.lst)
})

## Set up data for model input
cl1 <- makeCluster(40)
clusterExport(cl1, c('all.data', 'parms', 'prior.samp', 'sv.lst'))
clusterEvalQ(cl1,{
  library(tidyverse)
  library(tmbstan)
  library(TMB)
  library(lme4)
  library(rstan)
  
  basename <- "regniere_structured_onetemp"
  basename_loc <- paste0('code/', basename)
  cc <- compile(paste0(basename_loc, '.cpp'))
  try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
  dyn.load(dynlib(basename_loc))
  
  nList <- lme4:::namedList
  stages <- paste0('L', 2:6)
})

gr.lst1 <- parLapply(cl1, 1:20, function(i) {
  
  dd <- subset(all.data, prior.samp == i)
  dd$stagename <- dd$stage
  dd$stage <- as.numeric(factor(dd$stage)) - 1
  
  dd2 <- with(dd, nList(temp1,time1,time1d,stage,
                        nobs=as.integer(nobs)))
  dd2$t_block1 <- dd2$stage*7 + dd2$temp1/5 - 1
  dd2$use_prior <- 1
  
  nblock <- length(unique(dd2$stage))*length(unique(dd2$temp1))
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
  stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                   chains = 4, iter = 500, warmup = 350)
  post.df <- as.data.frame(stan1)
  write.csv(post.df, paste0('code/output/post_iter_onetemp', i, '.csv'))
  return(stan1)
})
stopCluster(cl1)

diag.lst <- lapply(gr.lst1, function(fit) {
  stan.array <- as.array(fit)
  Rhat <- max(apply(stan.array, 3, Rhat))
  ESS.bulk <- min(apply(stan.array, 3, ess_bulk))
  ESS.tail <- min(apply(stan.array, 3, ess_tail))
  
  g <- get_sampler_params(fit, inc_warmup = FALSE)
  g.df <- do.call('rbind', g)
  div <- sum(g.df[,'divergent__'])
  return(list('Rhat' = Rhat, 
              'ESS.bulk' = ESS.bulk,
              'ESS.tail' = ESS.tail,
              'divergences' = div))
})
diag.df <- bind_rows(diag.lst)
write.csv(diag.df, 'code/output/diagnostics_onetemp.csv')

post.lst <- lapply(gr.lst1, as.data.frame)
post.lst.red <- post.lst


curve.pars <- c('phi_rho', 'psi_rho', 'y0_rho', 
                'HA', 'TL', 'HL', 'TH', 'HH', 's_alpha')
ranks.lst <- lapply(1:20, function(x) {
  prior <- subset(priors.df, prior.samp == x)
  
  prior.cp <- prior[,curve.pars]
  prior.cp <- prior.cp[!duplicated(prior.cp),]
  
  posterior <- post.lst.red[[x]]
  post.cp <- posterior[,curve.pars]
  
  ranks.cp <- sapply(curve.pars, function(cp) {
    prior <- prior.cp[1,cp]
    post <- post.cp[,cp]
    if (cp == 'HL') {
      post <- -abs(post)
      prior <- -abs(prior)
    }
    rnk <- length(which(post < prior))/length(post)
    return(rnk)
  })
  
  post.seps <- posterior[,grep('s_eps', names(posterior))]
  prior.seps <- prior$s_eps
  
  ranks.seps <- sapply(1:length(unique(prior.seps)), function(i) {
    prior <- unique(prior.seps)[i]
    post <- post.seps[,i]
    rnk <- length(which(post < prior))/length(post)
    return(rnk)
  })
  names(ranks.seps) <- names(post.seps)
  
  post.sups <- posterior[,grep('s_upsilon', names(posterior))]
  prior.sups <- prior$s_upsilon
  
  ranks.sups <- sapply(1:length(unique(prior.sups)), function(i) {
    prior <- unique(prior.sups)[i]
    post <- post.sups[,i]
    rnk <- length(which(post < prior))/length(post)
    return(rnk)
  })
  names(ranks.sups) <- names(post.sups)
  
  return(c(ranks.cp, ranks.seps, ranks.sups))
})
ranks.df <- bind_rows(ranks.lst)
coverage <- apply(ranks.df, 2, function(x) {length(which(x > 0.05 & x < 0.95))/length(x)})
write.csv(ranks.df, 'code/output/ranks_onetemp.csv')



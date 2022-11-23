library(tidyverse)
library(tmbstan)
library(TMB)
library(lme4)

source('code/data_simulation.R')

##### Generate and Clean Data #####
set.seed(123)
ptm <- proc.time()
data.lst <- gen.pops(50)
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

##### Run Model #####
basename <- "regniere_structured_test"
basename_loc <- paste0('code/', basename)
cc <- compile(paste0(basename_loc, '.cpp'))
try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
dyn.load(dynlib(basename_loc))
options(mc.cores = parallel::detectCores())

## Params to initialize MakeADFun
parms <- list(
  'phi_rho' = 0.5,
  'psi_rho' = 0.5,
  'y0_rho' = 0.5,
  'HL' = rep(8.4, 5),
  'HH' = rep(9.3, 5),
  'HA' = rep(1.4, 5), 
  'TL' = rep(285.90432, 5),
  'TH' = rep(306.47530, 5),
  's_eps' = rep(0.3, 5),
  's_upsilon' = rep(0.08, 5),
  's_alpha' = 0.22)

chains <- 4
set.seed(124)
parm.lst <- prior.samp(chains)
parm.lst <- lapply(parm.lst, function(x) {
  x$HL <- abs(x$HL)
  x$upsilon <- rep(0, 35)
  x$alpha <- rep(0.05, 5)
  return(x)
})

## Set up data for model input
gr.lst <- list()
for (i in 1:length(data.lst)) {
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
  ##### Evaluate Gradient #####
  stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                   chains = 0, iter = 1400, warmup = 750,
                   control= list('max_treedepth' = 15, 'adapt_delta' = 0.9))
  ps <- priors.lst[[i]]
  ps <- ps[1:(length(ps)-1)]
  
  ps$HL <- abs(ps$HL)
  ps$upsilon <- rep(0, nblock)
  ps$alpha <- rep(0.05, nstage)
  
  gr <- grad_log_prob(stan1, unlist(ps))
  names(gr) <- names(unlist(ps))
  gr.lst <- append(gr.lst, list(gr))
}
gr.df <- bind_rows(gr.lst)

## Sample initial parameter values

# ptm <- proc.time()
# stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE,
#                  chains = chains, iter = 1500, warmup = 500)
# proc.time() - ptm
# 
# g <- get_sampler_params(stan1, inc_warmup = FALSE)
# g.df <- do.call('rbind', g)
# div <- sum(g.df[,'divergent__'])


## Compare priors with working starting points
ps1000 <- prior.samp(1000)
ps1000.df <- bind_rows(lapply(ps1000, as.data.frame))
ps1000.df$prior.samp <- rep(1:1000, each = 5)
ps1000.df$stage <- rep(stages, 1000)

dat.reps.lst <- lapply(1:length(data.lst), function(i) {
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
  stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                   chains = 0, iter = 1500, warmup = 500)
  gr.reps.lst <- lapply(1:1000, function(x) {
    ps <- ps1000[[x]]
    ps <- ps[1:(length(ps)-1)]
    
    ps$HL <- abs(ps$HL)
    ps$upsilon <- rep(0, nblock)
    ps$alpha <- rep(0.05, nstage)
    
    gr <- grad_log_prob(stan1, unlist(ps))
    names(gr) <- names(unlist(ps))
    return(gr)
  })
  gr.reps <- bind_rows(gr.reps.lst)
  gr.reps$dataset <- i
  return(gr.reps)
})
dat.reps <- bind_rows(dat.reps.lst)
dat.reps$prior <- rep(1:1000, length(data.lst))

diag.lst <- lapply(1:5, function(st) {
  nms <- names(dat.reps)[grep(st, names(dat.reps))]
  sub <- dat.reps[,c(nms, 'dataset', 'prior')]
  st.na <- apply(sub, 1, function(x) {any(is.na(x))})
  return(st.na)
})

diag.df <- as.data.frame(diag.lst)
names(diag.df) <- paste0('isna.L', 2:6)
diag.df$dataset <- dat.reps$dataset
diag.df$prior <- dat.reps$prior

##### Clean DF of Posterior Draws #####
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
  
  par.cols <- post.df[,c('HA', 'TL', 'HL', 'TH', 'HH')]
  par.cols.stg <- post.df[, grep(paste0('[', ind, ']'), names(post.df))[1:2]]
  names(par.cols.stg) <- c('s_eps', 's_upsilon')
  pars <- as.data.frame(cbind(q.pars,
                              rho25, 
                              par.cols,
                              par.cols.stg))
  pars$HL <- -pars$HL
  pars$TA <- sapply(1:nrow(pars), function(r) {
    row <- pars[r,]
    ta <- with(row, ta.fun(HL, HH, TL, TH))
    return(ta)
  })
  upsilon.df <- post.df[,grep('^upsilon', names(post.df))]
  end <- ind*7
  start <- end - 6
  udf.stg <- upsilon.df[,start:end]
  #udf.stg <- upsilon.df[,ag$block2[which(ag$stagename == stages[ind])] + 1]
  names(udf.stg) <- paste('upsilon.', seq(5, 35, by = 5), '.', sep = '')
  pars <- cbind(pars, udf.stg)
  
  # stg.df <- subset(dd2.df, stage == ind - 1)
  # stg.df <- stg.df[order(stg.df$temp1),]
  # blocks <- unique(stg.df$block) + 1
  # 
  # temps <- sort(unique(dd2$temp1))
  
  stage.df <- as.data.frame(cbind(pars, 'alpha' = alpha.mult, 
                                  's_alpha' = post.df$s_alpha))
  stage.df$stage <- stages[ind]
  return(stage.df)
})
p.df <- bind_rows(p.lst)

crv.lst <- lapply(stages, function(s) {
  df <- subset(p.df, stage == s)
  gc <- get_curves(df, spline = TRUE)
  gc$stage <- s
  return(gc)
})
crv.df <- bind_rows(crv.lst)

ggplot(data = crv.df) +
  geom_line(aes(x = temp, y = rate, group = index), alpha = 0.4) +
  facet_wrap(vars(stage)) +
  theme_minimal()

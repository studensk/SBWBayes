library(tidyverse)
library(tmbstan)
library(TMB)
library(lme4)

source("code/data_simulation.R")

##### Generate Datasets #####
set.seed(123)
ptm <- proc.time()
data.lst <- gen.pops(10)
proc.time() - ptm

priors.lst <- lapply(1:length(data.lst), function(i) {
  x <- data.lst[[i]]
  dat <- x[[1]]
  dat$prior.samp <- i
  return(dat)
})
all.data <- bind_rows(lapply(data.lst, function(x) {x[[2]]}))
names(all.data) <- c('temp1', 'temp2', 'stage', 'time1', 'time2', 'nobs', 'prior.samp')
all.data$block <- all.data$temp1/5
all.data$time2d <- pmax(0.25, all.data$time2 - 1)

## Correction for early moulters
early <- which(all.data$time2 == 0)
early.sub <- all.data[early,]
early.sub$time2 <- early.sub$time1
early.sub$time1 <- 0
early.sub$time2d <- pmax(early.sub$time2 - 1, 0.25)
early.sub$temp2 <- early.sub$temp1

all.data[early,] <- early.sub
all.data$time2 <- pmax(all.data$time2, 1)

##### Model Runs #####
basename <- "regniere_structured"
basename_loc <- paste0('code/', basename)
cc <- compile(paste0(basename_loc, '.cpp'))
try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
dyn.load(dynlib(basename_loc))
options(mc.cores = parallel::detectCores())

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
  's_upsilon' = rep(0.08, 5),
  's_alpha' = 0.22)

diagnostics <- list()
fit.lst <- list()

chains <- 4
set.seed(124)
parm.lst <- prior.samp(chains)

dd <- subset(all.data, prior.samp == 1)
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
  "s_eps", 's_upsilon', "s_alpha",
  'upsilon', 'alpha')
ff <- MakeADFun(data=dd2,
                parameters=as.list(parms[pnames]),
                DLL=basename,
                random=c('upsilon', 'alpha'),
                silent=TRUE)

parm.lst <- lapply(parm.lst, function(x) {
  x$upsilon <- rep(0, nblock)
  x$alpha <- rep(0.05, nstage)
  x$HL <- -x$HL
  return(x)
})
stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                 chains = chains, iter = 2000, warmup = 1000,
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

##### Prelim SBC #####
rankfun <- function(x, vec) {length(which(vec < x))/length(vec)}

l <- lapply(unique(priors.df$prior.samp), function(ind) {
  sub <- subset(priors.df, prior.samp == ind)
  sub$HL <- -sub$HL
  post <- fit.lst[[ind]]
  outer.vars <- c('phi_rho', 'psi_rho', 'y0_rho', 's_alpha')
  ranks <- sapply(outer.vars, function(y) {
    vec <- post[,y]
    x <- as.numeric(unique(sub[,y]))
    rf <- rankfun(x, vec)
    return(rf)
  })
  lp <- lapply(1:5, function(stg) {
    row <- sub[stg,]
    stgind <- paste0('[', stg, ']')
    vars <- names(post)[grep(stgind, names(post))]
    vars <- vars[1:(length(vars)-2)]
    ranks <- sapply(vars, function(y) {
      vec <- post[,y]
      nm <- strsplit(y, '\\[')[[1]][1]
      x <- as.numeric(sub[stg, nm])
      rf <- rankfun(x, vec)
      return(rf)
    })
  })
  ranks.all <- c(ranks, unlist(lp))
  return(ranks.all)
})
sbc.df <- as.data.frame(bind_rows(l))

##### Test gradients ######
stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                 chains = 0, iter = 1500, warmup = 500)
ps1000 <- prior.samp(1000)

test.lst <- lapply(1:1000, function(i) {
  ps.orig <- ps1000[[i]]
  ps <- ps.orig
  ps$upsilon <- rep(0, nblock)
  ps$alpha <- rep(0.05, nstage)
  ps$HL <- -ps$HL
  gr <- grad_log_prob(stan1, unlist(ps))
  dat <- data.frame('parameter' = names(unlist(ps)), 'value' = unlist(ps), gr)
  dat$run <- i
  return(dat)
})
test.lst2 <- lapply(1:1000, function(i) {
  ps.orig <- ps1000[[i]]
  ps <- ps.orig
  ps$upsilon <- rep(0, nblock)
  ps$alpha <- rep(0.05, nstage)
  ps$HL <- -ps$HL
  gr <- grad_log_prob(stan1, unlist(ps))
  w <- which(is.na(gr))
  stgs <- names(unlist(ps))[w][grep('HL', names(unlist(ps))[w])]
  hl.ind <- paste0('HL', 1:5)
  ps.orig$status <- rep('pass', 5)
  ps.orig$status[which(hl.ind %in% stgs)] <- 'fail'
  ps.df <- as.data.frame(ps.orig)
  ps.df$stage <- stages
  ps.df$run <- i
  return(ps.df)
})
test.df <- bind_rows(test.lst)
test.df2 <- bind_rows(test.lst2)
test.df$status <- sapply(test.df$gr, function(x) {ifelse(is.na(x), 'fail', 'pass')})

test.df2$rho25 <- sapply(1:nrow(test.df2), function(s) {
  stg <- test.df2$stage[s]
  w <- which(stages == stg)
  stage.n <- w - 3
  rho <- with(test.df2[s,], quadratic(stage.n, phi_rho, psi_rho, y0_rho))
  return(rho)
})
test.df2$TA <- mapply(function(hl, hh, tl, th) {
  ta.fun(hl, hh, tl, th)
}, test.df2$HL, test.df2$HH, test.df2$TL, test.df2$TH)

test.df.ord <- test.df2[order(test.df2$status),]
vals <- seq(0, 40, by = 1)
crv <- with(test.df.ord[1,], {calc_pred3(vals, rho25, HA, TL, HL, TH, HH, TA)})
plot(vals, log(crv), col = ifelse(test.df.ord$status[1] == 'fail', 'red', 'black'), 
     type = 'l', ylim = c(0, log(1.069695e+19)))
m.crv <- 0
for (i in 2:5000) {
  crv <- with(test.df.ord[i,], {calc_pred3(vals, rho25, HA, TL, HL, TH, HH, TA)})
  m.crv <- pmax(m.crv, max(crv))
  lines(vals, log(crv), col = ifelse(test.df.ord$status[i] == 'fail', 'red', 'black'))
}

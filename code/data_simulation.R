source("code/functions.R")
source('code/test_stan.R')

all.days.df <- read.csv('data/all_days_df.csv')
nList <- lme4:::namedList
stages <- paste0('L', 2:6)

##### Prior Sampling #####
# prior.samp <- function(chains) {
#   l <- lapply(1:chains, function(x) {
#     lst <- list(
#       'phi_rho' = rbeta(1, 4, 4),
#       'psi_rho' = rbeta(1, 4, 4),
#       'y0_rho' = rgamma(1, 5, scale = 0.1),
#       'HA' = rgamma(5, 5, scale = 0.2),
#       'TL' = rnorm(5, 284, 2),
#       'HL' = -rgamma(5, 6, scale = 2),
#       'TH' = rnorm(5, 304, 2),
#       'HH' = rgamma(5, 10, scale = 3),
#       's_eps' = exp(rnorm(5, -1.5, 0.1)),
#       's_upsilon' = exp(rnorm(5, -2.5, 0.05)),
#       's_alpha' = exp(rnorm(1, -1.5, 0.3)),
#       'alpha' = 0.05
#     )
#     return(lst)
#   })
#   return(l)
# }
prior.samp <- function(chains) {
  l <- lapply(1:chains, function(x) {
    tl <- exp(rnorm(1, 5.64, 0.0067))
    tdiff <- rgamma(1, 112, scale = 0.226) 
    th <- tl + tdiff
    lst <- list(
      'phi_rho' = rbeta(1, 4, 4),
      'psi_rho' = rbeta(1, 4, 4),
      'y0_rho' = rgamma(1, 8.3, scale = 0.046),
      'HA' = rgamma(1, 5.4, scale = 0.134),
      'TL' = tl,
      'HL' = -rgamma(1, 3.6, scale = 2.253),
      'TH' = th,
      'HH' = rgamma(1, 7.6, scale = 3.12),
      's_eps' = exp(rnorm(5, -1.5, 0.1)),
      's_upsilon' = exp(rnorm(5, -2.5, 0.05)),
      's_alpha' = exp(rnorm(1, -1.5, 0.3)),
      'alpha' = 0.05
    )
    return(lst)
  })
  return(l)
}
ps <- prior.samp(1000)

##### Data Generation #####
on.add <- subset(all.days.df, province == 'ON' & transfer)
on.add$time1 <- round(on.add$time1)
treat.times <- aggregate(data = on.add, time1 ~ temp1 + stage, 
                         function(x) {which.max(tabulate(x))})

ind.progress <- function(param.df, dev.df) {
  sustainable <- ifelse(dev.df$treatment[1] == dev.df$sustainable[1], TRUE, FALSE)
  days <- vector()
  eps <- rnorm(5, 0, param.df$s_eps)
  deltas <- exp(eps)
  names(deltas) <- stages
  dev <- 0
  i <- 1
  dev.lst <- list()
  extra.dev <- 0
  dev.df$r.treatment <- dev.df$r.treatment*deltas
  dev.df$r.sustainable <- dev.df$r.sustainable*deltas
  
  for (s in 1:length(stages)) {
    st <- stages[s]
    rdf <- subset(dev.df, stage == st)
    trt.time <- (1/rdf$r.treatment)*(1-extra.dev)
    trt.time <- pmax(trt.time, 0)
    
    if (sustainable) {
      obs.time <- ceiling(trt.time)
      rdf$time.treatment <- 0
      rdf$time.sustainable <- obs.time
      if (st != 'L6') {
        extra.time <- obs.time - trt.time
        rdf.new <- subset(dev.df, stage == stages[s+1])
        extra.dev <- extra.time*rdf.new$r.sustainable
      }
      dev.lst <- append(dev.lst, list(rdf))
      next
    }
    
    if (trt.time > rdf$time.treatment) {
      rem.prop <- 1 - extra.dev - rdf$r.treatment*rdf$time.treatment
      sus.time <- rem.prop/rdf$r.sustainable
      obs.time <- ceiling(sus.time)
      rdf$time.sustainable <- obs.time
      
      if (st != 'L6') {
        extra.time <- obs.time - sus.time
        rdf.new <- subset(dev.df, stage == stages[s + 1])
        extra.dev <- extra.time*rdf.new$r.sustainable
      }
    }
    
    else {
      rdf$time.treatment <- ceiling(trt.time)
      rdf$time.sustainable <- 0
      if (st != 'L6') {
        if (extra.dev > 1) {extra.time <- extra.time - 1/rdf$r.sustainable}
        else {extra.time <- ceiling(trt.time) - trt.time}
        rdf.new <- subset(dev.df, stage == stages[s + 1])
        extra.dev <- extra.time*rdf.new$r.sustainable
      }
    }
    dev.lst <- append(dev.lst, list(rdf))
  }
  data <- bind_rows(dev.lst)
  data$nobs <- 1
  return(data)
  
}

pop.dev <- function(t.treat, param.df, size = 250) {
  if (t.treat %in% c(5, 10, 30, 35)) {all.temps <- c('treatment' = t.treat, 'sustainable' = 20)}
  else {all.temps <- c('treatment' = t.treat, 'sustainable' = t.treat)}
  dev.lst <- lapply(stages, function(st) {
    s <- which(stages == st)
    params <- param.df[s,]
    times <- with(params, {
      rho25 <- quadratic(s-3, phi_rho, psi_rho, y0_rho)
      TA <- ta.fun(HL, HH, TL, TH)
      calc_pred3(all.temps, rho25, HA, TL, HL, TH, HH, TA)})
    rates <- 1/times
    names(rates) <- paste('r', names(all.temps))
    s.df <- as.data.frame(append(as.list(all.temps), as.list(rates)))
    s.df$stage <- st
    if (t.treat %in% c(5, 10, 30, 35))
    {s.df$time.treatment <- subset(treat.times, temp1 == all.temps['treatment'] & 
                                     stage == st)$time1}
    return(s.df)
  })
  dev.df <- bind_rows(dev.lst)
  l <- replicate(size, {
    ind.progress(param.df, dev.df)
  }, simplify = FALSE)
  df <- bind_rows(l)
  nms <- c('treatment', 'sustainable', 'stage', 'time.treatment', 'time.sustainable', 'nobs')
  df <- df[,nms]
  ag <- aggregate(data = df, nobs ~ ., sum)
  return(ag)
}

gen.pops <- function(N) {
  ps <- prior.samp(N)
  lst <- lapply(1:N, function(i) {
    p <- ps[[i]]
    t.lst <- lapply(seq(5, 35, by = 5), function(tmp) {
      pdat <- pop.dev(tmp, as.data.frame(p))
    })
    pd <- bind_rows(t.lst)
    pd$prior.samp <- i
    return(list('priors' = p, 'data' = pd))
  })
  return(lst)
}

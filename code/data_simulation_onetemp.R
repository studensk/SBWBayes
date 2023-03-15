source("code/functions.R")
source('code/test_stan.R')

all.days.df <- read.csv('data/all_days_df.csv')
nList <- lme4:::namedList
stages <- paste0('L', 2:6)

##### Prior Sampling #####

prior.samp <- function(chains) {
  l <- lapply(1:chains, function(x) {
    tl <- exp(rnorm(1, 5.64, 0.0067))
    tdiff <- rgamma(1, 112, scale = 0.226) 
    th <- tl + tdiff
    s_upsilon <- exp(rnorm(5, -2.5, 0.05))
    s_alpha <- exp(rnorm(1, -1.5, 0.3))
    sa2 <- (-s_alpha^2)/2
    alpha <- rnorm(5, sa2, s_alpha)
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
      's_upsilon' = s_upsilon,
      's_alpha' = s_alpha,
      #'upsilon' = rlnorm(35, 0, rep(s_upsilon, each = 7)),
      'upsilon' = exp(rcauchy(35, 0, rep(s_upsilon, 7))),
      'alpha' = alpha,
      'zeta' = exp(alpha)
    )
    return(lst)
  })
  return(l)
}

##### Data Generation #####

ind.progress <- function(param.df, dev.df) {
  days <- vector()
  eps <- rnorm(5, 0, param.df$s_eps)
  deltas <- exp(eps)
  names(deltas) <- stages
  dev <- 0
  i <- 1
  dev.lst <- list()
  extra.dev <- 0
  dev.df$r.treatment <- dev.df$r.treatment/deltas
  
  for (s in 1:length(stages)) {
    st <- stages[s]
    rdf <- subset(dev.df, stage == st)
    trt.time <- 1/rdf$r.treatment
    #rdf$time.treatment <- ceiling(trt.time)
    rdf$time.treatment <- trt.time
    dev.lst <- append(dev.lst, list(rdf))
  }
  data <- bind_rows(dev.lst)
  data$nobs <- 1
  return(data)
  
}

pop.dev <- function(t.treat, param.df, size = 250) {
  all.temps <- c('treatment' = t.treat)
  dev.lst <- lapply(stages, function(st) {
    s <- which(stages == st)
    params <- subset(param.df, stage == st & temp == t.treat)
    times <- with(params, {
      rho25 <- quadratic(s-3, phi_rho, psi_rho, y0_rho)*zeta
      TA <- ta.fun(HL, HH, TL, TH)
      cp1 <- calc_pred3(all.temps[1], rho25, HA, TL, HL, TH, HH, TA)*upsilon
      cvec <- unique(cp1)
      names(cvec) <- names(all.temps)
      return(cvec)})
    rates <- 1/times
    names(rates) <- paste('r', names(all.temps))
    s.df <- as.data.frame(append(as.list(all.temps), as.list(rates)))
    s.df$stage <- st
    return(s.df)
  })
  dev.df <- bind_rows(dev.lst)
  l <- replicate(size, {
    ind.progress(param.df, dev.df)
  }, simplify = FALSE)
  df <- bind_rows(l)
  nms <- c('treatment', 'stage', 'time.treatment', 'nobs')
  df <- df[,nms]
  df$cup <- rep(1:size, each = 5)
  #ag <- aggregate(data = df, nobs ~ ., sum)
  #return(ag)
  return(df)
}

gen.pops <- function(N, sz = 250) {
  ps <- prior.samp(N)
  lst <- lapply(1:N, function(i) {
    p <- ps[[i]]
    pdf <- as.data.frame(p)
    pdf$stage <- rep(stages, 7)
    pdf$temp <- rep(seq(5, 35, by = 5), each = 5)
    t.lst <- lapply(seq(5, 35, by = 5), function(tmp) {
      #psub <- subset(pdf, temp == tmp)
      pdat <- pop.dev(tmp, pdf, size = sz)
    })
    pd <- bind_rows(t.lst)
    pd$prior.samp <- i
    return(list('priors' = p, 'data' = pd))
  })
  return(lst)
}

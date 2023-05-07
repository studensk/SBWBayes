library(tidyverse)
library(tmbstan)
library(TMB)
library(lme4)
library(parallel)
library(reshape2)
library(cat)

source('code/noss_data_simulation.R')

##### Generate and Clean Data #####
set.seed(123)
ptm <- proc.time()

cl <- makeCluster(50)
clusterEvalQ(cl, {
  library(tidyverse)
  source('code/noss_data_simulation.R')
})
data.lst <- parLapply(cl,1:1000, function(x) {
  set.seed(x + 1)
  gp <- gen.pops(1)[[1]]
  gp$data$prior.samp <- x
  return(gp)
})
stopCluster(cl)

priors.lst <- lapply(1:length(data.lst), function(i) {
  x <- data.lst[[i]]
  dat <- x[[1]]
  dat$prior.samp <- i
  return(dat)
})
priors.df <- bind_rows(lapply(priors.lst, as.data.frame))
priors.df$stage <- rep(stages, length(priors.lst))

all.data <- bind_rows(lapply(data.lst, function(x) {x[[2]]}))
names(all.data) <- c('temp1', 'temp2', 'stage', 'time1', 'time2', 'nobs', 'cup', 'prior.samp')
#all.data$block <- all.data$temp1/5
all.data$block <- paste(all.data$stage, all.data$temp1, sep = '_')
lev.ord <- paste(rep(stages, each = 7), rep(seq(5, 35, by = 5), 5), sep = '_')
all.data$block <- factor(all.data$block)
all.data$block <- factor(all.data$block, levels = lev.ord)
all.data$ind <- paste(all.data$cup, all.data$temp1, all.data$prior.samp, sep = '_')
w1 <- which(all.data$temp1 %in% c(15, 20, 25) & all.data$time2 == 0)
all.data$time2[w1] <- 1
w2 <- which(all.data$temp1 %in% c(5, 10, 30, 35) & all.data$time1 == 0)
all.data$time1[w2] <- 1

all.data$time1.orig <- all.data$time1
all.data$time2.orig <- all.data$time2

all.data$l2 <- sapply(all.data$stage, function(x) {ifelse(x == 'L2', 1, 0)})
all.data$cur0 <- sapply(all.data$time2.orig, function(x) {ifelse(x == 0, 1, 0)})
all.data$prev0 <- c(0, all.data$cur0[1:(nrow(all.data) - 1)])

all.data$time1 <- all.data$time1.orig + (1-all.data$l2)*all.data$prev0
all.data$time2 <- all.data$time2.orig + (1-all.data$l2)*(1-all.data$prev0)
all.data$time1d <- all.data$time1.orig - all.data$cur0
all.data$time2d <- all.data$time2.orig - (1 - all.data$cur0)
all.data.orig <- all.data

all.data <- aggregate(data = all.data, nobs ~ temp1 + temp2 + stage + 
                        time1 + time2 + time1d + time2d + prior.samp + block,
                      sum)

all.data <- subset(all.data, time2 <= 60)
proc.time() - ptm

write.csv(all.data, 'data/noss_sim_data_reduced.csv')
write.csv(priors.df, 'data/noss_sim_priors_reduced.csv')

##### Run Model #####

## Params to initialize MakeADFun
parms <- list(
  'rho' = 0.045*(-1.05*(0:4)^2 + 4.22*(0:4) + 4.08),
  'HL' = 8.4,
  'HH' = 9.3,
  'HA' = 1.4, 
  'TL' = 285.90432,
  'TH' = 306.47530,
  's_eps' = rep(0.3, 5),
  's_upsilon' = rep(0.08, 5))

samps <- unique(all.data$prior.samp)

basename <- "noss_regniere_structured_reduced"
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
    index <- 0
    while(anyna) {
      print(ind)
      index <- index + 1
      set.seed(ind + (index-1)*4000)
      ps <- prior.samp(1)[[1]]
      dd2 <- with(dd, nList(temp1,time1,temp2,time2,time1d,time2d,stage,
                            nobs=as.integer(nobs)))
      dd2$t_block1 <- dd2$stage*7 + dd2$temp1/5 - 1
      dd2$t_block2 <- dd2$stage*7 + dd2$temp2/5 - 1
      sm1 <- 0:4
      dd2$k_vec <- -1.05*sm1^2 + 4.22*sm1 + 4.08
      dd2$use_prior <- 1
      
      nblock <- length(unique(dd2$stage))*length(unique(dd2$temp1))
      nstage <- length(unique(dd2$stage))
      
      parms$upsilon <- rep(0, nblock)
      
      pnames <- c("rho",
        "HA","TL","HL","TH","HH",
        "s_eps", 's_upsilon',
        'upsilon')
      ff <- MakeADFun(data=dd2,
                      parameters=as.list(parms[pnames]),
                      DLL=basename,
                      random=c('upsilon'),
                      silent=TRUE)
      mod <- tmbstan(ff, init = unlist(parms[1:(length(parms)-1)]), silent = TRUE, 
                     chains = 0)
      
      nm.ord <- c(unique(names(ff$par)), 'upsilon')
      
      ps$upsilon <- rep(0, nblock)
      ps$HL <- abs(ps$HL)
      
      eval <- unlist(ps[nm.ord])
      nms <- c(names(ff$par), rep('upsilon', nblock))
      names(eval) <- nms
      
      gr <- grad_log_prob(mod, eval)
      nms <- names(unlist(ps[1:(length(ps))]))
      print(nms[which(is.na(gr))])
      #print(gr)
      anyna <- any(is.na(gr))
      if (index > 20) {return(NULL)}
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
  
  basename <- "noss_regniere_structured_reduced"
  basename_loc <- paste0('code/', basename)
  cc <- compile(paste0(basename_loc, '.cpp'))
  try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
  dyn.load(dynlib(basename_loc))
  
  nList <- lme4:::namedList
  stages <- paste0('L', 2:6)
})

gr.lst1 <- parLapply(cl1, 751:1000, function(i) {
  
  dd <- subset(all.data, prior.samp == i)
  dd$stagename <- dd$stage
  dd$stage <- as.numeric(factor(dd$stage)) - 1
  
  dd2 <- with(dd, nList(temp1,time1,temp2,time2,time1d,time2d,stage,
                        nobs=as.integer(nobs)))
  dd2$t_block1 <- dd2$stage*7 + dd2$temp1/5 - 1
  dd2$t_block2 <- dd2$stage*7 + dd2$temp2/5 - 1
  sm1 <- 0:4
  dd2$k_vec <- -1.05*sm1^2 + 4.22*sm1 + 4.08
  dd2$use_prior <- 1
  
  nblock <- length(unique(dd2$stage))*length(unique(dd2$temp1))
  nstage <- length(unique(dd2$stage))
  
  parms$upsilon <- rep(0, nblock)
  parms$alpha <- rep(0, nstage)
  
  pnames <- c("rho",
              "HA","TL","HL","TH","HH",
              "s_eps", 's_upsilon',
              'upsilon')
  ff <- MakeADFun(data=dd2,
                  parameters=as.list(parms[pnames]),
                  DLL=basename,
                  random=c('upsilon'),
                  silent=TRUE)
  
  parm.lst <- sv.lst[[i]]
  parm.lst <- lapply(parm.lst, function(x) {
    x$HL <- abs(x$HL)
    x$upsilon <- rep(0, 35)
    return(x)
  })
  stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                   chains = 4, iter = 1250, warmup = 750,
                   control = list(max_treedepth = 15, adapt_delta = 0.99))
  post.df <- as.data.frame(stan1)
  write.csv(post.df, paste0('code/final/output/noss/post_iter', i, '.csv'))
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
write.csv(diag.df, 'code/final/diagnostics_reduced3.csv')

# post.lst.red <- lapply(post.lst, function(x) {
#   if (nrow(x) > 0) {
#     rows <- sample(1:nrow(x), 200)
#     return(x[rows,])
#   }
#   else {return(NULL)}
# })

w <- which(diag.df$divergences == 0)

curve.pars <- c('HA', 'TL', 'HL', 'TH', 'HH')
priors.df.new <- priors.df
priors.df.new$prior.samp <- as.numeric(factor(priors.df.new$prior.samp))

# rank.vec <- function(prior, posterior, varname) {
#   upr <- unique(prior[,varname])
#   post.var <- posterior[,grep(varname, names(posterior))]
#   ranks <- sapply(1:length(upr), function(i) {
#     pr <- upr[i]
#     post <- post.var[,i]
#     rnk <- length(which(post < pr))/length(post)
#     return(rnk)
#   })
#   names(ranks) <- names(post.var)
#   return(ranks)
# }
# 
# ranks.lst <- lapply(1:length(post.lst.red), function(x) {
#   prior <- subset(priors.df.new, prior.samp == w[x])
#   
#   prior.cp <- prior[,curve.pars]
#   prior.cp <- prior.cp[!duplicated(prior.cp),]
#   
#   posterior <- post.lst.red[[x]]
#   post.cp <- posterior[,curve.pars]
#   
#   ranks.cp <- sapply(curve.pars, function(cp) {
#     prior <- prior.cp[1,cp]
#     post <- post.cp[,cp]
#     if (cp == 'HL') {
#       post <- -abs(post)
#       prior <- -abs(prior)
#     }
#     rnk <- length(which(post < prior))/length(post)
#     return(rnk)
#   })
#   
#   ranks.stg.lst <- lapply(c('rho', 's_eps', 's_upsilon'), function(x) {
#     rank.vec(prior, posterior, x)
#   })
#   
#   return(c(ranks.cp, unlist(ranks.stg.lst)))
# })
# ranks.df <- bind_rows(ranks.lst)
# coverage <- apply(ranks.df, 2, function(x) {length(which(x > 0.05 & x < 0.95))/length(x)})
# write.csv(ranks.df, 'code/final/ranks_reduced1.csv')

post.lst.thin <- lapply(post.lst.red, function(x) {
  if (nrow(x) > 0) {
    rows <- seq(1, nrow(x), length.out = 35)
    return(x[rows,])
  }
  else {return(NULL)}
})

rank.vec <- function(prior, posterior, varname) {
  upr <- unique(prior[,varname])
  post.var <- posterior[,grep(varname, names(posterior))]
  ranks <- sapply(1:length(upr), function(i) {
    pr <- upr[i]
    post <- post.var[,i]
    rnk <- length(which(post < pr))
    return(rnk)
  })
  names(ranks) <- names(post.var)
  return(ranks)
}

ranks.lst <- lapply(1:length(post.lst.thin), function(x) {
  prior <- subset(priors.df.new, prior.samp == w[x])
  
  prior.cp <- prior[,curve.pars]
  prior.cp <- prior.cp[!duplicated(prior.cp),]
  
  posterior <- post.lst.thin[[x]]
  post.cp <- posterior[,curve.pars]
  
  ranks.cp <- sapply(curve.pars, function(cp) {
    prior <- prior.cp[1,cp]
    post <- post.cp[,cp]
    if (cp == 'HL') {
      post <- -abs(post)
      prior <- -abs(prior)
    }
    rnk <- length(which(post < prior))
    return(rnk)
  })
  
  ranks.stg.lst <- lapply(c('rho', 's_eps', 's_upsilon'), function(x) {
    rank.vec(prior, posterior, x)
  })
  
  return(c(ranks.cp, unlist(ranks.stg.lst)))
})
ranks.df <- bind_rows(ranks.lst)
coverage <- apply(ranks.df, 2, function(x) {length(which(x > 0.05 & x < 0.95))/length(x)})
write.csv(ranks.df, 'code/final/ranks_reduced1.csv')

par(mfrow = c(3, 3))
for (i in 1:ncol(ranks.df)) {
  hist(as.data.frame(ranks.df)[,i], main = names(ranks.df)[i], breaks = 35)
}

##### SBC on devrates #####
pr.lst <- lapply(w, function(ps) {
  sub <- subset(priors.df, prior.samp == ps)
  st.lst <- lapply(stages, function(x) {
    sub2 <- subset(sub, stage == x, select = -X)
    rw <- sub2[1,c('rho', 'HA', 'TL', 'HL', 'TH', 'HH', 's_eps', 's_upsilon')]
    
    ups.orig <- sub2$upsilon
    names(ups.orig) <- paste0('upsilon.', seq(5, 35, by = 5), '.')
    p.ups <- pnorm(log(ups.orig), 0, rw$s_upsilon)
    ups <- qnorm(p.ups)
    df <- as.data.frame(c(rw, ups, 'rho25' = rw$rho))
    df$stage <- x
    return(df)
  })
  st.df <- bind_rows(st.lst)
  st.df$prior.samp <- ps
  return(st.df)
})
pr.df <- bind_rows(pr.lst)

pr.gc <- get_curves(pr.df, temp = seq(5, 35, by = 5), spline = TRUE)
pr.gc$stage <- rep(pr.df$stage, each = 7)
pr.gc$index <- rep(pr.df$prior.samp, each = 7)
d.pr.gc <- dcast(data = pr.gc, index ~ temp + stage, value.var = 'rate')
names(d.pr.gc) <- c('index', paste0('rate', names(d.pr.gc)[2:length(d.pr.gc)]))
write.csv(d.pr.gc, paste0('code/final/output/noss/devrates/prior_rates.csv'))

rateCluster <- makeCluster(10)
clusterEvalQ(rateCluster, {
  library(tidyverse)
  library(reshape2)
  source('code/functions.R')
})
clusterExport(rateCluster, {
  c('stages', 'w')
})
rate.lst <- parLapply(rateCluster, w, function(x) {
  post.df <- read.csv(paste0('code/final/output/noss/post_iter', x, '.csv'))
  post.df <- post.df[,2:(length(post.df)-1)]
  
  p.lst <- lapply(1:5, function(ind) {
    cnms <- c('HA', 'TL', 'HL', 'TH', 'HH')
    par.cols <- post.df[,cnms]
    stg.cols <- post.df[, grep(paste0('.', ind, '.'), names(post.df))[1:3]]
    pars <- as.data.frame(cbind(par.cols, stg.cols))
    names(pars) <- c(cnms, 'rho25', 's_eps', 's_upsilon')
    pars$HL <- -abs(pars$HL)
    pars$TA <- sapply(1:nrow(pars), function(r) {
      row <- pars[r,]
      ta <- with(row, ta.fun(HL, HH, TL, TH))
      return(ta)
    })
    upsilons <- post.df[,grep('^upsilon', names(post.df))]
    inds <- (ind-1)*7 + 1:7
    ups.ind <- upsilons[,inds]
    names(ups.ind) <- paste0("upsilon.", seq(5, 35, by = 5), '.')
    
    stage.df <- as.data.frame(cbind(pars, ups.ind))
    stage.df$stage <- stages[ind]
    stage.df$index <- 1:nrow(stage.df)
    return(stage.df)
  })
  p.df <- bind_rows(p.lst)
  
  gc <- get_curves(p.df, temp = seq(5, 35, by = 5), spline = TRUE)
  gc$stage <- rep(p.df$stage, each = 7)
  gc$index <- rep(p.df$index, each = 7)
  d.gc <- dcast(data = gc, index ~ temp + stage, value.var = 'rate')
  names(d.gc) <- c('index', paste0('rate', names(d.gc)[2:length(d.gc)]))
  write.csv(d.gc, paste0('code/final/output/noss/devrates/ratedf', x, '.csv'))
  d.gc$over_ind <- x
  return(d.gc)
})
rate.df <- bind_rows(rate.lst)
stopCluster(rateCluster)

rate.lst.thin <- lapply(rate.lst, function(x) {
  if (nrow(x) > 0) {
    rows <- seq(1, nrow(x), length.out = L)
    return(x[rows,])
  }
  else {return(NULL)}
})
rate.df.thin <- subset(bind_rows(rate.lst.thin), select = -c(index, over_ind))
pr.rate.df <- subset(d.pr.gc, select = -index)
rate.df.thin$sample <- 'posterior'
pr.rate.df$sample <- 'prior'

all.rates <- bind_rows(pr.rate.df, rate.df.thin)
m.all.rates <- melt(all.rates, id.vars = 'sample')

ggplot(data = m.all.rates) +
  geom_boxplot(aes(x = sample, y = value)) +
  facet_wrap(vars(variable), scales = 'free', nrow = 7, ncol = 5) +
  theme_minimal()

r.lst <- lapply(1:length(w), function(r) {
  x <- w[r]
  prior <- unlist(d.pr.gc[r,grep('rate', names(d.pr.gc))])
  post <- rate.lst.thin[[r]]
  post <- post[,grep('rate', names(post))]
  ranks <- lapply(1:length(prior), function(i) {
    prr <- prior[i]
    rank <- length(which(post[,i] < prr))
    return(rank)
  })
  df <- as.data.frame(ranks)
  names(df) <- names(prior)
  return(df)
})
r.df <- bind_rows(r.lst)
coverage.rates <- apply(r.df, 2, function(x) {length(which(x > 0.05 & x < 0.95))/length(x)})

m.rdf <- melt(r.df)
m.rdf$temp <- sapply(as.character(m.rdf$variable), function(x) {
  str <- substr(x, 5, nchar(x))
  strsplit(str, split = '_')[[1]][1]
})
m.rdf$stage <- sapply(as.character(m.rdf$variable), function(x) {strsplit(x, split = '_')[[1]][2]})

ggplot() +
  geom_hline(yintercept = seq(9, 32, by = 0.1), alpha = 0.5, col = 'lightgrey') +
  geom_histogram(data = m.rdf, aes(x = value), bins = 27, col = 'black') +
  facet_wrap(vars(variable), scales = 'free', nrow = 7, ncol = 5) +
  theme_minimal() 

# for (s in stages) {
#   par(mfrow = c(3, 3))
#   dat <- r.df[,grep(s, names(r.df))]
#   for (i in 1:ncol(dat)) {
#     hist(as.data.frame(dat)[,i], main = names(dat)[i], breaks = 30)
#     abline(h = 9, col = 'red')
#     abline(h = 19)
#     abline(h = 32, col = 'red')
#   }
# }

## Power transformation params
post.lst.red.ind <- lapply(1:length(post.lst.red), function(i) {
  x <- post.lst.red[[i]]
  x$ind <- w[i]
  return(x)
})
## exp(sc.rat*(BC(y, lambda) - mu.post) + mu.prior)
post.lst.all <- bind_rows(post.lst.red.ind)
seps.par.lst <- lapply(1:5, function(x) {
  post <- post.lst.all[,paste0('s_eps.', x, '.')]
  lambda <- powerTransform(post)$lambda
  pt.post <- powertrans(post)
  
  mu.post <- median(pt.post)
  sc.post <- scale.est(pt.post)
  
  mu.prior <- -1.5
  sc.prior <- 0.1
  
  sc.rat <- sc.prior/sc.post
  df <- data.frame(lambda, mu.post, mu.prior, sc.rat)
  df$stage <- x
  return(df)
})
seps.par <- bind_rows(seps.par.lst)
write.csv(seps.par, 'code/final/s_eps_trans_pars.csv', row.names = FALSE)

seps.lst <- lapply(1:5, function(x) {
  post <- post.lst.all[,paste0('s_eps.', x, '.')]
  row <- seps.par[x,]
  trans <- with(row, {exp(sc.rat*(powertrans(post, lambda) - mu.post) + mu.prior)})
  return(trans)
})
seps.df <- bind_cols(seps.lst)
names(seps.df) <- paste0('s_eps.', 1:5, '.')

seps.df$prior.samp <- rep(w, each = 2000)
seps.ranks.lst <- lapply(1:5, function(x) {
  prior.df <- subset(priors.df, stage == stages[x])
  prior <- unique(prior.df$s_eps)
  ranks <- sapply(1:length(prior), function(i) {
    prr <- prior[i]
    post <- subset(seps.df, prior.samp == w[i])
    rank <- length(which(post[,x] < prr))/nrow(post)
    return(rank)
  })
  return(ranks)
})
names(seps.ranks.lst) <- paste0('s_eps.', 1:5, '.')

qb <- qbinom(c(0.005, 0.5, 0.995), length(w), 10/round(length(w)))
qbn <- data.frame('x' = c(-1, 36), 'ymin' = rep(qb[1], 2),
                  'ymax' = rep(qb[3], 2), 'med' = rep(qb[2], 2))
ggplot(data = qbn) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), fill = 'lightgrey') +
  geom_line(aes(x = x, y = med), col = 'darkgrey') +
  geom_histogram(data = m.rdf, aes(value), bins = round(length(w)/20)) +
  facet_wrap(vars(variable), scales = 'free') +
  theme_minimal()

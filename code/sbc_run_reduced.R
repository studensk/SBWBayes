library(tidyverse)
library(tmbstan)
library(TMB)
library(lme4)
library(parallel)

source('code/data_simulation.R')

##### Generate and Clean Data #####
set.seed(123)
ptm <- proc.time()
data.lst <- gen.pops(1000)
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

#all.data$time2d <- pmax(0, all.data$time2 - 1)
tdiff.lst <- lapply(unique(all.data$ind), function(u) {
  dat <- subset(all.data, ind == u)
  dat$time1d <- dat$time1
  dat$time2d <- dat$time2
  for (i in 1:5) {
    row <- dat[i,]
    if (i == 1) {
      if (row$time2 == 0) {
        row$time1d <- row$time1 - 1
        row$time2d <- row$time2
      }
      else {
        row$time2d <- row$time2 - 1
        row$time1d <- row$time1
      }
    }
    else {
      if (dat$time2[i-1] > 0) {
        #time2d <- time2 - 1
        if (row$time2 == 0) {
          row$time1d <- pmax(row$time1 - 1, 0)
          row$time2d <- 0
          row$time2 <- 1
        }
        else {
          row$time2d <- row$time2 - 1
          row$time2 <- row$time2 + 1 
          row$time1d <- row$time1
        }
      }
      else if (dat$time2[i-1] == 0) {
        row$time1 <- row$time1 + 1
        row$time1d <- row$time1 - 1
        if (row$time2 == 0) {
          row$time2d <- 0
          row$time1d <- pmax(row$time1d - 1, 0)
        }
        else {
          row$time2d <- row$time2 - 1
        }
      }
    }
    dat[i,] <- row
  }
  return(dat)
})
all.data <- bind_rows(tdiff.lst)
all.data <- aggregate(data = all.data, nobs ~ temp1 + temp2 + stage + 
                        time1 + time2 + time1d + time2d + prior.samp + block,
                      sum)

all.data <- subset(all.data, time2 <= 60)

write.csv(all.data, 'data/sim_data_reduced.csv')
write.csv(priors.df, 'data/sim_priors_reduced.csv')

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

basename <- "regniere_structured_reduced"
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
      dd2 <- with(dd, nList(temp1,time1,temp2,time2,time1d,time2d,stage,
                            nobs=as.integer(nobs)))
      dd2$t_block1 <- dd2$stage*7 + dd2$temp1/5 - 1
      dd2$t_block2 <- dd2$stage*7 + dd2$temp2/5 - 1
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
cl1 <- makeCluster(400)
clusterExport(cl1, c('all.data', 'parms', 'prior.samp', 'sv.lst'))
clusterEvalQ(cl1,{
  library(tidyverse)
  library(tmbstan)
  library(TMB)
  library(lme4)
  library(rstan)
  
  basename <- "regniere_structured_reduced"
  basename_loc <- paste0('code/', basename)
  cc <- compile(paste0(basename_loc, '.cpp'))
  try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
  dyn.load(dynlib(basename_loc))
  
  nList <- lme4:::namedList
  stages <- paste0('L', 2:6)
})

gr.lst1 <- parLapply(cl1, 1:1000, function(i) {
  
  dd <- subset(all.data, prior.samp == i)
  dd$stagename <- dd$stage
  dd$stage <- as.numeric(factor(dd$stage)) - 1
  
  dd2 <- with(dd, nList(temp1,time1,temp2,time2,time1d,time2d,stage,
                        nobs=as.integer(nobs)))
  dd2$t_block1 <- dd2$stage*7 + dd2$temp1/5 - 1
  dd2$t_block2 <- dd2$stage*7 + dd2$temp2/5 - 1
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
                   chains = 4, iter = 1500, warmup = 750,
                   control= list('max_treedepth' = 15, 'adapt_delta' = 0.99))
  post.df <- as.data.frame(stan1)
  write.csv(post.df, paste0('code/output/post_iter', i, '.csv'))
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

stopCluster(cl1)

diag.lst <- lapply(gr.lst, function(fit) {
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
write.csv(diag.df, 'code/output/diagnostics_reduced.csv')

post.lst <- lapply(gr.lst, as.data.frame)
rows <- sapply(post.lst, nrow)
w <- which(rows > 0)
post.lst.red <- lapply(post.lst, function(x) {
  if (nrow(x) > 0) {
    rows <- sample(1:nrow(x), 100)
    return(x[rows,])
  }
  else {return(NULL)}
})


curve.pars <- c('phi_rho', 'psi_rho', 'y0_rho', 
                'HA', 'TL', 'HL', 'TH', 'HH', 's_alpha')
priors.df.new <- priors.df
priors.df.new$prior.samp <- as.numeric(factor(priors.df.new$prior.samp))
ranks.lst <- lapply(w, function(x) {
  if (is.null(post.lst.red[[x]])) {return(NULL)}
  prior <- subset(priors.df.new, prior.samp == x)
  
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
    pr <- unique(prior.seps)[i]
    post <- post.seps[,i]
    rnk <- length(which(post < pr))/length(post)
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
write.csv(ranks.df, 'code/output/ranks_reduced.csv')

par(mfrow = c(3, 3)) 
for (i in 1:ncol(ranks.df)) {
  hist(as.data.frame(ranks.df)[,i], main = names(ranks.df)[i], breaks = 15)
}


##### Re-Generate Data #####
par.names <- sapply(names(post.lst[[1]]), function(x) {
  strsplit(x, split = '[.]')[[1]][1]
})
upn <- unique(par.names)[2:(length(unique(par.names))-1)]
ps.lst <- lapply(1:nrow(post.lst[[1]]), function(r) {
  row <- unlist(post.lst[[1]][r,])
  lap <- lapply(upn, function(x) {
    row[which(par.names == x)]
  })
  names(lap) <- upn
  lap <- append(lap, list(exp(lap$alpha)))
  names(lap)[length(lap)] <- 'zeta'
  lap$HL <- -lap$HL
  return(lap)
})

gen.pops2 <- function(N) {
  ps <- ps.lst[sample(1:length(ps.lst), N)]
  lst <- lapply(1:length(ps), function(i) {
    p <- ps.lst[[i]]
    pdf <- as.data.frame(p)
    pdf$stage <- rep(stages, 7)
    pdf$temp <- rep(seq(5, 35, by = 5), each = 5)
    t.lst <- lapply(seq(5, 35, by = 5), function(tmp) {
      #psub <- subset(pdf, temp == tmp)
      pdat <- pop.dev(tmp, pdf)
    })
    pd <- bind_rows(t.lst)
    pd$prior.samp <- i
    return(list('priors' = p, 'data' = pd))
  })
  return(lst)
}


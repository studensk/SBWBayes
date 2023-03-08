library(tmbstan)
library(TMB)
library(tidyverse)
library(lme4)
all.days.df <- read.csv('data/all_days_df.csv')
source('code/functions.R')
source('code/data_simulation.R')
nList <- lme4:::namedList

##### Stage Structured #####
basename <- "regniere_structured_reduced"
basename_loc <- paste0('code/', basename)
cc <- compile(paste0(basename_loc, '.cpp'))
try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
dyn.load(dynlib(basename_loc))
options(mc.cores = parallel::detectCores())

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
parms1 <- parms
diagnostics <- list()
chains <- 4
parm.lst <- prior.samp(chains)

p <- 'ON'
parms <- parms1
all.data <- subset(all.days.df, province == p & generation == 'F1')

all.data$block <- paste(all.data$stage, all.data$temp1, sep = '_')
lev.ord <- paste(rep(stages, each = 7), rep(seq(5, 35, by = 5), 5), sep = '_')
all.data$block <- factor(all.data$block)
all.data$block <- factor(all.data$block, levels = lev.ord)
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
all.data <- subset(all.data, select = c(nobs, temp1, temp2, stage, time1, 
                                        time2, time1d, time2d, block))

dd <- all.data
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

parm.lst <- lapply(parm.lst, function(x) {
  x$HL <- abs(x$HL)
  x$upsilon <- rep(0, 35)
  x$alpha <- rep(0.05, 5)
  return(x)
})
stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                 chains = chains, iter = 1250, warmup = 750,
                 control = list(max_treedepth = 15, adapt_delta = 0.999))
write.csv(as.data.frame(stan1), 'code/output/real_data_post.csv')

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

all.data$dd <- all.data$time1*all.data$temp1 + all.data$time2*all.data$temp2

## Clean model output for plotting and simulations

nms <- unique(sapply(names(post.df), function(x) {
  strsplit(x, split = '\\[')[[1]][1]
}))
nms <- nms[1:(length(nms)-1)]

posterior.lst <- lapply(1:nrow(post.df), function(r) {
  row <- post.df[r,]
  s_ups <- unlist(row[,grep('s_upsilon', names(row))])
  names(s_ups) <- NULL
  nm.lst <- lapply(nms, function(nm) {
    post <- unlist(row[,grep(paste0('^', nm), names(row))])
    names(post) <- NULL
    if (nm == 'HL') {post <- -abs(post)}
    else if (nm == 'upsilon') {
      post <- exp(qcauchy(pnorm(post, 0, 1), 0, rep(s_ups, 7)))
    }
    return(post)
  })
  names(nm.lst) <- nms
  nm.lst <- append(nm.lst, list('zeta' = exp(nm.lst$alpha)))
  return(nm.lst)
})

ppc <- makeCluster(40)
clusterEvalQ(ppc, {
  library(tidyverse)
  source('code/functions.R')
  source('code/data_simulation.R')
})
clusterExport(ppc, c('posterior.lst'))
gp.lst <- parLapply(ppc, 1:2000, function(x) {
  all.data <- gen.pops(N = 1, ps = list(posterior.lst[[x]]))[[1]]$data
  all.data$prior.samp <- x
  
  names(all.data) <- c('temp1', 'temp2', 'stage', 'time1', 'time2', 'nobs', 'cup', 'prior.samp')
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
  
  all.data <- aggregate(data = all.data, nobs ~ temp1 + temp2 + stage + 
                          time1 + time2 + time1d + time2d + prior.samp + block,
                        sum)
  write.csv(all.data, paste0('code/output/ppc/pop', x, '.csv'))
  return(all.data)
})
stopCluster(ppc)
gp.df <- bind_rows(gp.lst)
#gp <- gen.pops(N = 1, ps = posterior.lst)

gp.df$dd <- gp.df$temp1*gp.df$time1 + gp.df$temp2*gp.df$time2

gp.df.lst <- lapply(1:nrow(gp.df), function(r) {
  row <- gp.df[r,]
  rows <- dclone(gp.df[r,], n.clones = row$nobs)
  return(rows)
})
gp.df2 <- bind_rows(gp.df.lst)

ggplot(data = gp.df2) +
  geom_violin(aes(x = factor(temp1), y = dd), fill = 'grey50') +
  geom_point(data = all.data, aes(x = factor(temp1), y = dd, size = nobs), col = 'red') +
  scale_y_continuous(trans = 'log') +
  facet_wrap(vars(stage)) +
  theme_minimal()


#write.csv(gp.df, 'data/real_data_post.csv', row.names = FALSE)

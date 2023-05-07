library(tmbstan)
library(TMB)
library(tidyverse)
library(lme4)
all.days.df <- read.csv('data/all_days_df.csv')
source('code/functions.R')
source('code/noss_data_simulation.R')
nList <- lme4:::namedList

##### Stage Structured #####
basename <- "noss_regniere_structured_reduced"
basename_loc <- paste0('code/', basename)
cc <- compile(paste0(basename_loc, '.cpp'))
try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
dyn.load(dynlib(basename_loc))
options(mc.cores = parallel::detectCores())

parms <- list(
  'rho' = 0.045*(-1.05*(0:4)^2 + 4.22*(0:4) + 4.08),
  'HL' = 8.4,
  'HH' = 9.3,
  'HA' = 1.4, 
  'TL' = 285.90432,
  'TH' = 306.47530,
  's_eps' = rep(0.3, 5),
  's_upsilon' = rep(0.08, 5))
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
sm1 <- 0:4
dd2$k_vec <- -1.05*sm1^2 + 4.22*sm1 + 4.08

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

parm.lst <- lapply(parm.lst, function(x) {
  x$HL <- abs(x$HL)
  x$upsilon <- rep(0, 35)
  return(x)
})
stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                 chains = chains, iter = 1250, warmup = 750,
                 control = list(max_treedepth = 15, adapt_delta = 0.99))
write.csv(as.data.frame(stan1), 'code/output/noss_real_data_post.csv')

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
post.df.orig <- as.data.frame(stan1)

d.post.lst <- lapply(pnames, function(nm) {
  sub <- post.df.orig[,grep(paste0('^', nm), names(post.df.orig))]
  if (length(sub) == 5) {
    new <- unlist(lapply(1:nrow(sub), function(x) {unlist(sub[x,])}))
  }
  else if (length(sub) == 35) {
    lst <- lapply(1:7, function(i) {
      inds <- 7*(0:4) + i
      ups <- sub[,inds]
      df <- unlist(lapply(1:nrow(ups), function(x) {unlist(ups[x,])}))
      return(df)
    })
    new <- bind_cols(lst)
    names(new) <- paste0('upsilon.', 1:7, '.')
  }
  else {
    new <- rep(unlist(sub), each = 5)
  }
  return(as.data.frame(new))
})
d.post.df <- bind_cols(d.post.lst)
names(d.post.df)[1:8] <- pnames[1:8]
write.csv(d.post.df, 'data/noss_model_results.csv', row.names = FALSE)

## Power transformations
seps.par <- read.csv('code/final/s_eps_trans_pars.csv')
s_eps.lst <- lapply(1:5, function(x) {
  post <- post.df.orig[,paste0('s_eps[', x, ']')]
  row <- seps.par[x,]
  trans <- with(row, {exp(sc.rat*(powertrans(post, lambda) - mu.post) + mu.prior)})
  return(trans)
})
s_eps.df <- as.data.frame(s_eps.lst)
names(s_eps.df) <- paste0('s_eps[', 1:5, ']')
post.df <- post.df.orig
post.df[,grep('s_eps', names(post.df))] <- s_eps.df

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
      post <- exp(post*rep(s_ups, each = 7))
    }
    return(post)
  })
  names(nm.lst) <- nms
  return(nm.lst)
})

ppc <- makeCluster(40)
clusterEvalQ(ppc, {
  library(tidyverse)
  source('code/functions.R')
  source('code/noss_data_simulation.R')
})
clusterExport(ppc, c('posterior.lst'))
gp.lst <- parLapply(ppc, 1:1000, function(x) {
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
  write.csv(all.data, paste0('code/final/output/noss/ppc/pop', x, '.csv'))
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

all.data$dd <- all.data$temp1*all.data$time1 + all.data$temp2*all.data$time2

ggplot(data = gp.df2) +
  geom_violin(aes(x = factor(temp1), y = dd), fill = 'grey70', col = 'grey70', scale = 'width', bw = 10) +
  geom_point(data = all.data, aes(x = factor(temp1), y = dd, size = nobs), alpha = 0.5) +
  facet_wrap(vars(stage), scales = 'free') +
  labs(x = expression('Treatment Temperature ' (degree~C)), y = 'Degree Days', size = '# of Observations') +
  theme_minimal() +
  theme(legend.position = c(0.8, 0.3),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), size = 14),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 14),
        legend.title = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0), size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = 'bold'))

#write.csv(gp.df, 'data/real_data_post.csv', row.names = FALSE)

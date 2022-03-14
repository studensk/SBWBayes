library(tidyverse)
library(splines)
library(lme4)
library(reshape2)
grouped.data <- read.csv('data/grouped.csv')
source('code/functions.R')
nList <- lme4:::namedList
stages <- paste0("L", 2:6)

cv.results <- lapply(1:10, function(x) {
  obs.data <- subset(grouped.data, group_index != x)
  dd <- aggregate(data = obs.data, nobs ~ ., sum)
  
  dd$stagename <- dd$stage
  dd$stage <- as.numeric(factor(dd$stage)) - 1
  dd2 <- with(dd, nList(temp1,time1,temp2,time2,block,time2d,stage,
                        nobs=as.integer(nobs)))
  dd2$block <- as.numeric(factor(dd2$block)) - 1
  dd2.df <- as.data.frame(dd2)
  
  data <- read.csv(paste0('data/parametric_cv', x, '.csv'))
  post.df <- subset(data, select = -X)
  
  p.lst <- lapply(1:5, function(ind) {
    stg.ind <- (ind - 1) - 2
    
    alpha <- post.df[,paste0('alpha.', ind, '.')]
    sa2 <- -0.5*(post.df$s_alpha)^2
    alpha.mult <- exp(alpha*post.df$s_alpha + sa2)
    
    a <- (post.df$y0_rho/8)*(post.df$phi_rho + post.df$psi_rho - 2)
    b <- (post.df$y0_rho/4)*(post.df$psi_rho - post.df$phi_rho)
    c <- post.df$y0_rho
    rho <- a*stg.ind^2 + b*stg.ind + c
    rho25 <- rho*alpha.mult 
    q.pars <- post.df[,1:3]
    
    par.cols <- post.df[, grep(paste0('.', ind, '.'), names(post.df))[1:6]]
    pars <- as.data.frame(cbind(q.pars,
                                rho25, 
                                par.cols))
    names(pars) <- c('phi_rho', 'psi_rho', 'y0_rho', 'rho25', 
                     'HA', 'TL', 'HL', 'TH', 'HH', 's_eps')
    pars$HL <- -pars$HL
    pars$TA <- sapply(1:nrow(pars), function(r) {
      row <- pars[r,]
      ta <- with(row, ta.fun(HL, HH, TL, TH))
      return(ta)
    })
    
    stg.df <- subset(dd2.df, stage == ind - 1)
    stg.df <- stg.df[order(stg.df$temp1),]
    blocks <- unique(stg.df$block) + 1
    
    temps <- sort(unique(dd2$temp1))
    
    stage.df <- as.data.frame(cbind(pars, 'alpha' = alpha.mult, 
                                    's_alpha' = post.df$s_alpha))
    stage.df$stage <- stages[ind]
    return(stage.df)
  })
  p.df <- bind_rows(p.lst)
  p.df$group.index <- x
  return(p.df)
})
cv.df <- bind_rows(cv.results)
s.eps.mean <- aggregate(data = cv.df, s_eps ~ stage + group.index, mean)
#write.csv(cv.df, 'data/parametric_draws.csv', row.names = FALSE)

sim.pop.lst <- lapply(1:10, function(group) {
  stg.lst <- lapply(stages, function(stg) {
    data <- subset(cv.df, group.index == group & stage == stg)
    tempvec <- seq(5, 35, by = 5)
    delta.lst <- lapply(1:nrow(data), function(x) {rnorm(3500, 0, 1)})
    rate.est <- get_curves(data, temp = tempvec)
    temp.lst <- lapply(1:length(tempvec), function(x) {
      tmp <- tempvec[x]
      sub <- subset(rate.est, temp == tmp)
      all.rates <- unlist(lapply(1:length(delta.lst), function(ind) {
        sigma <- data$s_eps[ind]
        deltas <- exp(sigma*delta.lst[[ind]])
        rates <- 1/(deltas*sub$time[ind])
        return(rates)
      }))
      probs <- seq(0.001, 0.999, length.out = 999)
      q <- quantile(all.rates, probs)
      d <- density(all.rates)
      d$x <- pmax(d$x, 0)
      d$y <- d$y/sum(d$y)
      den <- approx(d$x, d$y, q)
      df <- data.frame('rate' = q, 'temp' = tempvec[x],
                       'quantile' = probs, 'density' = den$y)
      return(df)
    })
    random.pops <- bind_rows(temp.lst)
    random.pops$stage <- stg
    return(random.pops)
  })
  stg.df <- bind_rows(stg.lst)
  stg.df$group.index <- group
  return(stg.df)
}) 
sim.pop <- bind_rows(sim.pop.lst)
sim.pop$index <- paste0(sim.pop$stage, '_', sim.pop$group.index)
sim.pop$quantile <- round(sim.pop$quantile,3)

med.df <- subset(sim.pop, quantile == 0.5)
med.spline.lst <- lapply(unique(med.df$index), function(ind) {
  sub <- subset(med.df, index == ind)
  i.spline <- interpSpline(sub$temp, sub$rate)
  df <- as.data.frame(predict(i.spline, seq(5, 35, by = 0.1)))
  names(df) <- c('temp', 'rate')
  df$rate <- pmax(df$rate, 0)
  df$time <- 1/df$rate
  df$stage <- unique(sub$stage)
  df$group.index <- unique(sub$group.index)
  return(df)
})
med.spline <- bind_rows(med.spline.lst)

comp.lst <- lapply(1:10, function(group) {
  estimates <- subset(sim.pop, group.index == group)
  ag.obs <- subset(grouped.data, group_index == group)
  #ag.obs <- aggregate(data = observations, nobs ~ ., sum)
  days <- lapply(1:nrow(ag.obs), function(r) {
    row <- ag.obs[r,]
    stage.r <- row$stage
    temp.r1 <- row$temp1
    temp.r2 <- row$temp2
    rate1 <- subset(estimates, temp == temp.r1 & 
                      stage == stage.r & quantile %in% c(0.050, 0.5, 0.950))$rate
    rate2 <- subset(estimates, temp == temp.r2 & 
                      stage == stage.r & quantile %in% c(0.050, 0.5, 0.950))$rate
    dev <- row$time1*rate1
    rem <- 1-pmin(dev, 1)
    ndays <- sort(rem/rate2)
    names(ndays) <- c('q5', 'q50', 'q95')
    return(ndays)
  })
  days.df <- bind_rows(days)
  ag.obs <- bind_cols(ag.obs, days.df)
  return(ag.obs)
})
comp.df <- bind_rows(comp.lst)

overlap <- function(l1, u1, l2, u2) {
  !(l1 > u2 | u1 < l2)
}

in.predint <- mapply(function(obs.l, obs.u, lb, ub) {
  return(ifelse(overlap(obs.l, obs.u, lb, ub), 1, 0))
}, comp.df$time2d, comp.df$time2, comp.df$q5, comp.df$q95)

comp.df$in.predint <- in.predint

ag.group <- aggregate(data = comp.df, in.predint ~ group_index, mean)
ag.stage.temp <- aggregate(data = comp.df, in.predint ~ temp1 + stage, mean)

write.csv(comp.df, 'data/parametric_cv_results.csv')

library(tidyverse)
library(splines)
library(lme4)
library(reshape2)
library(parallel)

grouped.data <- read.csv('data/grouped.csv')
cv.df <- read.csv('data/group_df.csv')

cl <- makeCluster(10)
clusterExport(cl, varlist = c('cv.df'))
clusterEvalQ(cl,{
  library(tidyverse)
  library(splines)
  library(lme4)
  library(reshape2)
})

sim.pop.lst <- parLapply(cl, 1:10, function(grp) {
  source('code/functions.R')
  stages <- paste0("L", 2:6)
  stg.lst <- lapply(stages, function(stg) {
    data <- subset(cv.df, group == grp & stage == stg)
    tempvec <- seq(5, 35, by = 5)
    delta.vec <- rnorm(3500, 0, 1)
    rate.est <- get_curves(data, temp = tempvec)
    temp.lst <- lapply(1:length(tempvec), function(x) {
      tmp <- tempvec[x]
      sub <- subset(rate.est, temp == tmp)
      all.rates <- unlist(lapply(1:nrow(data), function(ind) {
        sigma <- data$s_eps[ind]
        deltas <- exp(sigma*delta.vec)
        rates <- 1/(delta.vec*sub$time[ind])
        return(rates)
      }))
      all.rates <- pmax(all.rates, 0)
      probs <- seq(0.001, 0.999, length.out = 999)
      q <- quantile(all.rates, probs, na.rm = TRUE)
      d <- density(all.rates, na.rm = TRUE)
      d$x <- pmax(d$x, 0)
      d$y <- d$y/sum(d$y, na.rm = TRUE)
      den <- approx(d$x, d$y, q)
      den$y[is.na(den$y)] <- 0
      df <- data.frame('rate' = q, 'temp' = tempvec[x],
                       'quantile' = probs, 'density' = den$y,
                       'nna' = length(which(is.na(sub$rate)))/nrow(sub))
      write.csv(df, paste0('data/cv_mle/ddf', tempvec[x], '_', stg, '_', grp, '.csv'))
      return(df)
    })
    random.pops <- bind_rows(temp.lst)
    random.pops$stage <- stg
    return(random.pops)
  })
  stg.df <- bind_rows(stg.lst)
  stg.df$group.index <- grp
  return(stg.df)
}) 
stopCluster(cl)
sim.pop <- bind_rows(sim.pop.lst)
sim.pop$index <- paste0(sim.pop$stage, '_', sim.pop$group.index)
sim.pop$quantile <- round(sim.pop$quantile,3)

#write.csv(sim.pop, 'data/cv_mle/sim_pop.csv', row.names = FALSE)

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

comp.lst <- lapply(1:10, function(grp) {
  print(grp)
  estimates <- subset(sim.pop, group.index == grp)
  ag.obs <- subset(grouped.data, group_index == grp)
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

write.csv(comp.df, 'data/mle_parametric_cv_results.csv')

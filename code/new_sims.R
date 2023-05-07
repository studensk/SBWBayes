library(tidyverse)
library(imputeTS)
library(smooth)
all.days.df <- read.csv('data/all_days_df.csv')
source('code/functions.R')

## Import weather data
rename <- plyr::rename
weather <- weathercan::weather_dl
add_weather <- weathercan::weather_interp
stations_search <- weathercan::stations_search
stations <- weathercan::stations

on.weather <- as.data.frame(weather(50460, '2018-01-01', '2020-12-31', interval = 'hour', 
                            string_as = NULL))[,c('prov', 'date', 'year', 'month', 
                                                  'day', 'hour', 'temp')]
on2020 <- subset(on.weather, year == 2020 & between(as.numeric(month), 4, 9))
on2020$temp <- sapply(1:nrow(on2020), function(r) {
  temp <- on2020$temp[r]
  new.temp <- ifelse(is.na(temp), mean(on2020$temp[r + c(-1, 1)]), temp)
  return(new.temp)
})

temps.on <- on2020$temp
temps.unique <- sort(unique(temps.on))

write.csv(on2020, 'data/on2020.csv', row.names = FALSE)

stages <- c('L2', 'L3', 'L4', 'L5', 'L6')
prov <- 'ON'

## Import results of stage-structured model runs
struc.df <- read.csv('data/noss_model_results.csv')
struc.df$index <- rep(1:(nrow(struc.df)/5), each = 5)
struc.df$HL <- -abs(struc.df$HL)

##### Simulations #####

## Function to simulate daily development of a population based on 
##   weather inputs
pop.dev <- function(temps, param.df, size = 100, all.dev = FALSE) {
  temps.unique <- sort(unique(temps))
  dev.lst <- lapply(stages, function(st) {
    all.temps <- seq(0, 40, by = 0.1)
    params <- subset(param.df, stage == st)
    pdf.full$rho25 <- pdf.full$rho
    rates <- get_curves(pdf.full, temp = all.temps, spline = TRUE)$rate
    s.df <- data.frame(all.temps, rates)
    names(s.df) <- c('temp', 'rate')
    
    rates <- approx(s.df$temp, s.df$rate, xout = pmax(temps.unique, 0))$y
    df <- data.frame('temps' = temps.unique, rates, 'stage' = st)
    return(df)
  })
  dev.df <- bind_rows(dev.lst)
  l <- lapply(1:size, function(ind) {
    hours <- vector()
    eps <- rnorm(5, 0, param.df$s_eps)
    deltas <- exp(eps)
    #if (ind == 1) {deltas <- rep(1, 5)}
    names(deltas) <- stages
    dev <- 0
    i <- 1
    dev.vec <- vector()
    for (s in stages) {
      rdf <- subset(dev.df, stage == s)
      rdf$rates <- rdf$rates*deltas[s]
      while(dev < 1) {
        temp <- temps[i]
        add.dev <- rdf$rates[rdf$temps == temp]/24
        dev <- dev + add.dev
        dev.vec <- c(dev.vec, dev)
        if (i == length(temps)) {break}
        i <- i + 1}
      if (s != 'L6') {
        next.prop <- (dev - 1)/add.dev
        next.stage <- stages[which(stages == s) + 1]
        r <- subset(dev.df, stage == next.stage)
        r$rates <- r$rates*deltas[next.stage]
        add.dev <- next.prop*r$rates[r$temps == temp]/24
      }
      dev <- add.dev
      hrs <- i - 1
      hours <- c(hours, hrs)
    }
    days <- hours/24
    names(days) <- stages
    if (all.dev) {
      ddv <- diff(dev.vec)
      ddv <- ddv[ddv >= 0]
      cddv <- cumsum(ddv)
      df <- data.frame('dev' = cddv, 'time' = 1:length(cddv), 
                       'individual' = ind)
      df <- subset(df, (time %% 24) == 0)
      df <- rbind(df, data.frame('dev' = 5, 'time' = max(df$time + 24),
                                 'individual' = ind))
      return(df)
    }
    return(days)
  })
  if (all.dev) {return(bind_rows(l))}
  df <- bind_rows(l)
  return(df)
}

## Wrapper to run pop.dev over multiple populations
sim.multi <- function(param.lst, temps, size = 100, all.dev = FALSE) {
  sim.lst <- lapply(1:length(param.lst), function(ind) {
    param <- param.lst[[ind]]
    pd <- pop.dev(temps, param, size, all.dev)
    pd$index <- ind
    return(pd)
  })
  sim.df <- bind_rows(sim.lst)
  return(sim.df)
}

## Simulate 1000 individuals from 1000 posterior draws
struc.lst <- lapply(unique(struc.df$index), function(x) {
  subset(struc.df, index == x)})
ssd <- sim.multi(struc.lst, temps.on, size = 100, all.dev = TRUE)

ssd$group <- mapply(function(i1, i2) {paste0(i1, '_', i2)}, 
                    ssd$index, ssd$individual)
ssd$date <- on2020[ssd$time,'date']

## Obtain quantiles for daily development
q.lst <- lapply(unique(ssd$date), function(d) {
  print(d)
  sub <- subset(ssd, date == d)
  sub.dev <- sub$dev
  len <- 100000 - length(sub.dev)
  if (len > 0)  {sub.dev <- c(sub.dev, rep(5, len))}
  quants.low <- c(0.005, 0.025, 0.05, 0.25, 0.5)
  quants.high <- 1 - quants.low
  q.low <- quantile(sub.dev, quants.low)
  q.high <- quantile(sub.dev, quants.high)
  df <- data.frame('dev.low' = q.low,
                   'dev.high' = q.high,
                   'base' = quants.low, 
                   'date' = d)
  return(df)
})
q.df <- bind_rows(q.lst)
pred.int <- q.df

write.csv(pred.int, 'data/noss_ribbon_df.csv', row.names = FALSE)

gc.lst <- lapply(unique(struc.df$index), function(ind) {
  sub <- subset(struc.df, index == ind)
  gc <- get_curves(sub, spline = TRUE)
  gc$stage <- sapply(gc$index, function(x) {stages[x]})
  return(gc)
})
gc.df <- bind_rows(gc.lst)

quants.low <- c(0.005, 0.025, 0.05, 0.25, 0.5)
quants.high <- 1 - quants.low
quants <- sort(unique(c(quants.low, quants.high)))
ag.gc.df <- aggregate(data = gc.df, rate ~ temp + stage, 
                      function(x) {
                          q <- quantile(x, quants)
                          return(q)
                        })
ag.gc.df <- bind_cols(ag.gc.df[,1:2], as.data.frame(ag.gc.df[,3]))
names(ag.gc.df)[3:length(ag.gc.df)] <- 
  c('rate.005', 'rate.025', 'rate.05', 'rate.25', 'rate.5', 
    'rate.75', 'rate.95', 'rate.975', 'rate.995')
write.csv(ag.gc.df, 'data/dev_curves.csv', row.names = FALSE)

## Obtain splines from posterior draws

ssd$stage <- sapply(ssd$dev, function(d) {
  ind <- floor(d) + 1
  return(stages[ind])
})
ind.dev  <- floor(ssd$dev) + 1
ssd$stage <- stages[ind.dev]

## Generate tiles for upper panel of figure 5
ssd$stage <- factor(ssd$stage)
tables <- lapply(unique(ssd$date), function(d) {table(subset(ssd, date == d)$stage)})
count.df <- bind_rows(tables)
count.mat <- apply(count.df, 1, function(x) {x/sum(x)})
count.df <- as.data.frame(t(count.mat))
ag.ss <- aggregate(data = gc.df, rate ~ temp + stage, median)
dates <- sort(unique(ssd$date))
lst <- lapply(1:length(dates), function(d) {
  prop <- as.numeric(count.df[d,])
  s <- sapply(unique(ag.ss$temp), function(tmp) {
    sub <- as.numeric(subset(ag.ss, temp == tmp)$rate)
    sp <- weighted.mean(sub, prop)
    return(sp)
  })
  return(s)
})
lst[[length(lst)]] <- lst[[length(lst) - 1]]
lst2 <- lapply(1:length(lst), function(i) {
  date <- dates[i]
  x <- lst[[i]]
  return(data.frame('dev' = x, 'temp' = unique(ag.ss$temp), 'date' = date))
})
df2 <- bind_rows(lst2)
l <- lapply(seq(-10, 4.9, by = 0.1), function(temp) {
  data.frame('dev' = 0, temp, date = dates)
})
dl <- bind_rows(l)
df3 <- bind_rows(df2, dl)
df3$temp <- round(df3$temp, 1)

df3 <- write.csv(df3, 'data/df_counts.csv', row.names = FALSE)


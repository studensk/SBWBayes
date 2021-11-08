library(tidyverse)
library(weathercan)
library(plyr)
library(dplyr)
library(splines)
library(lubridate)
library(RColorBrewer)
library(imputeTS)
library(colorspace)
library(reshape2)
library(ggpubr)
library(smooth)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
source('code/items.R')
all.days.df <- read.csv('data/all_days_df.csv')
source('code/functions.R')
source('code/test_stan.R')

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

stages <- c('L2', 'L3', 'L4', 'L5', 'L6')
prov <- 'ON'

## Import results of stage-structured model runs
struc.df <- read.csv('code/output/simulations/re_structured_ON_NEW.csv')
struc.df$index <- rep(1:(nrow(struc.df)/5), 5)
names(struc.df)[grep('upsilon.', names(struc.df))] <-
  paste0('upsilon[', seq(5, 35, by = 5), ']')

##### Simulations #####

## Function to simulate daily development of a population based on 
##   weather inputs
pop.dev <- function(temps, param.df, size = 100, all.dev = FALSE) {
  temps.unique <- sort(unique(temps))
  dev.lst <- lapply(stages, function(stage) {
    treatments <- seq(5, 35, by = 5)
    all.temps <- seq(0, 40, by = 0.1)
    params <- param.df[which(stages == stage),]
    ups <- params[,grep('upsilon', names(params))]
    ups <- subset(ups, select = -s_upsilon)
    u.ups <- pnorm(unlist(ups))
    c.ups <- qcauchy(u.ups, 1, params$s_upsilon)
    times <- with(params, calc_pred3(treatments, rho25, HA, TL, HL, TH, HH, TA))
    times.new <- c.ups*times
    rates <- 1/times.new
    i.spline <- interpSpline(treatments, rates)
    s.df <- as.data.frame(predict(i.spline, all.temps))
    names(s.df) <- c('temp', 'rate')
    first <- subset(s.df, temp <= 20)
    last <- subset(s.df, temp > 20)
    if (length(which(first$rate <= 0) > 0)) {
      f.w <- max(which(first$rate <= 0))
      first[1:f.w, 'rate'] <- 0
    }
    if (length(which(last$rate <= 0) > 0)) {
      l.w <- min(which(last$rate <= 0))
      last[l.w:nrow(last), 'rate'] <- 0
    }
    s.df <- bind_rows(first, last)
    
    rates <- approx(s.df$temp, s.df$rate, xout = pmax(temps.unique, 0))$y
    df <- data.frame('temps' = temps.unique, rates, stage)
    return(df)
  })
  dev.df <- bind_rows(dev.lst)
  l <- lapply(1:size, function(ind) {
    hours <- vector()
    eps <- rnorm(5, 0, param.df$s_eps)
    deltas <- exp(eps)
    if (ind == 1) {deltas <- rep(1, 5)}
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
ssd <- sim.multi(lapply(sample(unique(struc.df$index), 1000), function(x) {
  subset(struc.df, index == x)}), temps.on, size = 1000, all.dev = TRUE)

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

#pred.int <- read.csv('code/output/simulations/ribbon_df.csv')

## Aggregate temperature data
on2020.cut <- subset(on2020, date %in% unique(ssd$date))
maxtemp <- aggregate(data = on2020.cut, temp ~ date, max)$temp
mintemp <- aggregate(data = on2020.cut, temp ~ date, min)$temp
meantemp <- aggregate(data = on2020.cut, temp ~ date, mean)$temp
ontemp <- data.frame('date' = unique(on2020.cut$date), maxtemp, mintemp, meantemp)

## Obtain splines from posterior draws
struc.spl.lst <- lapply(stages, function(st) {
  sub <- subset(struc.df, stage == st)
  eps <- sub$s_eps
  q1 <- qlnorm(0.01, mean = 0, sd = eps)
  q99 <- qlnorm(0.99, mean = 0, sd = eps)
  sp <- get_curves(sub, temp = tempvec, spline = TRUE)
  sp$q1 <- sp$rate*q1
  sp$q99 <- sp$rate*q99
  sp$stage <- st
  return(sp)
})
struc.spl <- bind_rows(struc.spl.lst)
rm(struc.spl.lst)

## Obtain parametric curves from posterior draws
struc.crv.lst <- lapply(stages, function(st) {
  sub <- subset(struc.df, stage == st)
  eps <- sub$s_eps
  q1 <- qlnorm(0.01, mean = 0, sd = eps)
  q99 <- qlnorm(0.99, mean = 0, sd = eps)
  sp <- get_curves(sub, temp = tempvec)
  sp$q1 <- sp$rate*q1
  sp$q99 <- sp$rate*q99
  sp$stage <- st
  return(sp)
})
struc.crv <- bind_rows(struc.crv.lst)
rm(struc.crv.lst)

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
ag.ss <- aggregate(data = struc.spl, rate ~ temp + stage, median)
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


##### Plotting ####
## Figure 5, bottom
daily <- ggplot(q.df) + 
  geom_ribbon(aes(x = date, ymin = dev.low, ymax = dev.high, 
                  group = base, fill = base)) +
  scale_fill_continuous(low = 'grey90', high = 'black') + 
  scale_y_continuous(minor_breaks = seq(0, 5, by = 0.25), breaks=c(0:5),
                     labels=c("L2", 'L3', 'L4', 'L5', 'L6', 'Pupa')) +
  scale_x_continuous(minor_breaks = seq(min(unique(ssd$date)), max(unique(ssd$date)), 
                                        by = 'day'), 
                     breaks = seq(min(unique(ssd$date)), max(unique(ssd$date)), 
                                  by = 'week')) +
  geom_line(data = subset(q.df, base == 0.5), aes(x = date, y = dev.low)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), 
                                    size = 14),
        legend.position = 'bottom',
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13)) +
  labs(y = 'Stage', x = 'Date', fill = 'Lower Quantile')

## Figure 5, top
temps <- ggplot(data = df3) + 
  geom_tile(aes(x = date, y = temp, fill = dev)) + 
  scale_fill_continuous(low = 'white', high = 'black') + 
  theme_minimal() + 
  geom_line(data = ontemp, aes(x = date, y = meantemp)) + 
  labs(y = expression('Temperature ' (degree~C)), x = '', 
       fill = 'Expected Development Rate   ') + 
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    size = 14),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), 
                                    size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))

ggarrange(temps, daily, ncol = 1, nrow = 2, heights = c(1, 2), labels = c('(a)', '(b)'))
on.data <- subset(all.days.df, province == 'ON' & generation == 'F1')

paired.cols <- brewer.pal(12, 'Paired')
line.cols <- paired.cols[seq(2,length(paired.cols), by = 2)]
ribbon.cols <- paired.cols[seq(1, length(paired.cols) - 1, by = 2)]

ss.trt <- subset(struc.spl, temp %in% seq(5, 35, by = 5))

## Figure 3
violins <- ggplot(struc.crv) +
  geom_line(aes(x = temp, y = rate, group = index),
            alpha = 0.6, col = 'darkgrey') +
  facet_wrap(vars(stage))+
  theme_minimal() +
  geom_violin(data = ss.trt, 
              aes(x = temp, y = rate, group = temp), 
              width = 1.5,
              scale = 'width',
              fill = 'black') +
  labs(x = expression('Rearing Temperature ' (degree~C)),
       y = 'Development Rate') +
  theme(strip.text = element_text(face = 'bold'),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    size = 14),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), 
                                    size = 14),
        strip.text.x = element_text(size = 15, face = 'bold',
                                    margin = margin(t = 0, r = 0, b = 14, l = 0)),
        axis.text = element_text(size = 12),
        panel.spacing = unit(1, 'lines'))

## Figure 4
ribbons <- ggplot(data = pred.int) +
  geom_ribbon(aes(x = temp, ymin = rate.005, ymax = rate.995),
              alpha = 0.1) +
  geom_ribbon(aes(x = temp, ymin = rate.025, ymax = rate.975), 
              alpha = 0.25) +
  geom_ribbon(aes(x = temp, ymin = rate.25, ymax = rate.75),
              alpha = 0.5) +
  geom_line(aes(x = temp, y = rate.5)) +
  geom_point(data = subset(on.data, r.est1 == r.est2), 
             aes(x = temp1, y = r.est1, size = nobs), alpha = 0.5) + 
  geom_segment(data = on.data, 
               aes(x = temp1, y = r.est1, 
                   xend = temp1, yend = r.est2, size = nobs), alpha = 0.5) + 
  geom_segment(data = subset(on.data, !is.finite(r.est2)), 
               aes(x = temp1 + 1.25, xend = temp1 + 1.25, y = 1.15, yend = 1.25),
               arrow = arrow(length = unit(0.1, 'cm'))) +
  labs(size = '# of Observations',
       x = expression('Rearing Temperature ' (degree~C)),
       y = 'Development Rate') +
  facet_wrap(vars(stage)) +
  theme_minimal() +
  theme(legend.position = 'top',
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(face = 'bold'),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0),
                                    size = 17),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),
                                    size = 17),
        strip.text.x = element_text(size = 18, face = 'bold',
                                    margin = margin(t = 0, r = 0, b = 14, l = 0)),
        axis.text = element_text(size = 15),
        panel.spacing = unit(1, 'lines')) +
  ylim(0, 1.25)

## Melt & clean data for Figure 2
sdf <- struc.df
names(sdf)[1:12] <- c('phi', 'psi', 'y[0]', 'rho', 'H[A]', 'T[L]', 'H[L]',
                'T[H]', 'H[H]', 'sigma[epsilon]', 'sigma[upsilon]', 'T[A]')
names(sdf)[21] <- 'sigma[alpha]'
m.sdf <- melt(data = sdf, 
              id.vars = 'stage',
              measure.vars = c('rho', 'alpha',
                               'sigma[epsilon]', 'sigma[upsilon]'))

m.sdf$variable <- factor(m.sdf$variable, 
                         levels = c('rho', 'alpha', 
                                    'sigma[epsilon]', 'sigma[upsilon]'))

ups.df <- sdf[,grep('upsilon', names(sdf))]
ups.df$stage <- sdf$stage
sig <- ups.df[,1]
ups <- ups.df[,2:8]
ups.lst <- lapply(1:length(sig), function(ind) {
  sigma <- sig[ind]
  upsilon <- as.numeric(ups[ind,])
  u.upsilon <- pnorm(upsilon, 0, 1)
  c.upsilon <- qcauchy(u.upsilon, 1, sigma)
  cl <- as.list(c.upsilon)
  df <- data.frame(cl, 'stage' = ups.df$stage[ind])
  names(df)[1:7] <- colnames(ups)
  return(df)
})
ups.df <- bind_rows(ups.lst)

ups.df <- melt(data = ups.df,
               id.vars = 'stage', 
               measure.vars = c('upsilon[5]', 'upsilon[10]', 'upsilon[15]',
                                'upsilon[20]', 'upsilon[25]', 'upsilon[30]',
                                'upsilon[35]'))

levels(m.sdf$variable)[2] <- 'zeta'

ps.lst <- lapply(prior.samp(70000), function(lst) {
  as.data.frame(lst[c('phi_rho', 'psi_rho', 'y0_rho',
                      's_eps', 's_upsilon', 's_alpha')])[1,]
})
ps.df <- bind_rows(ps.lst)
names(ps.df)[4:6] <- c('sigma[epsilon]', 'sigma[upsilon]', 'sigma[alpha]')
ps.df$zeta <- sapply(ps.df$`sigma[alpha]`, function(sigma) {
  alpha <- rnorm(1, -(sigma^2)/2, sigma)
  zeta <- exp(alpha)
  return(zeta)
})
ps.df$stage <- ''
rho <- as.vector(t(mapply(function(phi, psi, y0) {
  quadratic(-2:2, phi, psi, y0)
}, ps.df$phi_rho, ps.df$psi_rho, ps.df$y0_rho)))
rho.df <- data.frame('stage' = rep(stages, each = length(rho)/5),
                     'variable' = 'rho',
                     'value' = rho)

m.ps.df <- melt(data = ps.df, 
                id.vars = 'stage',
                measure.vars = c('sigma[epsilon]', 'sigma[upsilon]', 'zeta'))
m.ps.df <- bind_rows(m.ps.df, rho.df)

m.sdf$Distribution <- "Posterior"
m.ps.df$Distribution <- "Prior"
m.sdf2 <- bind_rows(m.sdf, m.ps.df)

l.short <- c('rho', 'zeta',
             'sigma[epsilon]',
             'sigma[upsilon]')

l.long <- c('`Multiplicative Intercept` ~~~ (rho)', 
            '`Bias Reduction Multiplier` ~~~ (zeta)',
            '`Scale of Individual Variation` ~~~ (sigma[epsilon])',
            '`Scale of Bias Reduction` ~~~ (sigma[upsilon])')

m.sdf2$variable <- factor(m.sdf2$variable, levels = l.short)
m.sdf2$variable <- mapvalues(m.sdf2$variable, from = l.short, to = l.long)
m.sdf2$Distribution <- factor(factor(m.sdf2$Distribution),
                              levels = c('Prior', 'Posterior'))

## Figure 2
ggplot(data = m.sdf2) + 
  geom_boxplot(aes(x = stage, y = value, fill = Distribution)) + 
  scale_fill_manual(values = c('grey90', 'grey50')) +
  facet_wrap(vars(variable), scales = 'free', labeller = 'label_parsed', nrow = 2) + 
  theme_minimal() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),
                                    size = 17),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), 
                                    size = 17),
        strip.text.x = element_text(size = 16, face = 'bold',
                                    margin = margin(t = 0, r = 0, b = 14, l = 0)),
        axis.text = element_text(size = 13),
        panel.spacing = unit(1, 'lines'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  labs(x = 'Stage', y = 'Value')

onl2 <- subset(all.days.df, province == 'ON' & stage == 'L2')  
temp.reps <- unlist(lapply(1:nrow(onl2), function(r) {
  row <- onl2[r,]
  temp <- row$temp1
  count <- row$nobs
  return(rep(temp, count))
}))
tr.df <- data.frame('temp' = temp.reps, 
                    'Treatment' = sapply(temp.reps, function(tr) {
                      ifelse(tr %in% c(15, 20, 25), 'Constant', 'Transfer')
                    }))
## Figure 1
ggplot(tr.df, aes(x = factor(temp), fill = Treatment)) + 
  scale_fill_manual(values = c('grey50', 'grey80')) +
  geom_bar(col = 'black') + 
  theme_minimal() + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),
                                    size = 17),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), 
                                    size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)) +
  labs(x = expression('Rearing Temperature ' (degree~C)), y = 'Count') +
  geom_text(stat='count', aes(label=..count..), vjust=-1, size = 6) +  
  ylim(0, 250)


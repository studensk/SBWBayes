library(tidyverse)
library(weathercan)
library(plyr)
library(dplyr)
library(splines)
library(lubridate)
library(ggpubr)
library(reshape2)

all.days.df <- read.csv('data/all_days_df.csv')
source('code/functions.R')
source('code/test_stan.R')

on2020 <- read.csv('data/on2020.csv')
pred.int <- read.csv('data/ribbon_df.csv')
df3 <- read.csv('data/df_counts.csv')
gc.df <- read.csv('data/dev_curves.csv')
struc.df <- read.csv('data/model_results.csv')
sp.struc.df <- read.csv('data/re_structured_ON_NEW.csv')

on2020$date <- as.Date(on2020$date)
pred.int$date <- as.Date(pred.int$date)
df3$date <- as.Date(df3$date)

stages <- paste0('L', 2:6)

## Aggregate temperature data
on2020.cut <- subset(on2020, date %in% unique(df3$date))
maxtemp <- aggregate(data = on2020.cut, temp ~ date, max)$temp
mintemp <- aggregate(data = on2020.cut, temp ~ date, min)$temp
meantemp <- aggregate(data = on2020.cut, temp ~ date, mean)$temp
ontemp <- data.frame('date' = unique(on2020.cut$date), maxtemp, mintemp, meantemp)

prior.samp <- function(chains) {
  l <- lapply(1:chains, function(x) {
    lst <- list(
      'phi_rho' = rbeta(1, 4, 4),
      'psi_rho' = rbeta(1, 4, 4),
      'y0_rho' = rgamma(1, 10, scale = 0.1),
      'HA' = rgamma(5, 5, scale = 0.2),
      'TL' = rnorm(5, 284, 2),
      'HL' = rgamma(5, 6, scale = 2),
      'TH' = rnorm(5, 304, 2),
      'HH' = rgamma(5, 10, scale = 3),
      's_eps' = exp(rnorm(5, -1.5, 0.1)),
      's_upsilon' = exp(rnorm(1, -2.5, 0.05)),
      's_alpha' = exp(rnorm(1, -1.5, 0.3)),
      'upsilon' = 0,
      'alpha' = 0.05
    )
    return(lst)
  })
  return(l)
}

## Figure 5, bottom
daily <- ggplot(data = pred.int) + 
  geom_ribbon(aes(x = date, ymin = dev.low, ymax = dev.high, 
                  group = base, fill = base)) +
  scale_fill_continuous(low = 'grey90', high = 'black') + 
  scale_y_continuous(minor_breaks = seq(0, 5, by = 0.25), breaks=c(0:5),
                     labels=c("L2", 'L3', 'L4', 'L5', 'L6', 'Pupa')) +
  scale_x_continuous(minor_breaks = seq(min(unique(df3$date)), max(unique(df3$date)), 
                                        by = 'day'), 
                     breaks = seq(min(unique(df3$date)), max(unique(df3$date)), 
                                  by = 'week')) +
  geom_line(data = subset(pred.int, base == 0.5), aes(x = date, y = dev.low)) +
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

## Figure 4
on.data$r.est2.trim <- pmin(on.data$r.est2, 1.25)
ribbons <- ggplot(data = gc.df) +
  geom_ribbon(aes(x = temp, ymin = rate.005, ymax = rate.995),
              alpha = 0.1) +
  geom_ribbon(aes(x = temp, ymin = rate.025, ymax = rate.975), 
              alpha = 0.25) +
  geom_ribbon(aes(x = temp, ymin = rate.25, ymax = rate.75),
              alpha = 0.5) +
  geom_line(aes(x = temp, y = rate.5)) +
  geom_point(data = subset(on.data, r.est1 == r.est2), 
             aes(x = temp1, y = r.est1, 
                 size = nobs), alpha = 0.5) + 
  geom_segment(data = on.data, 
               aes(x = temp1, y = r.est1, 
                   xend = temp1, yend = r.est2.trim,
                   size = nobs), alpha = 0.5) + 
  geom_segment(data = subset(on.data, !is.finite(r.est2)), 
               aes(x = temp1 + 1.4, xend = temp1 + 1.4, y = 1.15, yend = 1.25),
               arrow = arrow(length = unit(0.1, 'cm'))) +
  labs(size = '# of Observations',
       x = expression('Rearing Temperature ' (degree~C)),
       y = 'Development Rate') +
  facet_wrap(vars(stage), scales = 'free') +
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
        panel.spacing = unit(1, 'lines'))# +
  #ylim(0, 1.25)

## Melt data for Figure 2
sdf <- sp.struc.df
names(sdf)[1:12] <- c('phi', 'psi', 'y[0]', 'rho', 'H[A]', 'T[L]', 'H[L]',
                      'T[H]', 'H[H]', 'sigma[epsilon]', 'sigma[upsilon]', 'T[A]')
names(sdf)[21] <- 'sigma[alpha]'
m.sdf <- melt(data = sdf, 
              id.vars = 'stage',
              measure.vars = c('rho', 'alpha',
                               'sigma[epsilon]', 
                               'sigma[upsilon]'))

m.sdf$variable <- factor(m.sdf$variable, 
                         levels = c('rho', 'alpha', 
                                    'sigma[epsilon]', 
                                    'sigma[upsilon]'))

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
            '`Non-Parametric Adjustment` ~~~ (zeta)',
            '`Scale of Individual Variation` ~~~ (sigma[epsilon])',
            '`Scale of Adjustment` ~~~ (sigma[upsilon])')

m.sdf2$variable <- factor(factor(m.sdf2$variable), levels = l.short)
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

## Supplementary Materials 
curve.parms <- c('HL', 'HA','HH', 'TL','TA','TH')
c.sdf <- subset(struc.df, select = c(curve.parms, 'stage'))
names(c.sdf)[1:length(curve.parms)] <- c('H[L]', 'H[A]', 'H[H]',
                                         'T[L]', 'T[A]', 'T[H]')
m.sdf <- melt(c.sdf, id.vars = 'stage', variable.name = 'parameter')

ggplot(data = m.sdf) +
  geom_violin(aes(x = stage, y = value), draw_quantiles = 0.5,
              fill = 'black', alpha = 0.4, scale = 'width') +
  facet_wrap(vars(parameter), labeller = 'label_parsed', scales = 'free') +
  theme_minimal() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),
                                    size = 15),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), 
                                    size = 15),
        strip.text.x = element_text(size = 14, face = 'bold',
                                    margin = margin(t = 0, r = 0, b = 14, l = 0)),
        axis.text = element_text(size = 11),
        panel.spacing = unit(1.5, 'lines')) +
  labs(x = 'Stage', y = 'Posterior Distribution')
  

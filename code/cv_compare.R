library(reshape2)
library(tidyverse)

preds <- read.csv('data/cv_results.csv')
preds$model <- 'Semi-Parametric'
p.preds <- read.csv('data/parametric_cv_results.csv')
p.preds$model <- 'Parametric'
mle.preds <- read.csv('data/mle_parametric_cv_results.csv')
mle.preds$model <- 'MLE'

all.preds <- bind_rows(preds, p.preds, mle.preds)

all.preds$sse1 <- abs(log(all.preds$time2/ceiling(all.preds$q50)))

ag.pint <- aggregate(data = all.preds, in.predint ~ temp1 + stage + model, mean)
ggplot(data = subset(ag.pint, model != 'MLE')) +
  geom_hline(yintercept = 0.9, size = 0.25) +
  geom_point(aes(x = temp1, y = in.predint, shape = model, 
                 col = model), size = 2.5) +
  scale_colour_manual(values = c('black', 'grey65')) +
  scale_x_continuous(breaks = seq(5, 35, by = 5)) +
  facet_wrap(vars(stage)) +
  theme_minimal() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size = 14),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), 
                                    size = 14),
        strip.text = element_text(size = 12, face = 'bold',
                                  margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10), 
        legend.position = c(5/6, 1/4),
        #panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing = unit(2, "lines")) +
  labs(x = 'Rearing Temperature', y = '% of Observations in Prediction Interval',
       shape = 'Model', colour = 'Model')


sse2.1 <- log(all.preds$time2/all.preds$q50)
sse2.2 <- log(all.preds$time2d/all.preds$q50)
all.preds$sse2 <-  mapply(function(s1, s2) {
  w <- which.min(abs(c(s1, s2)))
  s <- c(s1, s2)[w]
  s.corr <- ifelse(s2 <= 0 & s1 >= 0, 0, s)
  return(abs(s.corr))
}, sse2.1, sse2.2)

ag.sse <- aggregate(data = all.preds, cbind(I(sse1), I(sse2))
                     ~ model + temp1 + stage, mean)

names(ag.sse)[4:5] <- c('Ceiling', 'Best')
m.ag.sse <- melt(data = ag.sse, id.vars = c('model',  'temp1', 'stage'))
#m.ag.sse$value <- exp(m.ag.sse$value)
m.ag.sse$med <- 1
ggplot(data = subset(m.ag.sse, variable == 'Best' & model != 'MLE')) +
  geom_point(aes(x = temp1, y = value,
                 shape = model, col = model), size = 2.5) +
  theme_minimal() + 
  facet_wrap(vars(stage), scales = 'fixed') +
  scale_colour_manual(values = c('black', 'grey65')) +
  scale_x_continuous(breaks = seq(5, 35, by = 5)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), 
                                    size = 14),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0), 
                                    size = 14),
        strip.text.x = element_text(size = 12, face = 'bold',
                                  margin = margin(t = 0, r = 0, b = 10, l = 0)),
        strip.text.y = element_text(size = 12, face = 'bold',
                                    margin = margin(t = 0, r = 0, b = 0, l = 10)),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10), 
        #panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing = unit(2, "lines"),
        legend.position = c(5/6, 1/4)) +
  labs(x = 'Group Index', y = 'Mean Absolute log-Scale Error',
       shape = 'Model', col = 'Model')

ag <- aggregate(data = all.preds, exp(cbind(I(sse1), I(sse2)))
                ~ temp1 + stage + model, mean)

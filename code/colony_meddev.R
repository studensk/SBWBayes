library(tidyverse)
source('code/functions.R')
stages <- paste0('L', 2:6)

m.res.lst <- lapply(provs2, function(p) {
  file <- paste0('data/', p, '_model_results.csv')
  dat <- read.csv(file)
  dat$stage <- rep(stages, nrow(dat)/length(stages))
  dat$province <- p
  dat$HL <- -abs(dat$HL)
  names(dat)[names(dat == 'rho')] <- 'rho25'
  dat$province <- p
  return(dat)
})
m.res.df <- bind_rows(m.res.lst)


results.lst <- lapply(m.res.lst, function(dat) {
  st.lst <- lapply(stages, function(s) {
    st.dat <- subset(dat, stage == s)
    gc <- get_curves(st.dat, spline = TRUE)
    ag <- aggregate(data = gc, rate ~ temp, median)
    ag$province <- unique(dat$province)
    ag$stage <- s
    return(ag)
  })
  st.df <- bind_rows(st.lst)
  return(st.df)
})
results <- bind_rows(results.lst)
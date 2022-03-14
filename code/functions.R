library(dclone)
count <- dplyr::count

## Functions for imputing observed dev rates at lethal temps
fn <- function(x, df, col, r.opt) {
  1/(x*df$time1 + r.opt*(df[,col]))
}

mean.wrap <- function(x, df, ...) {
  mean(fn(x, df,...)) - 1
}

med.wrap <- function(x, df, ...) {
  median(fn(x, df, ...)) - 1
}

## Calculate TA from other parameters
ta.fun <- function(HL, HH, TL, TH) {
  c3 <- 0.0001987
  den <- c3*log(-HL/HH) + (HL/TL) - (HH/TH)
  tphi <- (HL - HH)/den
  return(tphi)
}

## Extract dev rate at 20 degrees C, given set of parameters
post20.fun <- function(post) {
  s <- sapply(1:nrow(post), function(r) {
    row <- post[r,]
    ta <- with(row, ta.fun(HL, HH, TL, TH))
    rate <- 1/(with(row, calc_pred3(20, rho25, HA, TL, HL, TH, HH, ta)))
    return(rate)
  })
  return(s)
}

## Ipmute dev rates at lethal temps using estimates of development at 
##  20 degrees C from posterior draws
r.est12.fun <- function(data, post) {
  data <- data[order(data$temp1),]
  post20 <- post20.fun(post)
  l <- lapply(seq(5, 35, by = 5), function(temp) {
    dtemp <- subset(data, temp1 == temp)
    if (temp %in% c(15, 20, 25)) {
      est1 <- 1/dtemp$time2
      est2 <- 1/dtemp$time2d
      est <- 1/(dtemp$time2-0.5)
      return(data.frame(est1, est2, est))
    }
    else {
      p20 <- median(post20)
      time.df <- subset(dtemp, 
                        select = c(time1, time2, time2.orig, nobs, transfer))
      transfer <- time.df$transfer
      time.lst <- lapply(1:nrow(time.df), function(r) {
        row <- time.df[r,1:3]
        if (row$time1 == 0) {
          row$time1 <- row$time2.orig
          row$time2.orig <- 0
        }
        nobs <- time.df[r,]$nobs
        new <- dclone(row, n.clones = nobs)
        return(new)
      })
      times <- bind_rows(time.lst)
      root <- max(uniroot(med.wrap, interval = c(0, 1),
                          df = times, r.opt = p20,
                          extendInt = 'yes', maxiter = 100000)$root, 0)
      mult <- sapply(transfer, function(x) {ifelse(x, p20, root)})
      delta1 <- 1/(dtemp$time1*root + dtemp$time2*mult)
      delta2 <- 1/(dtemp$time1*root + dtemp$time2d*mult)
      delta <- 1/(dtemp$time1*root + (dtemp$time2-0.5)*mult)
      est1 <- delta1*root
      est2 <- delta2*root
      est <- delta*root
      return(data.frame(est1, est2, est))
    }
  })
  r.est <- bind_rows(l)
  data <- cbind(data, r.est)
  return(data)
}

tempvec <- seq(0, 40, by = 0.1)

## Dev time function
calc_pred3 <- function(temp, rho25, HA, TL, HL, TH, HH, TA) {
  c3 <- 0.0001987
  tK = temp + 273
  num = rho25 * tK/TA * exp(HA/c3* (1/TA-1/tK))
  den1 = exp(HL/c3 *(1/TL-1/tK))
  den2 = exp(HH/c3 *(1/TH-1/tK))
  tau = 1/(num/(1 + den1 + den2))
  return(tau)
}

## Quadratic stage structure function
quadratic <- function(stage, phi, psi, y0) {
  a <- y0*(phi + psi - 2)/8
  b <- y0*(psi - phi)/4
  rho <- a*stage^2 + b*stage + y0
  return(rho)
}

## Obtain either parametric curves or bias-reduced splines from parameter vectors
get_curves <- function(data, temp = tempvec, spline = FALSE) {
  lst <- lapply(1:nrow(data), function(r) {
    row <- data[r,]
    ta <- with(row, ta.fun(HL, HH, TL, TH))
    if (spline) {
      temps <- seq(5, 35, by = 5)
      times <- with(row, calc_pred3(temps, rho25, HA, TL, HL, TH, HH, ta))
      ups <- unlist(row[,grep('upsilon.', names(row))])
      sigma <- row$s_upsilon
      u.ups <- pnorm(ups, 0, 1)
      c.ups <- qcauchy(u.ups, 1, sigma)
      times <- times*c.ups
      rates <- 1/times
      if (length(temp) != length(temps)) {
        i.spline <- interpSpline(temps, rates)
        df <- as.data.frame(predict(i.spline, temp))
        names(df) <- c('temp', 'rate')
        df$rate <- pmax(df$rate, 0)
        df$time <- 1/df$rate
      }
      else {df <- data.frame('temp' = temps,
                             'rate' = pmax(rates, 0),
                             'time' = 1/pmax(rates, 0))}
    }
    else {
      times <- with(row, calc_pred3(temp, rho25, HA, TL, HL, TH, HH, ta))
      rates <- 1/times
      df <- data.frame('temp' = temp, 'rate' = rates, 'time' = times)
    }
    df$index <- r
    return(df)
  })
  df2 <- bind_rows(lst)
  return(df2)
}

get_curves_opt <- function(data, temp = tempvec, spline = FALSE) {
  lst <- lapply(1:nrow(data), function(r) {
    row <- data[r,]
    ta <- with(row, ta.fun(HL, HH, TL, TH))
    if (spline) {
      temps <- seq(5, 35, by = 5)
      times <- with(row, calc_pred3(temps, rho25, HA, TL, HL, TH, HH, ta))
      ups <- unlist(row[,grep('upsilon.', names(row))])
      sigma <- row$s_upsilon
      u.ups <- pnorm(ups, 0, 1)
      c.ups <- qcauchy(u.ups, 1, sigma)
      times <- times*c.ups
      rates <- 1/times
      if (length(temp) != length(temps)) {
        i.spline <- interpSpline(temps, rates)
        df <- as.data.frame(predict(i.spline, temp))
        names(df) <- c('temp', 'rate')
        df$rate <- pmax(df$rate, 0)
        df$time <- 1/df$rate
      }
      else {df <- data.frame('temp' = temps,
                             'rate' = pmax(rates, 0),
                             'time' = 1/pmax(rates, 0))}
    }
    else {
      times <- with(row, calc_pred3(temp, rho25, HA, TL, HL, TH, HH, ta))
      rates <- 1/times
      w <- which.max(rates)
      df <- data.frame('temp' = temp, 'rate' = rates, 'time' = times,
                       'r.opt' = rates[w], 't.opt' = temp[w])
    }
    df$index <- r
    return(df)
  })
  df2 <- bind_rows(lst)
  return(df2)
}

prov.factor <- function(data, name = 'province') {
  data[,name] <- factor(data[,name], 
                        levels = c('IN', 'AB', 'QC', 'ON', 'NB2', 'IPU'),
                        labels = c('Northwest Territories', 'Alberta',
                                   'Quebec', 'Ontario', 'New Brunswick', 'Lab-Reared'))
  return(data)
}




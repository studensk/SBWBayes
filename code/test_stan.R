library(DPQ)

## Test stan code for stage-structured model
stan.test.structured <- function(data, params) {
  s <- data$stage + 1
  params$alpha <- rep(0.05, length(unique(data$stage)))
  params$upsilon <- rep(0, length(unique(data$block)))
  
  rho25 <- with(params, quadratic(data$stage - 2, phi_rho, psi_rho, y0_rho))
  rho25.new <- with(params, rho25[s]*(exp(alpha[s]*s_alpha - (s_alpha^2)/2)))
  
  params$HL <- -abs(params$HL)
  
  TA <- with(params, ta.fun(HL[s], HH[s], TL[s], TH[s]))
  
  tpred1 <- with(params, calc_pred3(data$temp1, rho25.new, HA[s], TL[s],
                                    HL[s], TH[s], HH[s], TA[s]))
  tpred2 <- with(params, calc_pred3(data$temp2, rho25.new, HA[s], TL[s],
                                    HL[s], TH[s], HH[s], TA[s]))
  
  upsilon <- params$upsilon[data$block]
  # u.upsilon <- pnorm(upsilon, 0, 1)
  # c.upsilon <- qcauchy(u.upsilon, 1, params$s_upsilon)
  c.upsilon <- exp(upsilon)
  
  tpred1 <- tpred1*c.upsilon
  tpred2 <- tpred2*c.upsilon

  tpred1 <- pmin(tpred1, 30)
  tpred2 <- pmin(tpred2, 30)
  
  epsm1 <- log(data$time1/tpred1 + data$time2d/tpred2)
  epsij <- log(data$time1/tpred1 + data$time2/tpred2)
  
  epsm1_std <- with(params, epsm1/s_eps[s])
  epsij_std <- with(params, epsij/s_eps[s])
  
  pnorm_ij <- pnorm(epsij_std, 0, 1, log.p = TRUE)
  pnorm_m1 <- pnorm(epsm1_std, 0, 1, log.p = TRUE)
  
  diff <- sapply(1:nrow(data), function(i) {
    x <- data$time2d[i]
    ifelse(x == 0, pnorm_ij[i], logspace.sub(pnorm_ij[i], pnorm_m1[i]))
  })
  
  data <- as.data.frame(data)
  df <- data.frame(data, tpred1, tpred2, epsm1, epsij, epsm1_std, 
                   epsij_std, pnorm_ij, pnorm_m1, diff)
  ll <- sum(df$nobs*df$diff)
  
  return(data.frame(data, tpred1, tpred2, epsm1, epsij, epsm1_std, 
               epsij_std, pnorm_ij, pnorm_m1, diff))
}

## Test stan code for model with no stage structure 
stan.test <- function(data, params) {
  params$HL <- -params$HL
  tpred1 <- with(params,calc_pred3(data$temp1, rho25, HA, TL, HL, TH, HH, 
                                   ta.fun(HL, HH, TL, TH)))
  tpred2 <- with(params, calc_pred3(data$temp2, rho25, HA, TL, HL, TH, HH,
                                    ta.fun(HL, HH, TL, TH)))
  upsilon <- params$upsilon[data$block + 1]
  u.upsilon <- pnorm(upsilon, 0, 1)
  c.upsilon <- qcauchy(u.upsilon, 1, params$s_upsilon)
  
  tpred1 <- with(params, tpred1*c.upsilon)
  tpred2 <- with(params, tpred2*c.upsilon)
  
  epsm1 <- log(data$time1/tpred1 + data$time2d/tpred2)
  epsij <- log(data$time1/tpred1 + data$time2/tpred2)
  
  epsm1_std <- with(params, epsm1/s_eps)
  epsij_std <- with(params, epsij/s_eps)
  
  pnorm_ij <- pnorm(epsij_std, 0, 1, log.p = TRUE)
  pnorm_m1 <- pnorm(epsm1_std, 0, 1, log.p = TRUE)
  
  diff <- log(exp(pnorm_ij) - exp(pnorm_m1))
  
  data <- as.data.frame(data)
  df <- data.frame(data, tpred1, tpred2, epsm1, epsij, epsm1_std, 
                   epsij_std, pnorm_ij, pnorm_m1, diff)
  ll <- sum(df$nobs*df$diff)
  
  return(df)
}

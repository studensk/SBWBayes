## Get data from kstudens/SBWPhenoData on github
sbw2 <- SBWPhenoData::DevTimeSBWv2
stages <- c('L2', 'L3', 'L4', 'L5', 'L6')
provs <- unique(sapply(names(sbw2), function(x) {strsplit(x, split = '_')[[1]][1]}))

## Conform data to correct input format
all.days <- lapply(1:length(sbw2), function(x) {
  dat <- sbw2[[x]]
  index <- strsplit(names(sbw2)[[x]], split = '_')[[1]]
  prov <- index[1]
  temp <- as.numeric(index[2])
  gen <- index[3]
  obs <- dat[,2:(ncol(dat)-1)]
  l <- as.list(obs)
  if (length(l) < 10) {
    time2 <- unlist(l)
    time1 <- rep(0, length(time2))
    stage <- rep(stages, each = length(l[[1]]))
    df <- data.frame(stage, time2, time1)
    df$temp2 <- temp
    df$transfer <- FALSE
  }
  else {
    odds <- seq(1, length(l)-1, by = 2)
    evens <- seq(2, length(l), by = 2)
    time1 <- unlist(l[odds])
    time2 <- unlist(l[evens])
    stage <- rep(stages, each = length(l[[1]]))
    df <- data.frame(stage, time2, time1)
    df$temp2 <- 20
    df$transfer <- TRUE
  }
  df$temp1 <- temp
  df$province <- prov
  df$generation <- gen
  df$block <- x
  return(df)
})
all.days.df <- bind_rows(all.days)
all.days.df$time2 <- sapply(all.days.df$time2, function(y) {
  if (y < 0.5) {
    y <- ifelse(y < 0.25, 0, 0.5)
  }
  return(y)
})
(all.days.df <- all.days.df
  %>% group_by(across())
  %>% count(name = 'nobs') 
  %>% mutate(time2d = pmax(0, time2 - 1)))

## Estimated dev times at 20 degrees C
r20 <- lapply(provs, function(p) {
  s <- sapply(stages, function(stg) {
    sub <- subset(all.days.df, province == p & temp1 == 20 & stage == stg)
    t.est <- median(unlist(lapply(1:nrow(sub), function(r) {
      rp <- rep(sub$time2[r], sub$nobs[r])
      return(rp)
    })))
    r.est <- 1/(t.est - 0.5)
    r1 <- 1/(t.est - 1)
    r2 <- 1/t.est
    r.est2 <- mean(c(r1, r2))
    r <- mean(c(r.est, r.est2))
    return(r)
  })
  names(s) <- stages 
  return(s)
})
names(r20) <- provs

##########

all.days.df <- all.days.df %>% 
  mutate(block = paste(generation, province, temp1, stage, sep = "_"))

all.days.df$block <- as.numeric(factor(all.days.df$block))

## Manage data in which development at lethal temps finishes early
early <- which(all.days.df$time2 == 0)
early.sub <- all.days.df[early,]
early.sub$time2 <- early.sub$time1
early.sub$time1 <- 0
early.sub$time2d <- pmax(early.sub$time2 - 1, 0)
early.sub$temp2 <- early.sub$temp1
early.sub$transfer <- FALSE

all.days.df[early,] <- early.sub

all.days.df <- subset(all.days.df,
                      (province == 'IPU' & generation == 'F0') |
                        (province != 'IPU' & generation == 'F1'))

## Implement double interval censoring for stages > L2
all.days.df$time2.orig <- all.days.df$time2
all.days.df$time2 <- all.days.df$time2.orig + 
  sapply(all.days.df$stage, function(x) {ifelse(x == 'L2', 0, 1)})

## Impute observed development rates
all.days.df$r.est2 <- sapply(1:nrow(all.days.df), function(r) {
  if (!all.days.df$transfer[r]) {return(1/(all.days.df$time2d[r]))}
  else {
    temp <- all.days.df$temp1[r]
    blck <- all.days.df$block[r]
    prov <- all.days.df$province[r]
    stg <- all.days.df$stage[r]
    time.df <- subset(all.days.df, temp1 == temp & block == blck & stage == stg, 
                      select = c(time1, time2, time2.orig, time2d, nobs, transfer))
    time.lst <- lapply(1:nrow(time.df), function(r2) {
      row <- time.df[r2,1:(length(time.df)-2)]
      nobs <- time.df[r2,]$nobs
      transfer <- time.df[r2,]$transfer
      if (!transfer) {
        t1 <- row$time1
        t2 <- row$time2d
        row$time1 <- t2
        row$time2d <- t1
      }
      new <- dclone(row, n.clones = nobs)
      return(new)
    })
    times <- bind_rows(time.lst)
    root <- max(uniroot(med.wrap, interval = c(0, 1),
                        df = as.data.frame(time.df), r.opt = r20[[prov]][stg],
                        col = 'time2d',
                        extendInt = 'yes', maxiter = 100000)$root, 0)
    est <- (root*all.days.df$time1[r] + 
              r20[[prov]][stg]*(all.days.df$time2d[r]))/root
    return(1/est)
  }
})

## Impute observed development rates
all.days.df$r.est1 <- sapply(1:nrow(all.days.df), function(r) {
  if (!all.days.df$transfer[r]) {return(1/(all.days.df$time2[r]))}
  else {
    temp <- all.days.df$temp1[r]
    blck <- all.days.df$block[r]
    prov <- all.days.df$province[r]
    stg <- all.days.df$stage[r]
    time.df <- subset(all.days.df, temp1 == temp & block == blck & stage == stg, 
                      select = c(time1, time2, time2.orig, nobs, transfer))
    time.lst <- lapply(1:nrow(time.df), function(r2) {
      row <- time.df[r2,1:3]
      nobs <- time.df[r2,]$nobs
      transfer <- time.df[r2,]$transfer
      if (!transfer) {
        t1 <- row$time1
        t2 <- row$time2.orig
        row$time1 <- t2
        row$time2 <- t1
      }
      new <- dclone(row, n.clones = nobs)
      return(new)
    })
    times <- bind_rows(time.lst)
    root <- max(uniroot(med.wrap, interval = c(0, 1),
                        df = as.data.frame(time.df), r.opt = r20[[prov]][stg],
                        col = 'time2.orig',
                        extendInt = 'yes', maxiter = 100000)$root, 0)
    est <- (root*all.days.df$time1[r] + 
              r20[[prov]][stg]*(all.days.df$time2[r]))/root
    return(1/est)
  }
})

#write.csv(all.days.df, 'data/all_days_df.csv')

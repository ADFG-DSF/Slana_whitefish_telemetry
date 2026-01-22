library(jagshelper)
library(jagsUI)

library(tidyverse)

tracker <- read_csv("FDS_2025/flat_data/tracker.csv") %>%
  filter(!is.na(Code)) %>%
  select("10/16/2024", "10/24/2024", "11/1/2024",  "12/27/2024", "3/17/2025",  "5/19/2025",
         "6/11/2025",  "7/22/2025",  "9/8/2025",   "9/24/2025",  "10/10/2025", "10/22/2025",
         "10/29/2025")

# # data check
# for(j in 1:ncol(tracker)) print(table(tracker[,j], useNA='ifany'))


Y <- matrix(nrow=nrow(tracker), ncol=ncol(tracker))
for(j in 1:ncol(tracker)) {
  Y[,j] <- ifelse(tracker[,j]=="A", 1,
                      ifelse(tracker[,j]=="M", 0, NA))
}
Yraw <- Y

for(i in 1:nrow(Y)) {
  for(j in 1:ncol(Y)) {
    if(!is.na(Y[i,j])) {
      if(Y[i,j]==0) {
        Y[i,(j:ncol(Y))] <- 0  # dead fish stay dead
      }
      if(Y[i,j]==1) {
        Y[i, 1:j] <- 1  # alive fish were alive before
      }
    }
  }
}


# specify model, which is written to a temporary file
hm_jags <- tempfile()
cat('model {
  for(i in 1:ni) {
    Y[i,1] ~ dbinom(p[1], 1)
    # Y[i,1] ~ dbinom(pcum[1], 1)
    for(j in 2:nj) {
      Y[i,j] ~ dbinom(p[j], Y[i, j-1])
      # Y[i,j] ~ dbinom(pcum[j], 1)
    }
  }

  for(j in 1:nj) {
    p[j] ~ dbeta(0.5, 0.5)
    # pcum[j] ~ dbeta(0.5, 0.5)
  }

}', file=hm_jags)


# bundle data to pass into JAGS
hm_data <- list(Y = Y,
                ni = nrow(Y),
                nj = ncol(Y))

# JAGS controls
niter <- 10000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  hm_jags_out <- jagsUI::jags(model.file=hm_jags, data=hm_data,
                              parameters.to.save=c("Y", "p","pcum"),
                              n.chains=ncores, parallel=T, n.iter=niter,
                              n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

nbyname(hm_jags_out)
plotRhats(hm_jags_out)
traceworstRhat(hm_jags_out, parmfrow = c(3, 3))


dates <- c("10/16/2024", "10/24/2024", "11/1/2024",  "12/27/2024", "3/17/2025",  "5/19/2025",
           "6/11/2025",  "7/22/2025",  "9/8/2025",   "9/24/2025",  "10/10/2025", "10/22/2025",
           "10/29/2025")

caterpillar(hm_jags_out, "p", xax=dates)
parmar <- par("mar")
par(mar=parmar + c(1,0,0,0))
caterpillar(hm_jags_out, "p", xax=dates, las=2, main="Survival probability")
par(mar=parmar)

thetab <- apply(Y, 2, table, useNA='ifany')
colnames(thetab) <- dates
thetab

modsurv <- hm_jags_out$mean$Y
par(mar=parmar + c(1,0,0,0))
plot(NA, xlim=c(1, ncol(Y)), ylim=c(0, 1), xaxt='n', ylab="Modeled survival", xlab="")
axis(side=1, at=seq_along(dates), labels=dates, las=2)
for(i in 1:nrow(Y)) {
  lines(jitter(modsurv[i,], 0.1), col=adjustcolor(1, alpha.f=.4))
}
points(colMeans(modsurv), pch=16)
lines(colMeans(modsurv), lty=2, lwd=3)
par(mar=parmar)

cumulsurv <- NA*hm_jags_out$sims.list$p
for(j in 1:ncol(cumulsurv)) {
  for(i in 1:nrow(cumulsurv)) { # this is a kludge but I don't care
    cumulsurv[i,j] <- prod(hm_jags_out$sims.list$p[i, 1:j])
  }
}
par(mar=parmar + c(1,0,0,0))
caterpillar(cumulsurv, xax=dates, las=2, main="Cumulative survival probability")
# points(thetab[2,]/colSums(thetab[1:2,]))
par(mar=parmar)





plotcor_jags(hm_jags_out, p="p")


# and what if we don't go hidden markov-y and just estimate from each pair of surveys
Y1 <- cbind(1,Y)
Ylist <- list()
for(i in 1:ncol(Y)) {
  Ylist[[i]] <- table(factor(Y1[,i], levels=0:1),
                      factor(Y1[,i+1], levels=0:1),
                      useNA='always')
}
propest <- sapply(Ylist, \(x) x[2,2]/sum(x[2,1:2]))
caterpillar(hm_jags_out, "p", xax=dates)
points(propest, pch="x")
bp <- 0.5  # beta prior hyperparams
betalo <- sapply(Ylist, \(x) qbeta(0.025, x[2,2]+bp, x[2,1]+bp))
betahi <- sapply(Ylist, \(x) qbeta(0.975, x[2,2]+bp, x[2,1]+bp))
betamid <- sapply(Ylist, \(x) qbeta(0.5, x[2,2]+bp, x[2,1]+bp))
points(betalo, pch="-")
points(betahi, pch="-")
points(betamid)





### trying the HM model again, but including fish-level predictors

source("FDS_2025/R/1_Slana_telemetry_fromobjectives.R")


tagging_data <- read_csv("FDS_2025/flat_data/tracker.csv") %>%
  filter(!is.na(Code)) %>%
  select(Length, Sex)

# specify model, which is written to a temporary file
hm_jags_byindiv <- tempfile()
cat('model {
  for(i in 1:ni) {
    Y[i,1] ~ dbin(p[i,1], 1)
    for(j in 2:whichlast[i]) {#
      Y[i,j] ~ dbin(p[i,j], Y[i, j-1])
    }
  }

  for(i in 1:ni) {
    for(j in 1:whichlast[i]) {
      # logit(p[i,j]) <- b0[j] + bLength#*Length[i]   # WHY IS THIS SUCH A PROBLEM???

      # p[i,j] <- exp(b0[j] + bLength*Length[i]) / (1 + exp(b0[j] + bLength*Length[i]))

      # p[i,j] <- exp(b0[j] + bLength*Length[i] + bDist*dist_p_survey[i] + bLake*p_lake[i]) /
      #      (1 + exp(b0[j] + bLength*Length[i] + bDist*dist_p_survey[i] + bLake*p_lake[i]))

# ---------------

      # p[i,j] <- exp(b0[j] + bLake*lake[i,j] + bLength*Length[i] + bDist*dist_p_survey[i]) /
      #      (1 + exp(b0[j] + bLake*lake[i,j] + bLength*Length[i] + bDist*dist_p_survey[i]))

      # p[i,j] <- exp(b0[j] + bLake*lake[i,j]) /
      #      (1 + exp(b0[j] + bLake*lake[i,j]))

      # p[i,j] <- exp(b0[j]) /
      #      (1 + exp(b0[j]))

      # -----------

      p[i,j] <- exp(b0[j, lake[i,j]+1]) /
           (1 + exp(b0[j, lake[i,j]+1]))
    }
  }

  for(j in 1:nj) {
    # b0[j] ~ dnorm(0, 1)
    b0[j,1] ~ dnorm(0, 1)
    b0[j,2] ~ dnorm(0, 1)
  }
  bLength ~ dnorm(0, 1)
  bDist ~ dnorm(0, 1)
  bLake ~ dnorm(0, 1)

}', file=hm_jags_byindiv)



# last known location?
lake_mat_df <- ptdata %>%
  select(Code, Survey, lake) %>%
  mutate(Code = factor(Code)) %>%
  filter(Survey >= 1) %>%
  pivot_wider(names_from = Survey, values_from = lake, id_expand = T) %>% data.frame
lake_mat <- as.matrix(lake_mat_df[,-1])*1
for(i in 1:nrow(lake_mat)) {  # this is a terrible kludge
  if(is.na(lake_mat[i,1])) {
    lake_mat[i,1] <- floor(median(lake_mat[i,], na.rm=TRUE))
  }
  for(j in 2:ncol(lake_mat)) {
    if(is.na(lake_mat[i,j])) {
      lake_mat[i,j] <- lake_mat[i,j-1]
    }
  }
}


# bundle data to pass into JAGS
# whichrows <- rowSums(!is.na(Y)) ==13
whichrows <- rowSums(!is.na(Y)) > 1 & !is.na(by_indiv$homerange)

Y1 <- Y[whichrows,]
cc <- function(x) x - mean(x, na.rm=TRUE)
hm_data <- list(Y = Y1,
                Length = cc(tagging_data$Length)[whichrows],
                homerange = cc(by_indiv$homerange)[whichrows],
                totaldist = cc(by_indiv$totaldist)[whichrows],
                dist_p_survey = cc(by_indiv$totaldist/(by_indiv$n_surveys-1))[whichrows],
                p_lake = cc(by_indiv$inside_lake/(by_indiv$inside_lake+by_indiv$outside_lake))[whichrows],
                lake = lake_mat[whichrows,],
                whichlast = apply(Y1, 1, \(x) max(which(!is.na(x)))),
                ni = nrow(Y1),
                nj = ncol(Y1))

# JAGS controls
niter <- 10000
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  hm_jags_byindiv_out <- jagsUI::jags(model.file=hm_jags_byindiv, data=hm_data,
                              parameters.to.save=c("Y", "b0", "p",
                                                   "bLength", "bDist", "bLake"),
                              n.chains=ncores, parallel=T, n.iter=niter,
                              n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

nbyname(hm_jags_byindiv_out)
plotRhats(hm_jags_byindiv_out)
traceworstRhat(hm_jags_byindiv_out, parmfrow = c(3, 3))

# caterpillar(hm_jags_byindiv_out, "b0")
par(mfrow=c(2,2))
# caterpillar(hm_jags_byindiv_out, "b0", transform = "expit")
caterpillar(hm_jags_byindiv_out, "b0", transform = "expit", column=1, xax=dates, las=2)
caterpillar(hm_jags_byindiv_out, "b0", transform = "expit", column=2, add=T, x=1:13+.2, col=2)
plotdens(hm_jags_byindiv_out, "bLength")
abline(v=0, lty=2)
plotdens(hm_jags_byindiv_out, "bDist")
abline(v=0, lty=2)
plotdens(hm_jags_byindiv_out, "bLake")
abline(v=0, lty=2)



# specify model, which is written to a temporary file
hm_jags_detect <- tempfile()
cat('model {
  for(i in 1:ni) {
    Y[i,1] ~ dbin(phi[lake[i,1]+1, 1], 1)
    X[i,1] ~ dbin(p[lake[i,1]+1, Y[i,1]+1], 1)
    for(j in 2:whichlast[i]) {#
      Y[i,j] ~ dbin(phi[lake[i,j]+1, j], Y[i, j-1])
      X[i,j] ~ dbin(p[lake[i,j]+1, Y[i,j]+1], 1)
    }
  }

  for(ip in 1:2) {
    for(jp in 1:2) {
      p[ip,jp] ~ dbeta(1,1)
    }
  }

  for(iphi in 1:2) {
    for(jphi in 1:nj) {
      phi[iphi,jphi] ~ dbeta(1,1)
    }
  }

}', file=hm_jags_detect)

whichrows <- rowSums(!is.na(Y)) > 1 & !is.na(by_indiv$homerange)

Y1 <- Y[whichrows,]
cc <- function(x) x - mean(x, na.rm=TRUE)
hm_data <- list(Y = Y1,
                X = 1*(!is.na(Yraw[whichrows,])),
                lake = lake_mat[whichrows,],
                whichlast = apply(Y1, 1, \(x) max(which(!is.na(x)))),
                ni = nrow(Y1),
                nj = ncol(Y1))

# JAGS controls
niter <- 10000
# ncores <- 3
ncores <- 8# min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  hm_jags_detect_out <- jagsUI::jags(model.file=hm_jags_detect, data=hm_data,
                                      parameters.to.save=c("Y", "p", "phi"),
                                      n.chains=ncores, parallel=T, n.iter=niter,
                                      n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

nbyname(hm_jags_detect_out)
plotRhats(hm_jags_detect_out)
traceworstRhat(hm_jags_detect_out, parmfrow = c(2,2))

par(mfrow=c(1,1))
caterpillar(hm_jags_detect_out, "phi", row=1, ylim=0:1)
caterpillar(hm_jags_detect_out, "phi", row=2, x=1:13+.2, col=2, add=T)
caterpillar(hm_jags_detect_out, "p", row=1, ylim=0:1)
caterpillar(hm_jags_detect_out, "p", row=2, x=1:2+.2, col=2, add=T)

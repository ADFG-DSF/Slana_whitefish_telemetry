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
      p[ip,jp] ~ dbeta(0.5, 0.5)
    }
  }

  for(iphi in 1:2) {
    for(jphi in 1:nj) {
      phi[iphi,jphi] ~ dbeta(0.5, 0.5)
    }
  }
  pvec[1] <- p[1,1]
  pvec[2] <- p[1,2]
  pvec[3] <- p[2,1]
  pvec[4] <- p[2,2]

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
niter <- 50000
# ncores <- 3
ncores <- 8# min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  hm_jags_detect_out <- jagsUI::jags(model.file=hm_jags_detect, data=hm_data,
                                      parameters.to.save=c("Y", "p", "phi","pvec"),
                                      n.chains=ncores, parallel=T, n.iter=niter,
                                      n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

nbyname(hm_jags_detect_out)
plotRhats(hm_jags_detect_out)
traceworstRhat(hm_jags_detect_out, parmfrow = c(2,2))

par(mfrow=c(1,1))
par(mar=c(6,4,4,1))
caterpillar(hm_jags_detect_out, "phi", row=1,
            ylim=0:1, xax=dates, las=2,
            main="Survival Probability")
caterpillar(hm_jags_detect_out, "phi", row=2, x=1:13+.2, col=2, add=T)
# lines(y=hm_jags_detect_out$q50$phi[1,], col=adjustcolor(4, alpha.f=.5), x=1:13, lwd=2)
# lines(y=hm_jags_detect_out$q50$phi[2,], col=adjustcolor(2, alpha.f=.5), x=1:13+.2, lwd=2)
legend("bottomleft", legend=c("River", "Lake"), col=c(4,2), lwd=3)
par(mar=parmar)


# checking to see if the data support inferences
str(hm_data)

# Y is survival
# X is detection
# lake is lake
lakemat <- rivermat <- matrix(nrow=5, ncol=13)
rownames(lakemat) <- rownames(rivermat) <- c("Alive_Detected", "Alive_NotDetected","Died_Detected", "Died_Prev","MIA")
lakemat[1,1] <- sum((hm_data$lake[,1]==1) & (hm_data$Y[,1]==1) & (hm_data$X[,1]==1), na.rm=TRUE)
lakemat[2,1] <- sum((hm_data$lake[,1]==1) & (hm_data$Y[,1]==1) & (hm_data$X[,1]==0), na.rm=TRUE)
lakemat[3,1] <- sum((hm_data$lake[,1]==1) & (hm_data$Y[,1]==0) & (hm_data$X[,1]==1), na.rm=TRUE)
lakemat[4,1] <- 0
lakemat[5,1] <- sum((hm_data$lake[,1]==1) & is.na(hm_data$Y[,1]) & (hm_data$X[,1]==0), na.rm=TRUE)
for(j in 2:13) {
  lakemat[1,j] <- sum((hm_data$lake[,j]==1) & (hm_data$Y[,j]==1) & (hm_data$X[,j]==1), na.rm=TRUE)
  lakemat[2,j] <- sum((hm_data$lake[,j]==1) & (hm_data$Y[,j]==1) & (hm_data$X[,j]==0), na.rm=TRUE)
  lakemat[3,j] <- sum((hm_data$lake[,j]==1) & (hm_data$Y[,j]==0) & (hm_data$Y[,j-1]==1) & (hm_data$X[,j]==1), na.rm=TRUE)
  lakemat[4,j] <- sum((hm_data$lake[,j]==1) & (hm_data$Y[,j-1]==0), na.rm=TRUE)
  lakemat[5,j] <- sum((hm_data$lake[,j]==1) & is.na(hm_data$Y[,j]) & (hm_data$X[,j]==0), na.rm=TRUE)
}
rivermat[1,1] <- sum((hm_data$lake[,1]==0) & (hm_data$Y[,1]==1) & (hm_data$X[,1]==1), na.rm=TRUE)
rivermat[2,1] <- sum((hm_data$lake[,1]==0) & (hm_data$Y[,1]==1) & (hm_data$X[,1]==0), na.rm=TRUE)
rivermat[3,1] <- sum((hm_data$lake[,1]==0) & (hm_data$Y[,1]==0) & (hm_data$X[,1]==1), na.rm=TRUE)
rivermat[4,1] <- 0
rivermat[5,1] <- sum((hm_data$lake[,1]==0) & is.na(hm_data$Y[,1]) & (hm_data$X[,1]==0), na.rm=TRUE)
for(j in 2:13) {
  rivermat[1,j] <- sum((hm_data$lake[,j]==0) & (hm_data$Y[,j]==1) & (hm_data$X[,j]==1), na.rm=TRUE)
  rivermat[2,j] <- sum((hm_data$lake[,j]==0) & (hm_data$Y[,j]==1) & (hm_data$X[,j]==0), na.rm=TRUE)
  rivermat[3,j] <- sum((hm_data$lake[,j]==0) & (hm_data$Y[,j]==0) & (hm_data$Y[,j-1]==1) & (hm_data$X[,j]==1), na.rm=TRUE)
  rivermat[4,j] <- sum((hm_data$lake[,j]==0) & (hm_data$Y[,j-1]==0), na.rm=TRUE)
  rivermat[5,j] <- sum((hm_data$lake[,j]==0) & is.na(hm_data$Y[,j]) & (hm_data$X[,j]==0), na.rm=TRUE)
}
lakemat
rivermat
colSums(lakemat+rivermat)

xlake <- colSums(lakemat[1:2,])
nlake <- colSums(lakemat[1:3,])
plake <- xlake/nlake
se_plake <- sqrt(plake*(1-plake)/(nlake-1))

xriver <- colSums(rivermat[1:2,])
nriver <- colSums(rivermat[1:3,])
priver <- xriver/nriver
se_priver <- sqrt(priver*(1-priver)/(nriver-1))



# make that last plot again, for the purpose of adding check points
par(mfrow=c(1,1))
par(mar=c(6,4,4,1))
caterpillar(hm_jags_detect_out, "phi", row=1,
            ylim=0:1, xax=dates, las=2,
            main="Survival Probability")
caterpillar(hm_jags_detect_out, "phi", row=2, x=1:13+.2, col=2, add=T)
# lines(y=hm_jags_detect_out$q50$phi[1,], col=adjustcolor(4, alpha.f=.5), x=1:13, lwd=2)
# lines(y=hm_jags_detect_out$q50$phi[2,], col=adjustcolor(2, alpha.f=.5), x=1:13+.2, lwd=2)
grid(nx=NA, ny=NULL)
legend("bottomleft", legend=c("River", "Lake"), col=c(4,2), lwd=3)
par(mar=parmar)

points(priver, col=4)
points(priver-2*se_priver, pch="-", col=4)
points(priver+2*se_priver, pch="-", col=4)
points(y=plake, x=1:13+.2, col=2)
points(y=plake-2*se_plake, pch="-", x=1:13+.2, col=2)
points(y=plake+2*se_plake, pch="-", x=1:13+.2, col=2)


# row is lake, column is survival
# caterpillar(hm_jags_detect_out, "p", row=1, ylim=0:1, x=1:2-.1,
#             xax=c("Dead","Alive"), main="Detection Probability")
# caterpillar(hm_jags_detect_out, "p", row=2, x=1:2+.1, col=2, add=T)
# legend("bottomleft", legend=c("River", "Lake"), col=c(4,2), lwd=3)
#
# caterpillar(hm_jags_detect_out, "p", column=1, ylim=0:1, x=1:2-.1,
#             xax=c("River", "Lake"), main="Detection Probability")
# caterpillar(hm_jags_detect_out, "p", column=2, x=1:2+.1, col=2, add=T)
# legend("bottomleft", legend=c("Dead","Alive"), col=c(4,2), lwd=3)

caterpillar(hm_jags_detect_out, "pvec", main="Detection Probability", col=c(4,4,2,2),
            xax=rep(c("Dead","Alive"), 2),
            xlab="River                                             Lake")
grid(nx=NA, ny=NULL)


## compiling a results table to write to csv
part1 <- data.frame(parameter="Survival Probability",
                    location="River",
                    date=dates,
                    Survival="",
                    Estimate=hm_jags_detect_out$q50$phi[1,],
                    SE=hm_jags_detect_out$sd$phi[1,],
                    CI95_lo=hm_jags_detect_out$q2.5$phi[1,],
                    CI95_hi=hm_jags_detect_out$q97.5$phi[1,])
part2 <- data.frame(parameter="Survival Probability",
                    location="Lake",
                    date=dates,
                    Survival="",
                    Estimate=hm_jags_detect_out$q50$phi[2,],
                    SE=hm_jags_detect_out$sd$phi[2,],
                    CI95_lo=hm_jags_detect_out$q2.5$phi[2,],
                    CI95_hi=hm_jags_detect_out$q97.5$phi[2,])

part3 <- data.frame(parameter="Detection Probability",
                    location=c(rep("River", 2), rep("Lake", 2)),
                    date="",
                    Survival=rep(c("Dead","Alive"), 2),
                    Estimate=hm_jags_detect_out$q50$pvec,
                    SE=hm_jags_detect_out$sd$pvec,
                    CI95_lo=hm_jags_detect_out$q2.5$pvec,
                    CI95_hi=hm_jags_detect_out$q97.5$pvec)

# write.csv(rbind(part1, part2, part3), file="FDS_2025/output_tables/HM_survival_detection.csv")

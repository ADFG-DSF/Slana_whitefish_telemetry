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
    for(j in 2:nj) {
      Y[i,j] ~ dbinom(p[j], Y[i, j-1])
    }
  }

  for(j in 1:nj) {
    p[j] ~ dbeta(1, 1)
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
                              parameters.to.save=c("Y", "p"),
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
caterpillar(hm_jags_out, "p", xax=dates, las=2)

thetab <- apply(Y, 2, table, useNA='ifany')
colnames(thetab) <- dates
thetab

modsurv <- hm_jags_out$mean$Y
plot(NA, xlim=c(1, ncol(Y)), ylim=c(0, 1))
for(i in 1:nrow(Y)) {
  lines(jitter(modsurv[i,], 0.1), col=adjustcolor(1, alpha.f=.4))
}

cumulsurv <- NA*hm_jags_out$sims.list$p
for(j in 1:ncol(cumulsurv)) {
  for(i in 1:nrow(cumulsurv)) { # this is a kludge but I don't care
    cumulsurv[i,j] <- prod(hm_jags_out$sims.list$p[i, 1:j])
  }
}
caterpillar(cumulsurv)
points(thetab[2,]/colSums(thetab[1:2,]))

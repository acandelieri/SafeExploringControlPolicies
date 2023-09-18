rm(list=ls()); graphics.off(); cat("\014")
source("DataGenerator.R")
library(DiceKriging)
library(plot3D)

dgt <- 4

max.inflow <- 10
eta <- 0.85

H.min <- 5
H.max <- 50
h0.mean <- 25
h0.sd <- 10

n.days <- 365


# simulation step 1 hour - horizon (i.e., T) 24 hours
mu.xi <- c(2.2,2,1.8,1.8,2.3,4,7.3,7.8,7.4,6.3,5.5,4.7,4.5,4.1,4.1,4.5,5.7,6,5.3,4.5,4.0,3.5,2.7,2.2) # water demand (mean at each step)
mu.xi <- mu.xi + 5
#sd.xi <- rep(2,24) # water demand (sd at each step)
sd.xi <- rep(0.10*diff(range(mu.xi)),24) # water demand (sd at each step)
policy <- safe.policy # safe policy
prices <- c(rep(1,6),rep(3,2),rep(2,10),rep(3,2),rep(2,2),rep(1,2)) # energy price tariff


cat("> Starting generation of data for",n.days,"days...\n")
knowledge <- tank.episodes.generator( N=n.days, # number of days (i.e., 1 year)
                                      H=c(H.min,H.max), # min-max tank level
                                      X=c(0,max.inflow), # min-max pump
                                      mu.h=h0.mean, # initial tank level normal distribution (mean)
                                      sd.h=h0.sd, # initial tank level normal distribution (sd)
                                      mu.xi=mu.xi, # water demand (mean at each step)
                                      sd.xi=sd.xi, # water demand (sd at each step)
                                      policy=policy, # safe policy
                                      prices=prices, # energy price tariff
                                      eta=eta, # pump efficiency coefficient
                                      seed=42 #  seed
                                      )

knowledge$s <- round(knowledge$s,dgt)
knowledge$a <- round(knowledge$a,dgt)
knowledge$s_ <- round(knowledge$s_,dgt)
knowledge$r <- round(knowledge$r,dgt)

stopifnot( length(which(knowledge$s<H.min))==0 | length(which(knowledge$s_<H.min))==0)
stopifnot( length(which(knowledge$s>H.max))==0 | length(which(knowledge$s_>H.max))==0)

cat("> ...Done!\n")

cat("> Saving data...\n")
save( "knowledge", file="generated_data.RData" )
cat("> ...DOne!\n")
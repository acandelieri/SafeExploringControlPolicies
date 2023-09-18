rm(list=ls()); graphics.off(); cat("\014")

zoomed.chart <- T
simulation <- T

m = 1470
c = 1005.4
M = 3600
R = 4.329 *10^-7

h.ref <- 18
h0 <- -30
dt <- 1/60
Tf <- 24
ts <- seq(0,Tf-dt,by=dt)


Kp=43; Ki=0.17; Kd=0; max.heat=70

# csi <- numeric(length(ts))+h0

# 1 hour
# csi <- h0 + 10 * cos(2*pi*ts/Tf) # +rnorm(length(ts),0,2)

# 24  hours
csi <- h0 + 10 * cos(2*pi*ts/Tf) # +rnorm(length(ts),0,2)

# plot( ts, csi, "l"  )

min.dev <- 0.5


house.sim.closed.loop <- function( max.heat=30, Kp, Ki, Kd, h.ref, h0, ts, dt, csi, m=1470,  c=1005.4, M=3600, R=4.329*10^-7 ) {
  h <- numeric(length(ts)+1)
  err <- numeric(length(ts))
  h[1] <- h0
  x_ <- numeric(length(ts))+NA
  
  for( i in 1:length(ts) ) {
    
    err[i] <- h.ref - h[i]
    
    # if( all(is.na(c(Kp,Ki,Kd))) ) {
    #   if( i==1 )
    #     cat("***** Simple controller applied! *****\n")
    #   if( err[i]>0.1 ) {
    #     x <- max.heat
    #   } else {
    #     if( err[i]<0.1 ) {
    #       x <- 0
    #     } else {
    #       if( i==1 ) {
    #         x <- 0
    #       } else {
    #         x <- x_[i-1]
    #       }
    #     }
    #   }
    # } else {
      
      x <- Kp * (err[i]) # Proportional
      x <- x + Ki * sum( err[1:i] ) # Integrative
      if(i>1)
        x <- x + Kd * ( err[i]-err[i-1] )/dt # Derivative
      
      x <- min( x, max.heat )
    # }
    
    x_[i] <- x
  
    dQg <- M*c*(x-h[i])
    dQl <- (h[i]-csi[i])/R
      
    h[i+1] <- h[i] + (dt/(m*c)) * ( dQg - dQl )
    
   
          
    if( i<20 ) {
      cat("\n****** t =",i,"*****\n")
      cat("err[t] =",err[i],"\n")
      cat("x =",x,"\n")
      cat("h[t] =",h[i],"\n")
      cat("dQg =",dQg,"\n")
      cat("dQl =",dQl,"\n")
      cat("h[t+dt] =",h[i+1],"\n")
    }
    
  }
  return( list(h=h,x=x_) )
}


res <- house.sim.closed.loop( max.heat=max.heat, Kp, Ki, Kd, h.ref, h0, ts, dt, csi )
h <- res$h
x <- res$x

# max elongation
elong = max(h) - h.ref
cat("\n\nMax elongation: ",elong,"\n")
elong.t = ts[which.max(h)]
cat(" @ t =",elong.t,"\n")
devs <- abs(h-h.ref)
cat("Min deviation from r.ref =", min(devs),"\n")

i <- 1
while( i<=length(devs) & max(devs[i:length(devs)])>min.dev )
  i <- i+1
min.dev.t <- ts[i]
cat("Deviation +/- ",min.dev,"C from h.ref @ t=",min.dev.t,"\n")

curr.mar <- par("mar")
par(mar=c(4.1,4.6,0.6,1.1))
par(mfrow=c(2,1))

if( zoomed.chart ) {
  plot( ts, h[-length(h)], "l", lwd=3, ylim=c(h.ref-2,h.ref+1),
        xlab="time [hours]", ylab="In-house temperature [°C]",
        cex.lab=1.5, cex.axis=1.5 ) #zoomed!
} else {
  plot( ts, h[-length(h)], "l", lwd=3, ylim=c(h0,h.ref+1),
        xlab="time [hours]", ylab="In-house temperature [°C]",
        cex.lab=1.5, cex.axis=1.5 )
}


# abline( h=h.ref+1, col="purple" )
# abline( h=h.ref-1, col="purple" )
# abline( h=h.ref+0.1, col="green3" )
# abline( h=h.ref-0.1, col="green3" )
# abline( h=h.ref, col="red", lty=2 )
# abline( v=ts[5*60], lty=2, col="grey" )

abline( h=h.ref, lty=2, lwd=3, col="blue" )
abline( h=h.ref+1, col="red", lty=1, lwd=2 )
abline( v=min.dev.t, lty=2, col="red", lwd=2 )
lines( ts, h[-length(h)], lwd=2, ylim=c(h0,h.ref+0.75) )

abline( h=h.ref+min.dev, lwd=2, col="green3" )
abline( h=h.ref-min.dev, lwd=2, col="green3" )
lines( c(-1.0,elong.t), rep(h.ref+elong,2), col="pink", lwd=2, lty=1 )

legend("bottom", legend=c("target temperature",
                          "allowed deviation from target",
                          "max allowed sovraelongation",
                          "max time for sovraelongation"),
       col=c("blue","green3","red", "red"), lty=c(2,1,1,2), lwd=3, cex=1.3 )
lines( ts, h[-length(h)], lwd=3 )


plot( ts, x, "s", lwd=3,
      xlab="time [hours]",
      ylab="Heater temperature [°C]",
      cex.lab=1.5, cex.axis=1.5 )
par(mfrow=c(1,1))
par(mar=curr.mar)

cat("Costo energetico:",sum(h*0.1),"euros\n")


if( !simulation )
  stop("simulation not required")
if( Tf!=24 )
  stop("Impossibile continuare!")

# data from github - prophet yosemite example R/Python
T.data <- read.delim2("temperature_data_yosemite.txt",sep=",",header=T)
T.data$y <- as.numeric(T.data$y)

last.index <- 17568
T.data <- T.data[1:last.index,]
T.data <- matrix(T.data$y,ncol=24*12,byrow=T)
ixs <- which( is.na(T.data),arr.ind=T )
T.data <- T.data[-ixs[,1],]
T.data <- T.data - 32 # conversion to Celsius

CSIs <- matrix(NA,nrow=nrow(T.data),ncol=24*60)
for( i in 1:nrow(T.data) ) {
  # interpolate to the minute
  tmp <- approx(seq(0,24-5/60,by=5/60),T.data[i,],n=60*24)
  CSIs[i,] <- tmp$y 
}

curr.mar <- par("mar")
par(mar=c(4.6,4.6,1.1,1.1) )
for( i in 1:nrow(T.data) ) {
  if( i==1 ) {
    plot( (1:ncol(CSIs))/60, CSIs[i,], type="l",
          ylim=c(min(CSIs),max(CSIs)), lwd=2,
          ylab="outside temperature [°C]", xlab="time [hours]",
          cex.lab=1.5, cex.axis=1.5, col="darkgray" )
  } else {
    lines( (1:ncol(CSIs))/60, CSIs[i,], lwd=2, col="darkgray" )
  }
}
par(mar=curr.mar)


Hs <- matrix(NA,nrow(CSIs),ncol(CSIs)+1 )
Xs <- matrix(NA,nrow(CSIs),ncol(CSIs) )

for( i in 1:nrow(CSIs) ) {
  csi <- CSIs[i,]
  h0 <- csi[1]
  res <- house.sim.closed.loop( max.heat=max.heat, Kp, Ki, Kd, h.ref, h0, ts, dt, csi )
  Hs[i,] <- res$h
  Xs[i,] <- res$x
}

curr.mar <- par("mar")
par(mar=c(4.6,4.6,1.1,1.1) )
for( i in 1:nrow(Hs) ) {
  if( i==1 ) {
    plot( ts, Hs[i,-ncol(Hs)], type="l", ylim=c(16,19),lwd=2, 
          xlab="time [hours]", ylab="In-house temperature [°C]",
          cex.lab=1.5, cex.axis=1.5, col="darkgray")
  } else {
    lines( ts, Hs[i,-ncol(Hs)], lwd=2, col="darkgray" )
  }
}
abline( h=h.ref, col="blue", lty=2, lwd=3 )
abline( h=h.ref+1, col="red", lty=1, lwd=3 )
abline( h=h.ref+min.dev, col="green3", lty=1, lwd=2 )
abline( h=h.ref-min.dev, col="green3", lty=1, lwd=2 )

abline( v=ts[30], col="red", lty=2, lwd=2 )
legend( "bottom", legend=c("target temperature",
                           "max allowed deviation form target",
                           "max allowed sovraelongation",
                           "max time for sovraelongation"),
        col=c("blue","green3","red","red"), lty=c(2,1,1,2), 
        lwd=3, cex=1.5 )

for( i in 1:nrow(Xs) ) {
  if( i==1 ) {
    plot( ts, Xs[i,], type="l", ylim=c(0,70),
          xlab="time [hours]", ylab="heater temperature [°C]",
          cex.axis=1.5, cex.lab=1.5, lwd=2, col="darkgray" )
  } else {
    lines( ts, Xs[i,], lwd=2, col="darkgray" )
  }
}

par( mar=curr.mar )

record <- list( Kp=Kp, Ki=Ki, Kd=Kd, Hs=Hs, Xs=Xs, CSIs=CSIs, h.ref=h.ref,
                m=m, c=c, M=M, R=R, Tf=Tf, dt=dt, ts=ts,
                house.sim.closed.loop=house.sim.closed.loop )
saveRDS(record,"RECORD.RDS")


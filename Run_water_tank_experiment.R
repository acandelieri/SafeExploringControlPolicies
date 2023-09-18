rm(list=ls()); cat("\014"); graphics.off()

library(DiceKriging)
library(plot3D)
library(nloptr)

set.seed(42)

# ******************************************************************************************
# System's technical parameters (must be the same as for data generation!!!)
# ******************************************************************************************
max.inflow = 10
H.min = 5
H.max = 50
prices = c(rep(1,6),rep(3,2),rep(2,10),rep(3,2),rep(2,2),rep(1,2)) # energy price tariff
eta = 0.85


# ******************************************************************************************
# Approach's parameters
# ******************************************************************************************
kernel.dyn = "matern3_2"  # kernel for the GP approximating the system's dynamics
kernel.cost = "gauss"     # kernel for the GP approximating the costs
beta.R = 1                # beta for the GP-LCB acquisition function (recommended value is 1)
beta.D = 4                # beta for the dealing with uncertainty on safety estimation

size = 500                # number of observations to consider for fitting each GP model
increase = 25             # number of observations iteratively added up to "size"

n.candidates = 10         # solutions to consider as canditates in the inner optimization

n.test.days = 30          # number of days sampled as test set 

visualize.approximation = F   # charts generation will slow down
visualize.optimization = F    # charts generation will slow down


retrainingAfterDays = 10   # number of days under new policy to retrain both the GP models



# ******************************************************************************************
# Functions
# ******************************************************************************************

# defining the objective function (regret) as approximated by the associated GP
obj.fun <- function( x, gp.cost, gp.dyn, t, h, beta.cost, beta.dyn, h.min, h.max ) {
  x = data.frame( t=t, h=h, a=x )
  pred = predict( gp.cost, x, "UK" )
  return( pred$mean-beta.cost*pred$sd )
}

# defining the safety constraint as approximated by the associated GP
safety.fun <- function( x, gp.dyn, gp.cost, t, h, beta.dyn, beta.cost, h.min, h.max ) {
  x = data.frame( t=t, h=h, a=x )
  pred = predict( gp.dyn, x, "UK" )
  
  constr1 = (pred$mean - beta.dyn * pred$sd) - h.min
  constr2 = h.max - (pred$mean + beta.dyn * pred$sd)
  
  return( c( constr1, constr2) ) # <-- negative if unsafe! (inequaliy constraint "g(x)<=0")
}





# ******************************************************************************************
# Main
# ******************************************************************************************

# Load experience
cat("> Loading safe experience...\n")
load("generated_data.RData")

# convert experience into a suitable dataframe
N = nrow(knowledge$s) * ncol(knowledge$s)
X = data.frame( t=numeric(N), # time
                h=numeric(N), # state
                a=numeric(N) # action
                )
h.next = numeric(nrow(X)) # next state
cost = numeric(nrow(X)) # regret

k = 1
for( i in 1:nrow(knowledge$s) ) {
  for( j in 1:ncol(knowledge$s) ) {
    X$t[k] = j
    X$h[k] = round(knowledge$s[i,j],4) # 4 decimal digits
    X$a[k] = round(knowledge$a[i,j],4) # 4 decimal digits
    
    h.next[k] = round(knowledge$s_[i,j],4) # 4 decimal digits
    cost[k] = round(knowledge$r[i,j],4) # 4 decimal digits
    k = k+1
  }
}

# use a subset of experience (incrementally sampled) to reduce computational burden of GP fitting
ixs.D = sample(1:nrow(X),increase,F) 
ixs.R = ixs.D
while( length(ixs.D)<size ) {
  
  cat("> Experience size for the GPs:",length(ixs.D),"\n")
  
  # fit the GP approximating the system's transition dynamic
  cat("> Fitting the GP approximating system's dynamics...\n")
  gp.dyn = km( design=X[ixs.D,], response=h.next[ixs.D], covtype=kernel.dyn, nugget.estim=T, control=list(trace=F) )
  yy = predict(gp.dyn,X[-ixs.D,],"UK")
  err = abs(yy$mean-h.next[-ixs.D])
  iii = order(err,decreasing=T)
  iii = iii[1:increase]
  ixs.D = c(ixs.D,setdiff(1:nrow(X),ixs.D)[iii])
  
  # fit the GP approximating the regret
  cat("> Fitting the GP approximating regret...\n")
  gp.cost = km( design=X[ixs.R,], response=cost[ixs.R], covtype=kernel.cost, nugget.estim=T, control=list(trace=F) )
  yy = predict(gp.cost,X[-ixs.R,],"UK")
  err = abs(yy$mean-h.next[-ixs.R])
  iii = order(err,decreasing=T)
  iii = iii[1:increase]
  ixs.R = c(ixs.R,setdiff(1:nrow(X),ixs.R)[iii])
  
}

# Plotting approximated transition dynamic (for each t=1,...,24) along with safe experience
if( visualize.approximation ) {
  par(mfrow=c(2,2))
  curr.mar <- par("mar")
  par(mar=c(4.6,4.6,2.6,1.6))
  cat("> Plotting approximations and experience...\n")
  for( tt in 1:24 ) {
    XX <- expand.grid( seq(5,H.max,length.out=30), seq(0,max.inflow,length.out=30) )
    XX <- cbind( rep(tt,30*30), XX )
  
    preds <- predict( gp.dyn, XX, "UK", checkNames=F )
  
    image2D( z=matrix(preds$mean,30,30,byrow=T),
             x=sort(unique(XX[,2])),
             y=sort(unique(XX[,3])),
             xlab="s(t-1)", ylab="a(t)",
             col=rev(cm.colors(100)),
             main=paste0("predicted s(t) for t=",tt),
             cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
    contour2D( z=matrix(preds$mean,30,30,byrow=T),
               x=sort(unique(XX[,2])),
               y=sort(unique(XX[,3])),
               nlevels=10, add=T, col="black" )
  
    points2D( x=X$h[which(X$t==tt)], y=X$a[which(X$t==tt)], pch=19, cex=1.5, col="blue", add=T )
  
    image2D( z=matrix(preds$sd,30,30,byrow=T),
             x=sort(unique(XX[,2])),
             y=sort(unique(XX[,3])),
             xlab="s(t-1)", ylab="a(t)",
             col=rev(cm.colors(100)),
             main=paste("uncertainty on s(t) prediction" ),
             cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
    contour2D( z=matrix(preds$mean,30,30,byrow=T),
               x=sort(unique(XX[,2])),
               y=sort(unique(XX[,3])),
               nlevels=10, add=T, col="black" )
  
    points2D( x=X$h[which(X$t==tt)], y=X$a[which(X$t==tt)], pch=19, cex=1.5, col="blue", add=T )
    
  }
  par(mar=curr.mar)
  par(mfrow=c(1,1))
}


# retieve the (actual) water demand from experience data
demand <- knowledge$s + knowledge$a - knowledge$s_

# # ***** uncomment if you want to visualize the water demand ******************
# curr.mar = par("mar"); par(mar=c(4.1,5.6,1.1,1.1) )
# plot( demand[1,], type="l", ylim=c(min(demand),1.1*max(demand)),
#       ylab=expression(paste("outflow [",m^3,"]")), xlab="hours",
#       cex.lab=2, cex.axis=2)
# for( i in 1:nrow(demand) )
#   lines( demand[i,] )
# par(mar=curr.mar)
# # ****************************************************************************



#*********************************************************************************
# EXPERIMENT 1 - Safe-exploring GP approach applied on a set of test days
#*********************************************************************************

cat("\n> ***** Starting safe-exploring from safe-experience approach *****\n")

actions = matrix(NA,n.test.days,ncol(knowledge$s))
forced = matrix(F,n.test.days,ncol(knowledge$s))
tank.levels = matrix(NA,n.test.days,ncol(knowledge$s))
incurred.cost = matrix(NA,n.test.days,ncol(knowledge$s))

test.days = sample(1:(nrow(X)/24),n.test.days)

for( day in test.days )  {

  cat("\014...day",which(test.days==day),"out of",n.test.days,"...\n")  
  tt = 1
  h = knowledge$s[day,1] 
    
  
  while( tt<=24 ) {
  
    tank.levels[which(test.days==day),tt] = h
    
    x = data.frame(t=tt,h=h)
    
    # initially (safe) solution
    opt.details = data.frame( x0=round(runif(n.candidates,0,max.inflow),4),
                              is.x0.safe=logical(n.candidates),
                              xf=numeric(n.candidates),
                              yf=numeric(n.candidates),
                              is.xf.safe=logical(n.candidates) )
    for( i in 1:n.candidates ) {
      is.safe = safety.fun( opt.details$x0[i], gp.dyn, gp.cost, tt, h, beta.D, beta.R, H.min, H.max )
      opt.details$is.x0.safe[i] = all(is.safe>0)
    
      # res <- nloptr( x0=opt.details$x0[i], eval_f=obj.fun,
      #                lb=0, ub=max.inflow,
      #                eval_g_ineq=safety.fun,
      #                #eval_g_ineq=full.safety.fun,
      #                opts=list( #algorithm="NLOPT_LN_ISRES",
      #                           #algorithm="NLOPT_LN_COBYLA",
      #                           algorithm="NLOPT_LN_AUGLAG_EQ",
      #                           tol_constraints_ineq=rep(0.0,2),
      #                           xtol_rel=10^-4,
      #                           maxeval=500,
      #                           print_level=0 ),
      #                gp.cost=gp.cost,
      #                gp.dyn=gp.dyn,
      #                t=tt,
      #                h=h,
      #                beta.cost=beta.R,
      #                beta.dyn=beta.D,
      #                h.min=H.min,
      #                h.max=H.max )
      
      if( opt.details$is.x0.safe[i] ) {
        res = cobyla( x0=opt.details$x0[i],
                      fn=obj.fun,
                      lower=0,
                      upper=max.inflow,
                      hin=safety.fun,
                      nl.info=FALSE,
                      control=list(xtol_rel = 1e-8, maxeval = 1000),
                      gp.cost=gp.cost,
                      gp.dyn=gp.dyn,
                      t=tt,
                      h=h,
                      beta.cost=beta.R,
                      beta.dyn=beta.D,
                      h.min=H.min,
                      h.max=H.max )
        
        opt.details$xf[i] = res$par
        opt.details$yf[i] = res$value
        is.safe = safety.fun( res$par, gp.dyn, gp.cost, tt, h, beta.D, beta.R, H.min, H.max )
        opt.details$is.xf.safe[i] = all(is.safe>0)
      } else {
        opt.details$xf[i] = NA
        opt.details$yf[i] = NA
        opt.details$is.xf.safe[i] = F
      }
      
    }
    
    ixs = which( opt.details$is.x0.safe & opt.details$is.xf.safe )
    opt.details = opt.details[ixs,]
    
    if( nrow(opt.details)==0 ) { # there are not safe solutions!
      cat("* forcing safe-by-design policy\n")
      a = min( c(max.inflow,H.max-h) )
      forced[which(test.days==day),tt] = T
    } else {
      a = opt.details$xf[which.min(opt.details$yf)]
    }
      
    a = round(a,4)
    cat("h=",h," + a=",a," --> h'=",h + a - demand[day,tt],"\n",sep="")
    
    actions[which(test.days==day),tt] = a

    
    # Just to visualize... ------------------------------------
    if( visualize.optimization ) {
      XX = expand.grid( t=tt, h=h, a=seq(0,max.inflow,by=0.01) )
      pp = predict( gp.cost, XX, "UK" )
      bb = predict( gp.dyn, XX, "UK" )
      par(mfrow=c(2,1))
      curr.mar = par("mar")
      par(mar=c(4.6,4.6,2.1,1.1))
      plot( seq(0,max.inflow,by=0.01), pp$mean, type="l", lwd=3,
            ylim=c(min(pp$mean-beta.R*pp$sd),max(pp$mean+beta.R*pp$sd)),
            xlab=expression(paste("inflow [",m^3,"]")), ylab="Expected cost [€]",
            main=paste0("day=",day," - t=",tt), col="blue",
            cex.main=1.5, cex.axis=1.5, cex.lab=1.5 )
      polygon( c(seq(0,max.inflow,by=0.01),rev(seq(0,max.inflow,by=0.01))),
               c(pp$mean-beta.R*pp$sd,rev(pp$mean+beta.R*pp$sd)),
               col=adjustcolor("blue",alpha.f=0.2), border=F ) 
      lines( seq(0,max.inflow,by=0.01), pp$mean-beta.R*pp$sd, lty=2, col="blue", lwd=3 )
      if(nrow(opt.details)==0) {
        abline(v=a,lwd=2,lty=2,col="red")
      } else {
        abline(v=a,lwd=2,lty=2)
      }
      plot( seq(0,max.inflow,by=0.01), bb$mean, type="l", lwd=3, ylim=c(H.min-5,H.max+5),
            xlab=expression(paste("inflow [",m^3,"]")),ylab="next (estimated) tank level",
            cex.axis=1.5, cex.lab=1.5, col="blue")
      polygon( c(seq(0,max.inflow,by=0.01),rev(seq(0,max.inflow,by=0.01))),
               c(bb$mean-beta.D*bb$sd,rev(bb$mean+beta.D*bb$sd)),
               col=adjustcolor("blue",alpha.f=0.2), border=F ) 
      lines( seq(0,max.inflow,by=0.01), bb$mean+beta.D*bb$sd, lty=2, col="blue", lwd=3 )
      lines( seq(0,max.inflow,by=0.01), bb$mean-beta.D*bb$sd, lty=2, col="blue", lwd=3 )
      abline( h=H.min, col="red", lwd=3 )
      abline( h=H.max, col="red", lwd=3 )
      if(nrow(opt.details)==0) {
        abline(v=a,lwd=2,lty=2,col="red")
      } else {
        abline(v=a,lwd=2,lty=2)
      }
      par(mar=curr.mar)
      par(mfrow=c(1,1))
    }
    # ---------------------------------------------------------
    
    h = h + a - demand[day,tt]
    incurred.cost[which(test.days==day),tt] <- prices[tt] * (1/eta) * a^(1/3)

    tt <- tt <- tt+1
  }
  
}

# Visualizing results

cat("\n")
I1 <- which( tank.levels>=H.max, arr.ind=T )
I1 <- unique(I1[,1])
cat("UP violations:",length(I1),"\n")

I2 <- which( tank.levels<=H.min, arr.ind=T )
I2 <- unique(I2[,1])
cat("DOWN violations:",length(I2),"\n")

safe.costs <- apply( knowledge$r[test.days,], 1, sum )
new.costs <- apply( incurred.cost, 1, sum )

plot( safe.costs, type="l", lwd=3, ylim=c(min(safe.costs,new.costs),max(safe.costs,new.costs)) )
points( safe.costs, pch=19 )
lines( new.costs, col="blue", lwd=3 )
points( new.costs, pch=19, col="blue" )
points( c(I1,I2), new.costs[c(I1,I2)], pch=1, col="red", lwd=2, cex=2 )

curr.mar <- par("mar"); par(mar=c(2.1,5.1,1.1,1.1) )
boxplot( safe.costs, new.costs, col=c("pink","skyblue"),
         names=c("safe initial policy","new policy"),
         ylab="Incurred cost [€/day]", cex.lab=2, cex.axis=2, lwd=2 ) 
par(mar=curr.mar)

cat("\nSafe policy:\t average cost=",round(mean(safe.costs),2),"€ (+/-",round(sd(safe.costs),2),"€)\n",sep="")
cat("New policy:\t average cost=",round(mean(new.costs),2),"€ (+/-",round(sd(new.costs),2),"€)\n",sep="")

psafe = density(safe.costs,from=min(safe.costs,new.costs),to=max(safe.costs,new.costs),n=100)
pnew = density(new.costs,from=min(safe.costs,new.costs),to=max(safe.costs,new.costs),n=100)
plot( psafe$x, psafe$y, type="l", col="red", lwd=3,
      xlim=range(psafe$x,pnew$x), ylim=range(psafe$y,pnew$y),
      xlab="Incurred cost [€/day]", ylab="Probability density",
      cex.lab=1.3, cex.axis=1.3 )
polygon( c(psafe$x,rev(psafe$x)), c(psafe$y,numeric(100)), col=adjustcolor("red",alpha.f=0.2), border=F )
polygon( c(pnew$x,rev(pnew$x)), c(pnew$y,numeric(100)), col=adjustcolor("blue",alpha.f=0.2), border=F )
lines( pnew$x, pnew$y, col="blue", lwd=3 )
legend("topleft",legend=c("Initial safe policy","New safe policy"),col=c("red","blue"),lwd=3,cex=1.3)





#**************************************************************************************
# EXPERIMENT 2 - Long-run: safe-exploring GP over a year from an initial tank level
#**************************************************************************************

cat("\n> ***** Starting Long-run experiment *****\n")

# sampling an initial tank level from the historical data 
h0 <- sample(knowledge$s[,1],1)
h = h0 

actions.safe = actions.sgp = NULL
forced = NULL
tank.safe = tank.sgp = h
cost.safe = cost.sgp = NULL

cost.upds = dyn.upds = 0

# IMPORTANT: the beta.D must be increased because the approach will use less water and the initial
# tank level for the next day could be significantly different from any historical ones (likewise
# for the tank level during other hours of the day). Thus, we must be more conservative in safety
# estimation. From empirical evidence, increasing beta.D from 4 to 5 is sufficient to avoid safety
# violations also in the long-run experiment.

beta.D = 5 


X_ = h.next_ = cost_ = NULL

for( day in 1:nrow(demand) )  {

  cat("\014...day",day,"out of",nrow(demand),"...\n")
  tt <- 1
  
  if( day==retrainingAfterDays+1 ) {
    cat("\n\n\t***** Retraining GPs *****\n\n")
    gp.dyn = km( design=rbind(gp.dyn@X,X_), response=c(gp.dyn@y,h.next_),
                 covtype=kernel.dyn, nugget.estim=T, control=list(trace=F) )
    gp.cost = km( design=rbind(gp.cost@X,X_), response=c(gp.cost@y,cost_),
                  covtype=kernel.cost, nugget.estim=T, control=list(trace=F) )
  }

  while( tt<=24 ) {

    h = tank.sgp[length(tank.sgp)]
    
    x = data.frame(t=tt,h=h)

    # initially (safe) solution
    opt.details = data.frame( x0=round(runif(n.candidates,0,max.inflow),4),
                              is.x0.safe=logical(n.candidates),
                              xf=numeric(n.candidates),
                              yf=numeric(n.candidates),
                              is.xf.safe=logical(n.candidates) )
    for( i in 1:n.candidates ) {
      is.safe = safety.fun( opt.details$x0[i], gp.dyn, gp.cost, tt, h, beta.D, beta.R, H.min, H.max )
      opt.details$is.x0.safe[i] = all(is.safe>0)

      if( opt.details$is.x0.safe[i] ) {
        res = cobyla( x0=opt.details$x0[i],
                      fn=obj.fun,
                      lower=0,
                      upper=max.inflow,
                      hin=safety.fun,
                      nl.info=FALSE,
                      control=list(xtol_rel = 1e-8, maxeval = 1000),
                      gp.cost=gp.cost,
                      gp.dyn=gp.dyn,
                      t=tt,
                      h=h,
                      beta.cost=beta.R,
                      beta.dyn=beta.D,
                      h.min=H.min,
                      h.max=H.max )

        opt.details$xf[i] = res$par
        opt.details$yf[i] = res$value
        is.safe = safety.fun( res$par, gp.dyn, gp.cost, tt, h, beta.D, beta.R, H.min, H.max )
        opt.details$is.xf.safe[i] = all(is.safe>0)
      } else {
        opt.details$xf[i] = NA
        opt.details$yf[i] = NA
        opt.details$is.xf.safe[i] = F
      }

    }

    ixs = which( opt.details$is.x0.safe & opt.details$is.xf.safe )
    opt.details = opt.details[ixs,]

    if( nrow(opt.details)==0 ) { # there is not any safe solution
      cat("* forcing safe-by-design policy\n")
      a = min( max.inflow, H.max-h )
      forced = c(forced,T)
    } else {
      a = opt.details$xf[which.min(opt.details$yf)]
      forced = c(forced,F)
    }

    a = round(a,4)
    cat("h=",h," + a=",a," --> h'=",h + a - demand[day,tt],"\n",sep="")

    actions.sgp = c(actions.sgp,a)

    h_ = h + a - demand[day,tt]
    tank.sgp = c(tank.sgp,h_)

    cost.sgp = c(cost.sgp, prices[tt] * (1/eta) * a^(1/3))

    if( day<=retrainingAfterDays ) {
      X_ = rbind( X_, data.frame( t=tt, h=h, a=a ) )
      h.next_ = c(h.next_,h_)
      cost_ = c(cost_,prices[tt] * (1/eta) * a^(1/3))
    }
    
    
    #*******************************
    # ORIGINAL SAFE POLICY
    #*******************************

    h = tank.safe[length(tank.safe)]

    a = min( max.inflow, H.max-h )
    a = round(a,4)

    actions.safe = c(actions.safe,a)
    h_ = h + a - demand[day,tt]

    tank.safe = c(tank.safe, h_)
    cost.safe = c(cost.safe, prices[tt] * (1/eta) * a^(1/3))

    #*******************************

    tt = tt+1

  }

}

# Visualizing results

cat("\n")
I1.safe <- which( tank.safe>=H.max )
I2.safe <- which( tank.safe<=H.min )
if( length(I1.safe)>0 || length(I2.safe)>0 )
  cat("\n!!!!! Safety violations occurred with  safe-by-design !!!!!\n")
I1.sgp <- which( tank.sgp>=H.max )
I2.sgp <- which( tank.sgp<=H.min )
if( length(I1.sgp)>0 || length(I2.sgp)>0 )
  cat("\n!!!!! Safety violations occurred with safe-exploration !!!!!\n")

# time-serie plot for tank level
plot( tank.safe, type="l", lwd=3, ylim=c(min(H.min,tank.safe,tank.sgp),max(H.max,tank.safe,tank.sgp)), col="red" )
points( tank.safe, pch=19, col="red" )
lines( tank.sgp, col="blue", lwd=3 )
points( tank.sgp, pch=19, col="blue" )
abline( h=c(H.min,H.max), lwd=3, lty=2 )

# probability density for tank level
psafe = density(tank.safe,from=min(H.min,tank.safe,tank.sgp),to=max(H.max,tank.safe,tank.sgp),n=100)
pnew = density(tank.sgp,from=min(H.min,tank.safe,tank.sgp),to=max(H.max,tank.safe,tank.sgp),n=100)
plot( psafe$x, psafe$y, type="l", col="red", lwd=3,
      xlim=range(H.min,H.max,psafe$x,pnew$x), ylim=c(1,1.3)*range(psafe$y,pnew$y),
      xlab="Tank level", ylab="Probability density",
      cex.lab=1.3, cex.axis=1.3 )
polygon( c(psafe$x,rev(psafe$x)), c(psafe$y,numeric(100)), col=adjustcolor("red",alpha.f=0.2), border=F )
polygon( c(pnew$x,rev(pnew$x)), c(pnew$y,numeric(100)), col=adjustcolor("blue",alpha.f=0.2), border=F )
lines( pnew$x, pnew$y, col="blue", lwd=3 )
abline(v=c(H.min,H.max),lwd=3,lty=2)
legend("top",legend=c("Initial safe policy","New safe policy"),col=c("red","blue"),lwd=3,cex=1.3)


# boxplot and probability density for costs
cost.safe = apply(matrix(cost.safe,ncol=24,byrow=T),1,sum)
cost.sgp = apply(matrix(cost.sgp,ncol=24,byrow=T),1,sum)

curr.mar <- par("mar"); par(mar=c(2.1,5.1,1.1,1.1) )
boxplot( cost.safe, cost.sgp, col=c("pink","skyblue"),
         names=c("safe initial policy","new policy"),
         ylab="Incurred cost [€/day]", cex.lab=2, cex.axis=2, lwd=2 ) 
par(mar=curr.mar)


psafe = density(cost.safe,from=min(cost.safe,cost.sgp),to=max(cost.safe,cost.sgp),n=100)
pnew = density(cost.sgp,from=min(cost.safe,cost.sgp),to=max(cost.safe,cost.sgp),n=100)
plot( psafe$x, psafe$y, type="l", col="red", lwd=3,
      xlim=range(psafe$x,pnew$x), ylim=c(1,1.3)*range(psafe$y,pnew$y),
      xlab="Incurred costs [€/day]", ylab="Probability density",
      cex.lab=1.3, cex.axis=1.3 )
polygon( c(psafe$x,rev(psafe$x)), c(psafe$y,numeric(100)), col=adjustcolor("red",alpha.f=0.2), border=F )
polygon( c(pnew$x,rev(pnew$x)), c(pnew$y,numeric(100)), col=adjustcolor("blue",alpha.f=0.2), border=F )
lines( pnew$x, pnew$y, col="blue", lwd=3 )
legend("topleft",legend=c("Initial safe policy","New safe policy"),col=c("red","blue"),lwd=3,cex=1.3)


# time series for daily costs
plot( cost.safe, type="l", lwd=3, col="red", ylim=range(cost.safe,cost.sgp),
      xlab="days", ylab="costs [€/day]",
      cex.axis=1.3, cex.lab=1.3)
lines( cost.sgp, lwd=3, col=adjustcolor("blue",alpha.f=0.5) )
legend("bottomright",legend=c("Initial safe policy","New safe policy"),col=c("red",adjustcolor("blue",alpha.f=0.5)),lwd=3,cex=1.3)



# boxplot and probability density for actions
actions.safe = apply(matrix(actions.safe,ncol=24,byrow=T),1,sum)
actions.sgp = apply(matrix(actions.sgp,ncol=24,byrow=T),1,sum)


curr.mar <- par("mar"); par(mar=c(2.1,5.1,1.1,1.1) )
boxplot( actions.safe, actions.sgp, col=c("pink","skyblue"),
         names=c("safe initial policy","new policy"),
         ylab="Water inflow", cex.lab=2, cex.axis=2, lwd=2 ) 
par(mar=curr.mar)


psafe = density(actions.safe,from=min(actions.safe,actions.sgp),to=max(actions.safe,actions.sgp),n=100)
pnew = density(actions.sgp,from=min(actions.safe,actions.sgp),to=max(actions.safe,actions.sgp),n=100)
plot( psafe$x, psafe$y, type="l", col="red", lwd=3,
      xlim=range(psafe$x,pnew$x), ylim=range(psafe$y,pnew$y),
      xlab="actions level", ylab="Probability density",
      cex.lab=1.3, cex.axis=1.3 )
polygon( c(psafe$x,rev(psafe$x)), c(psafe$y,numeric(100)), col=adjustcolor("red",alpha.f=0.2), border=F )
polygon( c(pnew$x,rev(pnew$x)), c(pnew$y,numeric(100)), col=adjustcolor("blue",alpha.f=0.2), border=F )
lines( pnew$x, pnew$y, col="blue", lwd=3 )
legend("topleft",legend=c("Initial safe policy","New safe policy"),col=c("red","blue"),lwd=3,cex=1.3)


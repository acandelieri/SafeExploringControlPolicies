rm(list=ls()); graphics.off(); cat("\014")

library(nloptr)
library(DiceKriging)

library(foreach)
library(doParallel)
registerDoParallel(cores=6)




n.candidates <- 100
max.heat <- 70
dt <- 1/60
kernel.D = "matern3_2"
kernel.R = "gauss"
beta.D <- 6
beta.R <- 1

h.ref <- 18
min.dev <- 0.5
intemp.min <- h.ref - min.dev
intemp.max <- h.ref + min.dev
time.min <- 30 * dt # 30 minutes
t.elong <- 19

m = 1470
c = 1005.4
M = 3600
R = 4.329 * 10^-7


# *****************************************************************************
# Functions 
# *****************************************************************************

compute.cost <- function( ts, dt, xs, tariff ) {
  COST <- 0
  for( i in 2:length(ts) ) {
    hour <- ts[i]
    ix <- which( tariff$time<=hour )
    ix <- rev(ix)[1]
    COST <- COST + tariff$cost[ix] * xs[i] * dt
  }
  return(COST)
}

# defining the objective function (regret) as approximated by the associated GP
obj.fun <- function( x, gp.cost, gp.dyn, t, aprev, h, beta.cost, beta.dyn, intemp.min, intemp.max, time.min, t.elong ) {
  x <- data.frame( aprev=aprev, h=h, a=x )
  pred <- predict( gp.cost, x, "UK" )
  return( pred$mean-beta.cost*pred$sd )
}

# defining the safety constraint as approximated by the associated GP
safety.fun <- function( x, gp.cost, gp.dyn, t, aprev, h, beta.cost, beta.dyn, intemp.min, intemp.max, time.min, t.elong ) {
  x <- data.frame( aprev=aprev, h=h, a=x )
  pred <- predict( gp.dyn, x, "UK" )
  if( t<=time.min ) {
    constr1 <- (pred$mean - beta.dyn * pred$sd) - intemp.min
    constr2 <- 1
    constr3 <- t.elong - (pred$mean + beta.dyn * pred$sd)
  } else {
    constr1 <- (pred$mean - beta.dyn * pred$sd) - intemp.min
    constr2 <- intemp.max - (pred$mean + beta.dyn * pred$sd)
    constr3 <- 1
  }
  pred <- predict( gp.cost, x, "UK" )
  constr.aux = pred$mean - beta.cost * pred$sd 
  return( c( constr1, constr2, constr3, constr.aux ) ) # <-- negative if unsafe!
}

aux.fun <- function( x, gp.dyn, aprev, h, h.ref, tol ) {
  x = data.frame( aprev=aprev, h=h, a=x )
  yp_ = predict( gp.dyn, x, "UK" )
  aux.fun = (yp_$mean - h.ref)^2
}

aux.constr <- function( x, gp.dyn, aprev, h, h.ref, tol ) {
  x = data.frame( aprev=aprev, h=h, a=x )
  yp_ = predict( gp.dyn, x, "UK" )
  aux.constr = tol - yp_$sd
}






# *****************************************************************************
# Data: safe experience
# *****************************************************************************

cat("> Reading data...\n")
record <- readRDS("RECORD.RDS")



tariff <- list( time=c(0,7,10,17,21), # inizio fasce
                cost=c(1,10,5,10,1) # prezzo associato
)

cat("> Computing costs on experience...\n")

daily.cost <- numeric(nrow(record$Xs)) 
for( i in 1:length(daily.cost) )
  daily.cost[i] <- compute.cost( record$ts, record$dt, record$Xs[i,], tariff )



forEachRes <- foreach( day = 1:length(daily.cost) ) %dopar% {
  library(nloptr)
  library(DiceKriging)
  
  cat("\014")
  cat("[ DAY :",day,"]\n\n")
  
  Hs = record$Hs[day,1]
  CSIs = record$CSIs[day,]
  As = NULL
  
  info = NULL
  
  incurred.costs = NULL
  
  for( i in 1:length(CSIs) ) {
    
    h = Hs[i] # current in-house temperature
    
    
    # Modelling system's dynamics
    gp.dyn = NULL
    if( i == 1 ) { aprev = 0 } else { aprev=record$Xs[-day,i-1] }
    try( expr = (gp.dyn = km( design = data.frame(aprev=aprev, h=record$Hs[-day,i], a=record$Xs[-day,i]),
                              response = record$Hs[-day,i+1], multistart=5, 
                              covtype = kernel.D, nugget.estim = T,
                              control = list(trace=0) ) ),
         silent = T )
    
    
    # Modelling action-related cost
    t = (i-1)*dt
    ix = which( tariff$time<=t )
    ix <- rev(ix)[1]
    gp.cost = NULL
    try( expr = (gp.cost = km( design = data.frame(aprev=aprev, h=record$Hs[-day,i], a=record$Xs[-day,i]),
                              response = record$Xs[-day,i] * tariff$cost[ix] * dt,
                              multistart = 5,
                              covtype = kernel.R, nugget.estim = T,
                              control = list(trace=0) ) ),
         silent = T )
  
  
  
    #******************************************************************
    # Solving the optimal control problem 
    #******************************************************************
    
    solveAuxProbl = F
    
    
    if( !is.null(gp.dyn) && !is.null(gp.cost) ) {
      
      # both system's dynamics and action-related cost models are available
      
      
      # solving the constrained optimization problem
      
      x0s = runif( n.candidates, 0, max.heat )
      
      # evaluating safety for initial candidate solutions
      safe.sols = logical(n.candidates)
      for( kk in 1:length(x0s) ) {
        is.safe = safety.fun( x=x0s[kk], gp.cost=gp.cost, gp.dyn=gp.dyn, t=(i-1)*dt,
                              aprev=record$Xs[day,i], h=h,
                              beta.cost=beta.R, beta.dyn=beta.D,
                              intemp.min=intemp.min, intemp.max=intemp.max,
                              time.min=time.min, t.elong= t.elong )
        safe.sols[kk] = all(is.safe>0)
      }
      
      
      if( length(which(safe.sols==T))==0 ) { # no safe candidates!
        
        solveAuxProbl = T  
        cat("* No safe candidate solutions for optimization!\n")
        
      } else { # at least one initial safe candidate solution is available!
        
        cat("> Optimizing cost constrained to safety...\n")
        
        best.res = NULL
        
        for( j in which(safe.sols==T) ) {
          res = cobyla( x0=x0s[j], fn=obj.fun,
                        lower=0, upper=max.heat,
                        hin=safety.fun,
                        gp.cost=gp.cost,
                        gp.dyn=gp.dyn,
                        t=t, aprev=record$Xs[day,i-1], h=h,
                        beta.cost=beta.R, beta.dyn=beta.D,
                        intemp.min=intemp.min, intemp.max=intemp.max,
                        time.min=time.min, t.elong=t.elong )
          
          # check for safety of the optimal solution!
          is.safe = safety.fun(res$par,gp.cost,gp.dyn,t,record$Xs[day,i],h,
                               beta.R,beta.D,intemp.min,intemp.max,time.min,t.elong)
          
          # updating the best safe solution so far
          if( all(is.safe>0) && (is.null(best.res) || best.res$value>res$value) ) {
            best.res = res
          } 
          
        }
        
        if( is.null(best.res) ) { # no safe solutions after constrained optimization!
          solveAuxProbl = T 
          cat("* No safe solutions after constrained optimization!\n")
        } else {
          info = c(info,"optimized")
        }
        
      }
      
    } else {
      
      if( is.null(gp.dyn) ) { # no approximations for safety are available
        cat("* No GP modeling the system's dynamics!\n")
      } else {
        solveAuxProbl = T
      }
    }
    
    
    if( solveAuxProbl ) {
      
      cat("> Optimizing auxiliary problem...\n")
  
      rm("best.res")
      
      if( t<time.min ) {
        h_ = h.ref
        htol = 0.2
      } else {
        h_ = h
        htol = 0.005
      }
      
      # searching for initial feasible solutions of the auxiliary problem
      x0s = runif(100,0,max.heat)
      feas0 = aux.constr( x0s, gp.dyn, record$Xs[day,i-1], h, h_, htol )
      x0s = x0s[which(feas0>0)]
      
      if( length(x0s)>0 ) {
        best.res = NULL
        
        for( j in 1:length(x0s) ) {
          res = cobyla( x0=x0s[j], fn=aux.fun,
                        lower=0, upper=max.heat,
                        hin=aux.constr,
                        gp.dyn=gp.dyn,
                        aprev=record$Xs[day,i-1], h=h, h.ref=h_, tol=htol )
          
          # check for feasibility of the optimal solution!
          isFeas = aux.constr( x=res$par, gp.dyn=gp.dyn, aprev=record$Xs[day,i-1], h=h, h.ref=h_, tol=htol )
          
          # update best feasible solution so far
          if( is.null(best.res) || best.res$value>res$value )
            best.res = res
        }
        
        if( is.null(best.res) ) {
          rm("best.res")
          cat("* Auxiliary problem cannot  be solved!\n")
        } else {
          info = c(info,"auxiliary")
        }
      } else {
        cat("* No initial feasible solutions for the auxiliary problem!\n")
      }
    }
    
    
    if( !exists("best.res") ) {
      cat("** Ad-hoc solution is applied...\n")
      
      # actual action
      a = record$Xs[day,i]
      if( h>h.ref & h<intemp.max ) {
        a = 0.5*a
      } else {
        if( h>intemp.max )
          a = 0.3*a
      }
      
      info = c(info,"ad-hoc")
      
    } else  {
      a = best.res$par
      rm(best.res)
    }
    
    
    
    a = round(a,4)
    
    
    dQg = M * c * (a-h)     # gain
    dQl = (h-CSIs[i]) / R   # loss
    
    h.next = h + (dt/(m*c)) * (dQg-dQl)
    
    cat("i=",i," --> Time ",t,": in-temp=",h," + a=",a," --> next in-temp=",h.next,"\n",sep="")
    
    As = c(As,a)
    
    ix = which( tariff$time<=t )
    ix = rev(ix)[1]

        incurred.costs = c(incurred.costs,tariff$cost[ix] * a * dt)
    
    Hs = c(Hs,h.next)
    
    
  }

  
  # returning for each
  
  forEachRes = list( As=As, Hs=Hs, costs=incurred.costs )
  
}

stopImplicitCluster()

save.image(paste0("heating_workspace_",as.character(Sys.Date()),".RData"))






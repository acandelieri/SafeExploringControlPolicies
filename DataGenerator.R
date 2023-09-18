tank.episodes.generator <- function( N, H, X, mu.h, sd.h, mu.xi, sd.xi, policy=safe.policy, prices, eta, seed=42 ) {
  
  # checks
  stopifnot( length(H)==2 & H[1]<H[2] )
  stopifnot( length(X)==2 & X[1]<X[2] )
  stopifnot( length(mu.h)==1 )
  stopifnot( length(sd.h)==1 )
  stopifnot( length(mu.xi)==length(sd.xi) )
  stopifnot( length(mu.xi)==length(prices) )
  stopifnot( length(eta)==1)
  
  
  # ********************************************************
  # generating episodes 
  # ********************************************************
  
  set.seed(seed)  
  
  # only 1 initial state, then every new episode starts from the final state of the previous
  h <- runif(1,H[1],H[2]) 
  
  # uncertainy on xi, for every episode and epoch
  T. <- length(mu.xi)
  Xi <- matrix(NA,N,T.)
  for( t in 1:T. ) 
    Xi[,t] <- rnorm(N,mu.xi[t],sd.xi[t])
  
  
  # ********************************************************
  # applying policy and storing observations
  # ********************************************************
  
  s <- matrix(NA,N,T.)
  a <- matrix(NA,N,T.)
  s_ <- matrix(NA,N,T.)
  r <- matrix(NA,N,T.)
  
  for( i in 1:N )
    for( t in 1:T. ) {
      s[i,t] <- h
      x <- do.call( policy, list( h=h, h.max=H[2], x.max=X[2] ) )
      a[i,t] <- x
      h <- h + x - Xi[i,t]
      if( h<0 ) {
        cat("Giorno:",i,"\n")
        cat("Istante temporale:",t,"\n")
        cat("Livello iniziale:",s[i,t],"\n")
        cat("Erogato:",a[i,t],"\n")
        cat("Domanda:",Xi[i,t],"\n")
        cat("Livello finale:",h,"\n")
        stop("ERRORE GRAVE: SERBATOIO VUOTO!")
      }
      s_[i,t] <- h
      r[i,t] <- prices[t] * (1/eta) * x^(1/3)
      if( is.na(r[i,t]) ) {
        print(prices[t])
        print(1/eta)
        print(x)
        print(x^(1/3))
        print(r[i,t])
        stop("ERRORE!!!!!")
      }
        
    }
  
  return( list( s=s, a=a, s_=s_, r=r ) )
}

# safe (but not-optimal) policy 
safe.policy <- function( h, h.max, x.max ) {
  safe.policy <- min(x.max, max(0,h.max - h) )
}
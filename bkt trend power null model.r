model{
  
  # Abundance model
  
    for(i in 1:nSites){
      for(j in 1:nYears){
        N[i,j] ~ dpois(lambda[i,j])
        log(lambda[i,j]) <- mu + trend*(j-1) + site.ran[i] + year.ran[j] + eps[i,j]
      }
    }
  
    ## priors
    mu ~ dnorm(0, 0.01)      # overall intercept
    trend ~ dnorm(0, 0.01)   # linear trend
  
    for (i in 1:nSites){
      site.ran[i] ~ dnorm(0, tau.site)  # Random site effect
    }
    tau.site <- pow(sd.site, -2)
    sd.site ~ dunif(0,5)
    sd2.site <- pow(sd.site, 2)
  
    for (j in 1:nYears){
      year.ran[j] ~ dnorm(0, tau.year)  # Random year effect
    }
    tau.year <- 1/(sd.year*sd.year) 
    sd.year ~ dunif(0,5)
    sd2.year <- pow(sd.year, 2)
  
    for(i in 1:nSites){
      for(j in 1:nYears){
        eps[i,j] ~ dnorm(0, tau)  # Over-dispersion
      }
    }
    tau <- pow(sigma, -2)
    sigma ~ dunif(0, 3)
    sigma2 <- pow(sigma, 2)
  
  # Detection model
    
    for(i in 1:nSites){
      for(j in 1:nYears){
        y[i,j,1] ~ dbin(p[i,j], N[i,j])
        y[i,j,2] ~ dbin(p[i,j], N[i,j]-y[i,j,1])
        y[i,j,3] ~ dbin(p[i,j], N[i,j]-y[i,j,1]-y[i,j,2])
        
        p[i,j] <- 1/(1 + exp(-lp.lim[i,j]))
        lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
        lp[i,j] <- p.mu + p.b*p.X[i,j] + p.eps[i,j] 
      }
    }
  
    ## priors  
    p.mean ~ dunif(0,1)
    p.mu <- log(p.mean/(1-p.mean))
    p.b ~ dnorm(0, 0.37)
    
    for(i in 1:nSites){
      for(j in 1:nYears){
        p.eps[i,j] ~ dnorm(0, p.tau)  # Over-dispersion
      }
    }
    p.tau <- pow(p.sigma, -2)
    p.sigma ~ dunif(0, 1)
    p.sigma2 <- pow(p.sigma, 2)
}
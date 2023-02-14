TransData <- function(x) {
  
  x_hat <- log(x+1)
  
  return(x_hat)
  
}

DeTransData <- function(x_hat) {
  
  x <- exp(x_hat)-1
  
  return(x)
  
}


GetMinus2LL <- function(est,dataset,C)  {
  
  est["sd"]   <- est[1]
  est["IC50"] <- est[2]


  
  sigmoidal <- function(S) 100/(1+10^(S-log(IC50)))
  
  # assign est values to constants
  IC50 <- est["IC50"]


    
  # perform model fit for measured concentrations
  fit <- sigmoidal(S=C)

    # calculate SSR
  SSR <- sum((dataset$Viability - fit)^2)
    
  #--- DERIVE -2 log(L) WITH a^2 AS VARIANCE FOR IID (INDEPENDENT IDENTICALLY DISTRUBUTED) NORMAL DISTRUBUTED ERROR
  
  # constant
  const <- length(dataset$Viability) * log(2*pi*est["sd"]^2)
  
  # minus 2 log likelihood
  minus2LL <- 1/est["sd"]^2 * SSR + const  
  
  # rv <- sample(1:100,1,replace=TRUE)
  # if(rv>90){
  #  print(minus2LL)
  #  print(est)
  # }
  
  if (is.na(minus2LL)) {
    minus2LL <- Inf
    
    print("Minus2LL returned NaN for est:")
    print(est)
    print(SSR)

    break
  }
  
  return(minus2LL)
}

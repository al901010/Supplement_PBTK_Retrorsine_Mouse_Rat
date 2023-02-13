TransData <- function(x) {
  
  x_hat <- log(x+0.00000001)
  
  return(x_hat)
  
}

DeTransData <- function(x_hat) {
  
  x <- exp(x_hat)-0.00000001
  
  return(x)
  
}



GetPred3 <- function(est,dataset,t) {
  
  
  est["sd_w"] <- est[1]
  est["k_w"]  <- est[2]
  est["sd_c"] <- est[3]
  est["k_c"]  <- est[4]
  
  
  # Create event table
  ev <- eventTable(amount.units ="percent", time.units="hours")
  
  # Add time span to event table
  ev$add.sampling(t)
  
  # Define initial conditions
  inits <- c(S_w=1,P_w=0,
             S_c=1,P_c=0)
  
  #solve ODE system
  out <- solve(mMonoexponential,est,ev,inits)
  
  #transform predicitons
  # create new rows with predicitons
  
  # warm
  new_rows_w <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_w)              <- colnames(dataset)
  new_rows_w$Condition              <- "warm"
  new_rows_w$Time                   <- t
  new_rows_w$Concentration          <- "0.7"
  new_rows_w$DepletionRate1         <- est["k_w"]
  new_rows_w$ModelSD                <- est["sd_w"]
  new_rows_w$Fit                    <- TransData(out$S_w) +  rnorm(n=length(t),mean=0,sd=est["sd_w"])
  
  # cold
  new_rows_c <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_c)              <- colnames(dataset)
  new_rows_c$Condition              <- "cold"
  new_rows_c$Time                   <- t
  new_rows_c$Concentration          <- "0.7"
  new_rows_c$DepletionRate1         <- est["k_c"]
  new_rows_c$ModelSD                <- est["sd_c"]
  new_rows_c$Fit                    <- TransData(out$S_c) +  rnorm(n=length(t),mean=0,sd=est["sd_c"])
  
  
  new_rows <- rbind(new_rows_w,new_rows_c)
  
  return(new_rows)
}

GetPred2 <- function(est,dataset,t) {
  
  est["sd_w"] <- est[1]
  est["k_w"]  <- est[2]
  est["sd_c"] <- est[3]
  est["k_c"]  <- est[4]
  
  
  # Create event table
  ev <- eventTable(amount.units ="percent", time.units="hours")
  
  # Add time span to event table
  ev$add.sampling(t)
  
  # Define initial conditions
  inits <- c(S_w=1,P_w=0,
             S_c=1,P_c=0)
  
  #solve ODE system
  out <- solve(mMonoexponential,est,ev,inits)
  
  #transform predicitons
  # create new rows with predicitons
  
  # warm
  new_rows_w <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_w)              <- colnames(dataset)
  new_rows_w$Condition              <- "warm"
  new_rows_w$Time                   <- t
  new_rows_w$Concentration          <- "0.7"
  new_rows_w$DepletionRate1         <- est["k_w"]
  new_rows_w$ModelSD                <- est["sd_w"]
  new_rows_w$Fit                    <- TransData(out$S_w)
  
  # cold
  new_rows_c <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_c)              <- colnames(dataset)
  new_rows_c$Condition              <- "cold"
  new_rows_c$Time                   <- t
  new_rows_c$Concentration          <- "0.7"
  new_rows_c$DepletionRate1         <- est["k_c"]
  new_rows_c$ModelSD                <- est["sd_c"]
  new_rows_c$Fit                    <- TransData(out$S_c)
  
 
  new_rows <- rbind(new_rows_w,new_rows_c)
  
  return(new_rows)
}


GetPred <- function(est,dataset,t) {
  
  est["sd_w"] <- est[1]
  est["k_w"]  <- est[2]
  est["sd_c"] <- est[3]
  est["k_c"]  <- est[4]
  
 
  # Create event table
  ev <- eventTable(amount.units ="percent", time.units="hours")
  
  # Add time span to event table
  ev$add.sampling(t)
  
  
  # Define initial conditions
  inits <- c(S_w=1,P_w=0,
             S_c=1,P_c=0)
  
  #solve ODE system
  out <- solve(mMonoexponential,est,ev,inits)
  
  #transform predicitons
  #transfer predictions for substrate into dataset
  dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0)]  <- rep(TransData(out$S_w[which(out$time==0)]),each=length(dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0)]))
  dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0.0417)]  <- rep(TransData(out$S_w[which(out$time==0.0417)]),each=length(dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0.0417)]))
  dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0.0833)]  <- rep(TransData(out$S_w[which(out$time==0.0833)]),each=length(dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0.0833)]))
  dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0.1670)]  <- rep(TransData(out$S_w[which(out$time==0.1670)]),each=length(dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0.1670)]))
  dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0.5000)]  <- rep(TransData(out$S_w[which(out$time==0.5000)]),each=length(dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0.5000)]))
  dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==1)]  <- rep(TransData(out$S_w[which(out$time==1)]),each=length(dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==1)]))

  dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==0)]  <- rep(TransData(out$S_c[which(out$time==0)]),each=length(dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==0)]))
  dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==0.0417)]  <- rep(TransData(out$S_c[which(out$time==0.0417)]),each=length(dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==0.0417)]))
  dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==0.0833)]  <- rep(TransData(out$S_c[which(out$time==0.0833)]),each=length(dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==0.0833)]))
  dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==0.1670)]  <- rep(TransData(out$S_c[which(out$time==0.1670)]),each=length(dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==0.1670)]))
  dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==0.5000)]  <- rep(TransData(out$S_c[which(out$time==0.5000)]),each=length(dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==0.5000)]))
  dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==1)]  <- rep(TransData(out$S_c[which(out$time==1)]),each=length(dataset$Fit[which(dataset$Condition=="cold" & dataset$Time==1)]))

  #dataset$Fit[which(dataset$Condition=="warm" & dataset$Time==0)]  <- rep(TransData(out$S_w[which(out$time==0)]),each=6)
  
  
  # dataset$Fit[which(dataset$Condition=="warm")]  <- rep(TransData(out$S_w),each=18)
  # dataset$Fit[which(dataset$Condition=="cold")]  <- rep(TransData(out$S_c),each=18)
  
  return(dataset)
}





GetMinus2LL <- function(est,dataset,t)  {
  
  est["sd_w"] <- est[1]
  est["k_w"]  <- est[2]
  est["sd_c"] <- est[3]
  est["k_c"]  <- est[4]
  
  #--- SOLVE ODE SYSTEM ---------------------------------------------------------
  dataset_pred <- GetPred(est,dataset,t)
  
  #--- CALCULATE SSR ------------------------------------------------------------
 
   # Remove columns with NormalizedConcentration is NA
  dataset_pred_narm <- dataset_pred[complete.cases(dataset_pred[ ,"NormalizedConcentration"]),]
  
  SSR_w <-  sum((dataset_pred_narm$NormalizedConcentration[which(dataset_pred_narm$Condition=="warm")] - dataset_pred_narm$Fit[which(dataset_pred_narm$Condition == "warm")])^2)
  SSR_c <-  sum((dataset_pred_narm$NormalizedConcentration[which(dataset_pred_narm$Condition=="cold")] - dataset_pred_narm$Fit[which(dataset_pred_narm$Condition == "cold")])^2)
  
  SSR <- SSR_w+SSR_c
  
  #--- DERIVE -2 log(L) WITH a^2 AS VARIANCE FOR IID (INDEPENDENT IDENTICALLY DISTRUBUTED) NORMAL DISTRUBUTED ERROR
  
  # constant
  const_w <- length(dataset_pred$Time[which(dataset_pred$Condition=="warm")]) * log(2*pi*est["sd_w"]^2)
  const_c <- length(dataset_pred$Time[which(dataset_pred$Condition=="cold")]) * log(2*pi*est["sd_c"]^2)

  
  # minus 2 log likelihood
  minus2LL_w <- 1/est["sd_w"]^2 * SSR_w + const_w 
  minus2LL_c <- 1/est["sd_c"]^2 * SSR_c + const_c
  
  minus2LL <- minus2LL_w + minus2LL_c
  
  
  # rv <- sample(1:100,1,replace=TRUE)
  # if(rv>95){
  #  print(minus2LL)
  #  print(est)
  # }
  
  if (is.na(minus2LL)) {
    minus2LL <- Inf
    
    print("Minus2LL returned NaN for est:")
    print(est)
    print(SSR)
    print(dataset_pred)
    
    break
  }
  
  return(minus2LL)
}
















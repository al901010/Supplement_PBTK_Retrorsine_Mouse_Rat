TransData <- function(x) {
  
  x_hat <- log(x+0.00000001)
  
  return(x_hat)
  
}

DeTransData <- function(x_hat) {
  
  x <- exp(x_hat)-0.00000001
  
  return(x)
  
}




#prediction with random error
GetPred3 <- function(est,dataset,t) {
  
  est["sd"]     <- est[1]
  est["EC50"]   <- est[2] 
  est["Vmax"]   <- est[3]
  est["Km"]     <- est[4]

  
  
  # Create event table
  ev <- eventTable(amount.units = NA, time.units="min")
  
  # Add time span to event table
  ev$add.sampling(t)
  
  
  # Define initial conditions
  inits <- c(S_1=1,S_15=1,S_50=1,S_200=1,P_1=0,P_15=0,P_50=0,P_200=0)
  
  
  #solve ODE system
  out <- mModel %>% rxSolve(est,ev,inits)
  
  #transform predicitons
  # create new rows with predicitons
  # add random error N~(0, sd^2)

  # 1 
  new_rows_1 <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_1)              <- colnames(dataset)
  new_rows_1$Species                <- "mouse"
  new_rows_1$Time                   <- t
  new_rows_1$Concentration          <- "1"
  new_rows_1$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="mouse"&dataset$Concentration=="1")])
  new_rows_1$ModelSD                <- est["sd"]
  new_rows_1$InhibConc              <- est["EC50"]
  new_rows_1$MaxDepletionRate       <- est["Vmax"]
  new_rows_1$MichaelisMentenConstant<- est["Km"]
  new_rows_1$Fit                    <- TransData(out$S_1) + rnorm(n=length(t),mean=0,sd=est["sd"])
  
  # 15 
  new_rows_15 <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_15)              <- colnames(dataset)
  new_rows_15$Species                <- "mouse"
  new_rows_15$Time                   <- t
  new_rows_15$Concentration          <- "15"
  new_rows_15$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="mouse"&dataset$Concentration=="15")])
  new_rows_15$ModelSD                <- est["sd"]
  new_rows_15$InhibConc              <- est["EC50"]
  new_rows_15$MaxDepletionRate       <- est["Vmax"]
  new_rows_15$MichaelisMentenConstant<- est["Km"]
  new_rows_15$Fit                    <- TransData(out$S_15) + rnorm(n=length(t),mean=0,sd=est["sd"])
  
  #50 mouse
  new_rows_50 <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_50)              <- colnames(dataset)
  new_rows_50$Species                <- "mouse"
  new_rows_50$Time                   <- t
  new_rows_50$Concentration          <- "50"
  new_rows_50$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="mouse"&dataset$Concentration=="50")])
  new_rows_50$ModelSD                <- est["sd"]
  new_rows_50$InhibConc              <- est["EC50"]
  new_rows_50$MaxDepletionRate       <- est["Vmax"]
  new_rows_50$MichaelisMentenConstant<- est["Km"]
  new_rows_50$Fit                    <- TransData(out$S_50) + rnorm(n=length(t),mean=0,sd=est["sd"])
  
  
  # 200 mouse
  new_rows_200 <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_200)              <- colnames(dataset)
  new_rows_200$Species                <- "mouse"
  new_rows_200$Time                   <- t
  new_rows_200$Concentration          <- "200"
  new_rows_200$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="mouse"&dataset$Concentration=="200")])
  new_rows_200$ModelSD                <- est["sd"]
  new_rows_200$InhibConc              <- est["EC50"]
  new_rows_200$MaxDepletionRate       <- est["Vmax"]
  new_rows_200$MichaelisMentenConstant<- est["Km"]  
  new_rows_200$Fit                    <- TransData(out$S_200) + rnorm(n=length(t),mean=0,sd=est["sd"])
  
 
  
  
  new_rows <- rbind(new_rows_1,new_rows_15,new_rows_50,new_rows_200)
  
  return(new_rows)
}



GetPred2 <- function(est,dataset,t) {
  
  
  est["sd"]     <- est[1]
  est["EC50"]   <- est[2] 
  est["Vmax"]   <- est[3]
  est["Km"]     <- est[4]
  

  # Create event table
  ev <- eventTable(amount.units = NA, time.units="min")
  
  # Add time span to event table
  ev$add.sampling(t)
  
 
  # Define initial conditions
  inits <- c(S_1=1,S_15=1,S_50=1,S_200=1,P_1=0,P_15=0,P_50=0,P_200=0)
  

  #solve ODE system
  out <- mModel %>% rxSolve(est,ev,inits)
  
  #transform predicitons
  # create new rows with predicitons
 
  # 1 mouse
  new_rows_1 <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_1)              <- colnames(dataset)
  new_rows_1$Species                <- "mouse"
  new_rows_1$Time                   <- t
  new_rows_1$Concentration          <- "1"
  new_rows_1$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="mouse"&dataset$Concentration=="1")])
  new_rows_1$ModelSD                <- est["sd"]
  new_rows_1$InhibConc              <- est["EC50"]
  new_rows_1$MaxDepletionRate       <- est["Vmax"]
  new_rows_1$MichaelisMentenConstant<- est["Km"]
  new_rows_1$Fit                    <- TransData(out$S_1) 
  
  # 15 
  new_rows_15 <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_15)              <- colnames(dataset)
  new_rows_15$Species                <- "mouse"
  new_rows_15$Time                   <- t
  new_rows_15$Concentration          <- "15"
  new_rows_15$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="mouse"&dataset$Concentration=="15")])
  new_rows_15$ModelSD                <- est["sd"]
  new_rows_15$InhibConc              <- est["EC50"]
  new_rows_15$MaxDepletionRate       <- est["Vmax"]
  new_rows_15$MichaelisMentenConstant<- est["Km"]
  new_rows_15$Fit                    <- TransData(out$S_15) 
  
  #50 mouse
  new_rows_50 <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_50)              <- colnames(dataset)
  new_rows_50$Species                <- "mouse"
  new_rows_50$Time                   <- t
  new_rows_50$Concentration          <- "50"
  new_rows_50$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="mouse"&dataset$Concentration=="50")])
  new_rows_50$ModelSD                <- est["sd"]
  new_rows_50$InhibConc              <- est["EC50"]
  new_rows_50$MaxDepletionRate       <- est["Vmax"]
  new_rows_50$MichaelisMentenConstant<- est["Km"]
  new_rows_50$Fit                    <- TransData(out$S_50)
  
  
  # 200 mouse
  new_rows_200 <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
  colnames(new_rows_200)              <- colnames(dataset)
  new_rows_200$Species                <- "mouse"
  new_rows_200$Time                   <- t
  new_rows_200$Concentration          <- "200"
  new_rows_200$ActualConcentration    <- unique(dataset$ActualConcentration[which(dataset$Species=="mouse"&dataset$Concentration=="200")])
  new_rows_200$ModelSD                <- est["sd"]
  new_rows_200$InhibConc              <- est["EC50"]
  new_rows_200$MaxDepletionRate       <- est["Vmax"]
  new_rows_200$MichaelisMentenConstant<- est["Km"]  
  new_rows_200$Fit                    <- TransData(out$S_200) 

  new_rows <- rbind(new_rows_1,new_rows_15,new_rows_50,new_rows_200)
  
  return(new_rows)
}






GetPred <- function(est,dataset,t) {

  est["sd"]     <- est[1]
  est["EC50"]   <- est[2] 
  est["Vmax"]   <- est[3]
  est["Km"]     <- est[4]
  
  # Create event table
  ev <- eventTable(amount.units = NA, time.units="min")
  
  # Add time span to event table
  ev$add.sampling(t)
  
  
  # Define initial conditions
  inits <- c(S_1=1,S_15=1,S_50=1,S_200=1,P_1=0,P_15=0,P_50=0,P_200=0)
  
  
  #solve ODE system
  out <- mModel %>% rxSolve(est,ev,inits)
  
  #transform predicitons
  #transfer predictions for substrate into dataset
  dataset$Fit[which(dataset$Concentration==1 & dataset$Species=="mouse")]  <- rep(TransData(out$S_1),each=2)
  dataset$Fit[which(dataset$Concentration==15 & dataset$Species=="mouse")] <- rep(TransData(out$S_15),each=2)
  dataset$Fit[which(dataset$Concentration==50 & dataset$Species=="mouse")] <- rep(TransData(out$S_50),each=2)
  dataset$Fit[which(dataset$Concentration==200 & dataset$Species=="mouse")]<- rep(TransData(out$S_200),each=2)

  return(dataset)
}





GetMinus2LL <- function(est,dataset,t)  {
  
  est["sd"]     <- est[1]
  est["EC50"]   <- est[2] 
  est["Vmax"]   <- est[3]
  est["Km"]     <- est[4]
  
  #--- SOLVE ODE SYSTEM ---------------------------------------------------------
  dataset_pred <- GetPred(est,dataset,t)
  
  #--- CALCULATE SSR ------------------------------------------------------------
  SSR_m <-  sum((dataset_pred$NormalizedConcentration[which(dataset_pred$Species=="mouse")] - dataset_pred$Fit[which(dataset_pred$Species == "mouse")])^2)
 
  SSR <- SSR_m

  #--- DERIVE -2 log(L) WITH a^2 AS VARIANCE FOR IID (INDEPENDENT IDENTICALLY DISTRUBUTED) NORMAL DISTRUBUTED ERROR
  
  # constant
  const_m <- length(dataset_pred$Time[which(dataset_pred$Species=="mouse")]) * log(2*pi*est["sd"]^2)

  # minus 2 log likelihood
  minus2LL_m <- 1/est["sd"]^2 * SSR_m + const_m

  minus2LL <- minus2LL_m 
  
  
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












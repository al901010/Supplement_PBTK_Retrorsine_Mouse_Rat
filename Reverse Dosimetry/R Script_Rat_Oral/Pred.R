
###############################################################################
# This function performs model predictions using the RxODE package and 
# includes model predictions into the dataset generated in MAIN.R
# Authors: Anja Lehmann, Christoph Hethey                                                       
# Date: 2023-02-13                                            
###############################################################################


###############################################################################
# PREDICTION FOR ESTIMATION INCLUDES ONLY OBSERVED TIMEPOINTS
###############################################################################

PredRxODE <- function(speciestype,
                      t,
                      runtable) {


  #--- SOLVE ODE SYSTEM FOR RUNTABLE AND PROCESS SOLUTION -----------------------
  
  # Initialise empty lists of length = nrow(runtable)
  results_ev   <- vector("list",nrow(runtable))
  results_wide <- vector("list",nrow(runtable))
  results_long <- vector("list",nrow(runtable))
  results_sum  <- vector("list",nrow(runtable))
  

  # FOR LOOP THROUGH ROWS OF RUNTABLE
  for(row in 1:nrow(runtable)) {
    
    # Speciestype related data
    theta       <- GetParameterVector(drug="Retrorsine",speciestype,
                                      method_partcoeff="RodgersAndRowland")

    # Create event table
    ev <- eventTable(amount.units = "mol", time.units="h")
    
    # Add time span to event table
    ev$add.sampling(t)
    
    # Dosing value [nmol] is converted to [mol] for simulation
    dosing_value <- runtable[row,"DOSE"]/10^9*get(paste("theta",
                                                        speciestype,sep="_"))["bw"]
    
    # Add dosing amount [mol] to event table according to administration route
    adm <- runtable[row,"ADM"]
    switch(adm,
           "ip" = {
             ev$add.dosing(
               dose=dosing_value,
               nbr.doses=1,
               dosing.to=2)
           },    
           "iv" = {
             ev$add.dosing(
               dose=dosing_value,
               nbr.doses=1,
               dosing.to=3)
           },
           "po" = {
             ev$add.dosing(
               dose=dosing_value,
               nbr.doses=1,
               dosing.to=1)
           })
    
    # Define initial conditions
    inits <- c(Alum=0,Aper=0,Aven=0,Aart=0,Aadi=0,Abon=0,Abra=0,Agut=0,Ahea=0,
               Akid=0,Alivvi=0,Alivc=0,Alun=0,Amus=0,Aski=0,Aspl=0,DHRGSH=0,
               DHRPROT=0,DHRDNA=0,Auri=0,Abil=0,bag=0)
    
    
    # Solve ODE system 
    wide <- mRetrorsine %>% rxSolve(theta,ev,inits,method="lsoda",atol=1e-12,
                                    rtol=1e-10) 
    
    # Check balance
    sum <- rowSums(wide) - parse_number(as.character(wide$time)) 
    
    # Convert units [mol] in [nmol]
    wide[,2: ncol(wide)] <-  wide[,2: ncol(wide)]*10^9
    
    # Determine Aliv and Apla
    wide$Aliv <- wide$Alivc + wide$Alivvi
    wide$Apla <- (wide$Aven + wide$Aart)*(1-get(paste("theta",speciestype,sep="_"))["hct"])/get(paste("theta",speciestype,sep="_"))["BP"]
    
    # Convert wide to long with tidyr package with new key column "compartment" 
    # and new value column "value"
    long <- gather(wide, compartment, value, 2:ncol(wide), factor_key=TRUE)
    
    # Save results in list
    results_ev[[row]]   <- ev
    results_sum[[row]]  <- sum
    results_wide[[row]] <- wide
    results_long[[row]] <- long
    
    # Identify Cmax ([µM]) of liv (Amax[nmol],V[L])
    Cmax_liv <- max(wide$Alivvi)/get(paste("theta",speciestype,sep="_"))["Vlivvi"]/1000
    Cmax_u_liv <- Cmax_liv*get(paste("theta",speciestype,sep="_"))["fuliv"]
    Cmax_u_liv <- signif(Cmax_u_liv, digits = 3)
    runtable[row,]$CMAXULIV <- Cmax_u_liv 
    
  }
  
  # # View dosing in event tables
  # print("Dosing")
  # print(results_ev[[1]]$get.dosing())
  # print(results_ev[[2]]$get.dosing())
  # print(results_ev[[3]]$get.dosing())
  # print(results_ev[[4]]$get.dosing())
  # # print(results_ev[[5]]$get.dosing())
  # 
  # # View head of balance check and results in wide format
  # print("Balance check")
  # print(head(results_sum))
  # print("Predictions wide-format")
  # print(head(results_wide))
  
  # Assign variable in the global environment as it is needed for plotting
  results_long <<- results_long 

  
  return(runtable)
  
}



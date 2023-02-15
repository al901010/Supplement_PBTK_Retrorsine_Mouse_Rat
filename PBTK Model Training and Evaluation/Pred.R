###############################################################################
# This function performs model predictions using the RxODE package and 
# includes model predictions into the dataset generated in MAIN.R
# Authors: Anja Lehmann, Christoph Hethey                                                       
# Date: 2023-02-13                                            
###############################################################################

###############################################################################
# PREDICTION FOR ESTIMATION INCLUDES ONLY OBSERVED TIMEPOINTS
###############################################################################

PredRxODE <- function(est,
                      t,
                      dataset, 
                      runtable, 
                      add_timepoints,
                      add_error) {

  est["a_RTS_m_pla"]      <- est[1]
  est["a_RTS_r_bil"]      <- est[2]
  est["a_RTS_m_uri"]      <- est[3]
  est["a_RTS_r_uri"]      <- est[4]
  est["a_DHRGSH_m_liv"]   <- est[5]
  est["a_DHRPROT_m_liv"]  <- est[6]
  est["a_DHRDNA_m_liv"]   <- est[7]
  est["lambda1DHRGSH"]    <- est[8]
  est["lambda1DHRPROT"]   <- est[9]
  est["lambda1DHRDNA"]    <- est[10]
  est["lambda2DHRDNA"]    <- est[11]
  est["kDHRDNA"]          <- est[12]
  est["fDHRGSH"]          <- est[13] 
  est["fDHRPROT"]         <- est[14] 
  est["fDHRDNA"]          <- est[15] 
  est["CLcanefr"]         <- est[16]
  est["ka"]               <- est[17]
  
  
  

  #--- SOLVE ODE SYSTEM FOR RUNTABLE AND PROCESS SOLUTION -----------------------
  

  # Initialise empty lists of length = nrow(runtable)
  results_ev   <- vector("list",nrow(runtable))
  results_wide <- vector("list",nrow(runtable))
  results_long <- vector("list",nrow(runtable))
  results_sum  <- vector("list",nrow(runtable))
  

  
  # FOR LOOP THROUGH ROWS OF RUNTABLE
  for(row in 1:nrow(runtable)) {
    
    # Speciestype related data
    speciestype <- runtable[row,"SPECIESTYPE"]
    theta       <- GetParameterVector(drug="Retrorsine",speciestype,method_partcoeff="RodgersAndRowland")
    theta       <- c(est,theta) # Combine estimates and theta
    
    if(speciestype=="mouse9weeks"){
      #theta["CLuri"]     <- theta["CLurim"]
      theta["CLcanef"]   <- 0 

      
      
    } 
    
    if(speciestype=="rat10weeks"){
      #theta["CLuri"]     <- theta["CLurir"]
      theta["CLcanef"]   <- theta["CLcanefr"]

      
    } 
    
    # Create event table
    ev <- eventTable(amount.units = "mol", time.units="h")
    
    # Add time span to event table
    ev$add.sampling(t)
    
    # Dosing value [nmol] is converted to [mol] for simulation
    dosing_value <- runtable[row,"DOSE"]/10^9*get(paste("theta",speciestype,sep="_"))["bw"]
    
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
               Akid=0,Alivvi=0, Alivc=0,Alun=0,Amus=0,Aski=0,Aspl=0,DHRGSH=0,
               DHRPROT=0,DHRDNA=0,Auri=0,Abil=0,bag=0)
    
    # Solve ODE system 
    wide <- mRetrorsine %>% rxSolve(theta,ev,inits,method="lsoda",atol=1e-12,
                                    rtol=1e-10) 
    
    # Check balance
    sum <- rowSums(wide) - parse_number(as.character(wide$time)) 
    
    # Convert units [mol] in [nmol]
    wide[,3: ncol(wide)] <-  wide[,3: ncol(wide)]*10^9
    
    # Determine Aliv and Apla
    wide$Aliv <- wide$Alivc + wide$Alivvi
    wide$Apla <- (wide$Aven + wide$Aart)*(1-get(paste("theta",speciestype,sep="_"))["hct"])/get(paste("theta",speciestype,sep="_"))["BP"]

    #--- TRANSFORM PREDICTIONS ------------------------------------------------
    #--- OPTIONAL: ADD RANDOM ERROR -------------------------------------------
    
    if(add_error==TRUE){
      
      if(speciestype=="mouse9weeks"){
        wide$Alum   <- TransData(wide$Alum)    
        wide$Aper   <- TransData(wide$Aper)     
        wide$Aven   <- TransData(wide$Aven)   
        wide$Aart   <- TransData(wide$Aart)    
        wide$Aadi   <- TransData(wide$Aadi)    
        wide$Abon   <- TransData(wide$Abon)    
        wide$Abra   <- TransData(wide$Abra)   
        wide$Agut   <- TransData(wide$Agut)    
        wide$Ahea   <- TransData(wide$Ahea)    
        wide$Akid   <- TransData(wide$Akid)   
        wide$Alivvi <- TransData(wide$Alivvi)  
        wide$Alivc  <- TransData(wide$Alivc)   
        wide$Alun   <- TransData(wide$Alun)    
        wide$Amus   <- TransData(wide$Amus)    
        wide$Aski   <- TransData(wide$Aski)    
        wide$Aspl   <- TransData(wide$Aspl)    
        wide$DHRGSH <- TransData(wide$DHRGSH)  + rnorm(n=length(t),mean=0,sd=est["a_DHRGSH_m_liv"])
        wide$DHRPROT<- TransData(wide$DHRPROT) + rnorm(n=length(t),mean=0,sd=est["a_DHRPROT_m_liv"])
        wide$DHRDNA <- TransData(wide$DHRDNA)  + rnorm(n=length(t),mean=0,sd=est["a_DHRDNA_m_liv"])
        wide$Auri   <- TransData(wide$Auri)    + rnorm(n=length(t),mean=0,sd=est["a_RTS_m_uri"])
        wide$Abil   <- TransData(wide$Abil)    
        wide$bag    <- TransData(wide$bag )
        wide$Aliv   <- TransData(wide$Aliv)    
        wide$Apla   <- TransData(wide$Apla)    + rnorm(n=length(t),mean=0,sd=est["a_RTS_m_pla"])

      } 
      
      if(speciestype=="rat10weeks"){
        wide$Alum   <- TransData(wide$Alum)    
        wide$Aper    <- TransData(wide$Aper)     
        wide$Aven   <- TransData(wide$Aven)   
        wide$Aart   <- TransData(wide$Aart)    
        wide$Aadi   <- TransData(wide$Aadi)    
        wide$Abon   <- TransData(wide$Abon)    
        wide$Abra   <- TransData(wide$Abra)   
        wide$Agut   <- TransData(wide$Agut)    
        wide$Ahea   <- TransData(wide$Ahea)    
        wide$Akid   <- TransData(wide$Akid)   
        wide$Alivvi <- TransData(wide$Alivvi)  
        wide$Alivc  <- TransData(wide$Alivc)   
        wide$Alun   <- TransData(wide$Alun)    
        wide$Amus   <- TransData(wide$Amus)    
        wide$Aski   <- TransData(wide$Aski)    
        wide$Aspl   <- TransData(wide$Aspl)    
        wide$DHRGSH <- TransData(wide$DHRGSH)  
        wide$DHRPROT<- TransData(wide$DHRPROT)  
        wide$DHRDNA <- TransData(wide$DHRDNA)  
        wide$Auri   <- TransData(wide$Auri)    + rnorm(n=length(t),mean=0,sd=est["a_RTS_r_uri"])
        wide$Abil   <- TransData(wide$Abil)    + rnorm(n=length(t),mean=0,sd=est["a_RTS_r_bil"])
        wide$bag    <- TransData(wide$bag )
        wide$Aliv   <- TransData(wide$Aliv)    
        wide$Apla   <- TransData(wide$Apla)   

        
      } 

    }
    else{
      wide$Alum   <- TransData(wide$Alum)
      wide$Aper    <- TransData(wide$Aper)
      wide$Aven   <- TransData(wide$Aven)
      wide$Aart   <- TransData(wide$Aart)
      wide$Aadi   <- TransData(wide$Aadi)
      wide$Abon   <- TransData(wide$Abon)
      wide$Abra   <- TransData(wide$Abra)
      wide$Agut   <- TransData(wide$Agut)
      wide$Ahea   <- TransData(wide$Ahea)
      wide$Akid   <- TransData(wide$Akid)
      wide$Alivvi <- TransData(wide$Alivvi)
      wide$Alivc  <- TransData(wide$Alivc)
      wide$Alun   <- TransData(wide$Alun)
      wide$Amus   <- TransData(wide$Amus)
      wide$Aski   <- TransData(wide$Aski)
      wide$Aspl   <- TransData(wide$Aspl)
      wide$DHRGSH <- TransData(wide$DHRGSH)
      wide$DHRPROT<- TransData(wide$DHRPROT)
      wide$DHRDNA <- TransData(wide$DHRDNA)
      wide$Auri   <- TransData(wide$Auri)
      wide$Abil   <- TransData(wide$Abil)
      wide$bag    <- TransData(wide$bag )
      wide$Aliv   <- TransData(wide$Aliv)
      wide$Apla   <- TransData(wide$Apla)

      }
    
    # Convert wide to long with tidyr package with new key column "compartment" 
    # and new value column "value"
    long <- gather(wide, compartment, value, 2:ncol(wide), factor_key=TRUE)
    
    # Save results in list
    results_ev[[row]]   <- ev
    results_sum[[row]]  <- sum
    results_wide[[row]] <- wide
    results_long[[row]] <- long
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

  
  #--- INCLUDE PREDICTIONS IN DATASET ------------------------------------------
  
  # Initialize new columns 
  dataset$YPRED <- NA
  dataset$YPREDUNIT <- NA
  
  # Extract ids
  ids <- unique(dataset$ID)
  
  # Loop through ids
  for (row in ids) {
    
    # Identify runid for id
    runid <-  unique(dataset$RUNID[which(dataset$ID==row)])
    # Identify rownumbers of timepoints from results_wide[[runid]]$time that 
    # match timepoints of observed data
    is_timepoints <- match(dataset$X[which(dataset$ID == row & 
                                                        is.na(dataset$DOSE))],
                           results_wide[[runid]]$time)
    
    # Transfer predictions for specific compartment at these timepoints to 
    # dataset column YPRED
    YPRED <- rep(NA, length(is_timepoints))

    if(row=="Yang_2017_retrorsine_mouse_plasma")       {YPRED         <- results_wide[[runid]]$Apla[is_timepoints]    }
    if(row=="Chu_1991_retrorsine_mouse_urine")         {YPRED         <- results_wide[[runid]]$Auri[is_timepoints]    }
    if(row=="Yang_2017_proteinadducts_mouse_liver")    {YPRED         <- results_wide[[runid]]$DHRPROT[is_timepoints] }
    if(row=="Yang_2017_GSHconjugates_mouse_liver")     {YPRED         <- results_wide[[runid]]$DHRGSH[is_timepoints]  }
    if(row=="Zhu_2017_DNAadducts_mouse_liver")         {YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]  }
    if(row=="White_1977_retrorsine_rat_bile")          {YPRED         <- results_wide[[runid]]$Abil[is_timepoints]    }
    if(row=="Chu_1991_retrorsine_rat_urine")           {YPRED         <- results_wide[[runid]]$Auri[is_timepoints]    }
    if(row=="Yang_2018_retrorsine_mouse_plasma")       {YPRED         <- results_wide[[runid]]$Apla[is_timepoints]    }
    if(row=="Yang_2018_retrorsine_mouse_liver")        {YPRED         <- results_wide[[runid]]$Aliv[is_timepoints]    }
    if(row=="Yang_2018_GSHconjugates_mouse_liver")     {YPRED         <- results_wide[[runid]]$DHRGSH[is_timepoints]  }
    if(row=="Yang_2018_GSH_mouse_liver")               {YPRED         <- results_wide[[runid]]$GSH[is_timepoints]     }
    if(row=="Zhu_2017_DNAadducts_mouse_liver_10")      {YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]  }
    if(row=="Zhu_2017_DNAadducts_mouse_liver_20")      {YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]  }
    if(row=="Zhu_2017_DNAadducts_mouse_liver_40")      {YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]  }
    if(row=="Zhu_2017_DNAadducts_mouse_liver_60")      {YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]  }
    if(row=="Yang_2018_proteinadducts_mouse_liver")    {YPRED         <- results_wide[[runid]]$DHRPROT[is_timepoints] }
    if(row=="Li_2022_retrorsine_mouse_serum")          {YPRED         <- results_wide[[runid]]$Apla[is_timepoints]    }
    if(row=="Li_2022b_retrorsine_mouse_serum")         {YPRED         <- results_wide[[runid]]$Apla[is_timepoints]    }
    
    
    dataset$YPRED[which(dataset$ID==row & !(is.na(dataset$Y)))] <- YPRED
  }
  
  #add units
  dataset$YPREDUNIT[which(is.na(dataset$DOSE))] <- "nmol"
  
  #--- INCLUDE MORE TIMEPOINTS FOR PLOTTING OF YPRED in dataset
  
  if(add_timepoints==TRUE){
    
    #loop over ids
    for (row in ids) {
      # New rows to be inserted in dataframe
      new_rows             <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
      colnames(new_rows)   <- colnames(dataset)
      new_rows$ID          <- row
      new_rows$X           <- t 
      new_rows$YTYPE       <- unique(dataset$YTYPE[which((dataset$ID==row)&(is.na(dataset[["DOSE"]])))])
      new_rows$SPECIESTYPE <- unique(dataset$SPECIESTYPE[which((dataset$ID==row)&(is.na(dataset[["DOSE"]])))])
      new_rows$XUNIT       <- unique(dataset$XUNIT[which((dataset$ID==row)&(is.na(dataset[["DOSE"]])))])
      new_rows$YPREDUNIT   <- unique(dataset$YPREDUNIT[which((dataset$ID==row)&(is.na(dataset[["DOSE"]])))])
      
      # Identify runid for id
      runid <-  unique(dataset$RUNID[which(dataset$ID==row)])
      
      # Identify rownumbers of timepoints from results_wide[[runid]]$time that match timepoints of t
      is_timepoints <- match(t,results_wide[[runid]]$time)
      
      # Transfer predictions for specific compartment at these timepoints to new_rows columns YPRED
      if(row=="Yang_2017_retrorsine_mouse_plasma")       {new_rows$YPRED         <- results_wide[[runid]]$Apla[is_timepoints]      }
      if(row=="Chu_1991_retrorsine_mouse_urine")         {new_rows$YPRED         <- results_wide[[runid]]$Auri[is_timepoints]      }
      if(row=="Yang_2017_proteinadducts_mouse_liver")    {new_rows$YPRED         <- results_wide[[runid]]$DHRPROT[is_timepoints]   }
      if(row=="Yang_2017_GSHconjugates_mouse_liver")     {new_rows$YPRED         <- results_wide[[runid]]$DHRGSH[is_timepoints]    }
      if(row=="Zhu_2017_DNAadducts_mouse_liver")         {new_rows$YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]    }
      if(row=="White_1977_retrorsine_rat_bile")          {new_rows$YPRED         <- results_wide[[runid]]$Abil[is_timepoints]      }
      if(row=="Chu_1991_retrorsine_rat_urine")           {new_rows$YPRED         <- results_wide[[runid]]$Auri[is_timepoints]      }
      if(row=="Yang_2018_retrorsine_mouse_plasma")       {new_rows$YPRED         <- results_wide[[runid]]$Apla[is_timepoints]      }
      if(row=="Yang_2018_retrorsine_mouse_liver")        {new_rows$YPRED         <- results_wide[[runid]]$Aliv[is_timepoints]      }
      if(row=="Yang_2018_GSHconjugates_mouse_liver")     {new_rows$YPRED         <- results_wide[[runid]]$DHRGSH[is_timepoints]    }
      if(row=="Zhu_2017_DNAadducts_mouse_liver_10")      {new_rows$YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]    }
      if(row=="Zhu_2017_DNAadducts_mouse_liver_20")      {new_rows$YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]    }
      if(row=="Zhu_2017_DNAadducts_mouse_liver_40")      {new_rows$YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]    }
      if(row=="Zhu_2017_DNAadducts_mouse_liver_60")      {new_rows$YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]    }
      if(row=="Yang_2018_proteinadducts_mouse_liver")    {new_rows$YPRED         <- results_wide[[runid]]$DHRPROT[is_timepoints]   }
      if(row=="Li_2022_retrorsine_mouse_serum")          {new_rows$YPRED         <- results_wide[[runid]]$Apla[is_timepoints]   }
      if(row=="Li_2022b_retrorsine_mouse_serum")         {new_rows$YPRED         <- results_wide[[runid]]$Apla[is_timepoints]   }
      
      
      #add new rows to dataset
      dataset <- rbind(dataset, new_rows)
    }
    
    
  }

    return(dataset)
  
}




###############################################################################
# PREDICT FRACTION METABOLIZED
###############################################################################

PredRxODE2 <- function(est,
                      t,
                      dataset, 
                      runtable, 
                      add_timepoints,
                      add_error) {
  
  est["a_RTS_m_pla"]      <- est[1]
  est["a_RTS_r_bil"]      <- est[2]
  est["a_RTS_m_uri"]      <- est[3]
  est["a_RTS_r_uri"]      <- est[4]
  est["a_DHRGSH_m_liv"]   <- est[5]
  est["a_DHRPROT_m_liv"]  <- est[6]
  est["a_DHRDNA_m_liv"]   <- est[7]
  est["lambda1DHRGSH"]    <- est[8]
  est["lambda1DHRPROT"]   <- est[9]
  est["lambda1DHRDNA"]    <- est[10]
  est["lambda2DHRDNA"]    <- est[11]
  est["kDHRDNA"]          <- est[12]
  est["fDHRGSH"]          <- est[13] 
  est["fDHRPROT"]         <- est[14] 
  est["fDHRDNA"]          <- est[15] 
  est["CLcanefr"]         <- est[16]
  est["ka"]               <- est[17]
  
  #--- SOLVE ODE SYSTEM FOR RUNTABLE AND PROCESS SOLUTION -----------------------
  
  # Initialise empty lists of length = nrow(runtable)
  results_ev   <- vector("list",nrow(runtable))
  results_wide <- vector("list",nrow(runtable))
  results_long <- vector("list",nrow(runtable))
  results_sum  <- vector("list",nrow(runtable))
  
  # FOR LOOP THROUGH ROWS OF RUNTABLE
  for(row in 1:nrow(runtable)) {
    
    # Speciestype related data
    speciestype <- runtable[row,"SPECIESTYPE"]
    theta       <- GetParameterVector(drug="Retrorsine",speciestype,method_partcoeff="RodgersAndRowland")
    theta       <- c(est,theta) # Combine estimates and theta
    
    if(speciestype=="mouse9weeks"){
      theta["CLcanef"]   <- 0 
    } 
    
    if(speciestype=="rat10weeks"){
      theta["CLcanef"]   <- theta["CLcanefr"]
    } 
    
    # Create event table
    ev <- eventTable(amount.units = "mol", time.units="h")
    
    # Add time span to event table
    ev$add.sampling(t)
    
    # Dosing value [nmol] is converted to [mol] for simulation
    dosing_value <- runtable[row,"DOSE"]/10^9*get(paste("theta",speciestype,sep="_"))["bw"]
    
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
               Akid=0,Alivvi=0, Alivc=0,Alun=0,Amus=0,Aski=0,Aspl=0,DHRGSH=0,
               DHRPROT=0,DHRDNA=0,Auri=0,Abil=0,AMetgut=0,AMetlivother=0)
    
    # Solve ODE system 
    wide <- mRetrorsine2 %>% rxSolve(theta,ev,inits,method="lsoda",atol=1e-12,
                                     rtol=1e-10) 
    
    # Check balance
    sum <- rowSums(wide) - parse_number(as.character(wide$time)) 
    # Convert units [mol] in [nmol]
    wide[,3: ncol(wide)] <-  wide[,3: ncol(wide)]*10^9
    
    # Determine Aliv and Apla
    wide$Aliv <- wide$Alivc + wide$Alivvi
    wide$Apla <- (wide$Aven + wide$Aart)*(1-get(paste("theta",speciestype,sep="_"))["hct"])/get(paste("theta",speciestype,sep="_"))["BP"]
    wide$AMetliv <- wide$AMetlivother + wide$DHRGSH + wide$DHRPROT + wide$DHRDNA
    
    #--- TRANSFORM PREDICTIONS ------------------------------------------------
    #--- OPTIONAL: ADD RANDOM ERROR -------------------------------------------
    
    if(add_error==TRUE){
      
      if(speciestype=="mouse9weeks"){
        wide$Alum   <- TransData(wide$Alum)    
        wide$Aper   <- TransData(wide$Aper)     
        wide$Aven   <- TransData(wide$Aven)   
        wide$Aart   <- TransData(wide$Aart)    
        wide$Aadi   <- TransData(wide$Aadi)    
        wide$Abon   <- TransData(wide$Abon)    
        wide$Abra   <- TransData(wide$Abra)   
        wide$Agut   <- TransData(wide$Agut)    
        wide$Ahea   <- TransData(wide$Ahea)    
        wide$Akid   <- TransData(wide$Akid)   
        wide$Alivvi <- TransData(wide$Alivvi)  
        wide$Alivc  <- TransData(wide$Alivc)   
        wide$Alun   <- TransData(wide$Alun)    
        wide$Amus   <- TransData(wide$Amus)    
        wide$Aski   <- TransData(wide$Aski)    
        wide$Aspl   <- TransData(wide$Aspl)    
        wide$DHRGSH <- TransData(wide$DHRGSH)  + rnorm(n=length(t),mean=0,sd=est["a_DHRGSH_m_liv"])
        wide$DHRPROT<- TransData(wide$DHRPROT) + rnorm(n=length(t),mean=0,sd=est["a_DHRPROT_m_liv"])
        wide$DHRDNA <- TransData(wide$DHRDNA)  + rnorm(n=length(t),mean=0,sd=est["a_DHRDNA_m_liv"])
        wide$Auri   <- TransData(wide$Auri)    + rnorm(n=length(t),mean=0,sd=est["a_RTS_m_uri"])
        wide$Abil   <- TransData(wide$Abil)    
        wide$Aliv   <- TransData(wide$Aliv)    
        wide$Apla   <- TransData(wide$Apla)    + rnorm(n=length(t),mean=0,sd=est["a_RTS_m_pla"])
        wide$AMetlivother <- TransData(wide$AMetlivother)
        wide$AMetliv <- TransData(wide$AMetliv)
        wide$AMetgut <- TransData(wide$AMetgut)
      } 
      
      if(speciestype=="rat10weeks"){
        wide$Alum   <- TransData(wide$Alum)    
        wide$Aper    <- TransData(wide$Aper)     
        wide$Aven   <- TransData(wide$Aven)   
        wide$Aart   <- TransData(wide$Aart)    
        wide$Aadi   <- TransData(wide$Aadi)    
        wide$Abon   <- TransData(wide$Abon)    
        wide$Abra   <- TransData(wide$Abra)   
        wide$Agut   <- TransData(wide$Agut)    
        wide$Ahea   <- TransData(wide$Ahea)    
        wide$Akid   <- TransData(wide$Akid)   
        wide$Alivvi <- TransData(wide$Alivvi)  
        wide$Alivc  <- TransData(wide$Alivc)   
        wide$Alun   <- TransData(wide$Alun)    
        wide$Amus   <- TransData(wide$Amus)    
        wide$Aski   <- TransData(wide$Aski)    
        wide$Aspl   <- TransData(wide$Aspl)    
        wide$DHRGSH <- TransData(wide$DHRGSH)  
        wide$DHRPROT<- TransData(wide$DHRPROT)  
        wide$DHRDNA <- TransData(wide$DHRDNA)  
        wide$Auri   <- TransData(wide$Auri)    + rnorm(n=length(t),mean=0,sd=est["a_RTS_r_uri"])
        wide$Abil   <- TransData(wide$Abil)    + rnorm(n=length(t),mean=0,sd=est["a_RTS_r_bil"])
        wide$Aliv   <- TransData(wide$Aliv)    
        wide$Apla   <- TransData(wide$Apla)   
        wide$AMetlivother <- TransData(wide$AMetlivother)
        wide$AMetliv <- TransData(wide$AMetliv)
        wide$AMetgut <- TransData(wide$AMetgut)
      } 
      
    }
    else{
      wide$Alum   <- TransData(wide$Alum)
      wide$Aper    <- TransData(wide$Aper)
      wide$Aven   <- TransData(wide$Aven)
      wide$Aart   <- TransData(wide$Aart)
      wide$Aadi   <- TransData(wide$Aadi)
      wide$Abon   <- TransData(wide$Abon)
      wide$Abra   <- TransData(wide$Abra)
      wide$Agut   <- TransData(wide$Agut)
      wide$Ahea   <- TransData(wide$Ahea)
      wide$Akid   <- TransData(wide$Akid)
      wide$Alivvi <- TransData(wide$Alivvi)
      wide$Alivc  <- TransData(wide$Alivc)
      wide$Alun   <- TransData(wide$Alun)
      wide$Amus   <- TransData(wide$Amus)
      wide$Aski   <- TransData(wide$Aski)
      wide$Aspl   <- TransData(wide$Aspl)
      wide$DHRGSH <- TransData(wide$DHRGSH)
      wide$DHRPROT<- TransData(wide$DHRPROT)
      wide$DHRDNA <- TransData(wide$DHRDNA)
      wide$Auri   <- TransData(wide$Auri)
      wide$Abil   <- TransData(wide$Abil)
      wide$Aliv   <- TransData(wide$Aliv)
      wide$Apla   <- TransData(wide$Apla)
      wide$AMetlivother <- TransData(wide$AMetlivother)
      wide$AMetliv <- TransData(wide$AMetliv)
      wide$AMetgut <- TransData(wide$AMetgut)
    }
    
    # Convert wide to long with tidyr package with new key column "compartment" 
    # and new value column "value"
    long <- gather(wide, compartment, value, 2:ncol(wide), factor_key=TRUE)
    
    # Save results in list
    results_ev[[row]]   <- ev
    results_sum[[row]]  <- sum
    results_wide[[row]] <- wide
    results_long[[row]] <- long
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
  
  
  #--- INCLUDE PREDICTIONS IN DATASET -------------------------------------------
  
  # Initialize new columns 
  dataset$YPRED <- NA
  dataset$YPREDUNIT <- NA
  
  # Extract ids
  ids <- unique(dataset$ID)
  
  # Loop through ids
  for (row in ids) {
    
    # Identify runid for id
    runid <-  unique(dataset$RUNID[which(dataset$ID==row)])
    # Identify rownumbers of timepoints from results_wide[[runid]]$time that match timepoints of observed data
    is_timepoints <- match(dataset$X[which(dataset$ID == row & 
                                             is.na(dataset$DOSE))],
                           results_wide[[runid]]$time)
    
    # Transfer predictions for specific compartment at these timepoints to dataset column YPRED
    YPRED <- rep(NA, length(is_timepoints))
    
    if(row=="Yang_2017_retrorsine_mouse_plasma")       {YPRED         <- results_wide[[runid]]$Apla[is_timepoints]    }
    if(row=="Chu_1991_retrorsine_mouse_urine")         {YPRED         <- results_wide[[runid]]$Auri[is_timepoints]    }
    if(row=="Yang_2017_proteinadducts_mouse_liver")    {YPRED         <- results_wide[[runid]]$DHRPROT[is_timepoints] }
    if(row=="Yang_2017_GSHconjugates_mouse_liver")     {YPRED         <- results_wide[[runid]]$DHRGSH[is_timepoints]  }
    if(row=="Zhu_2017_DNAadducts_mouse_liver")         {YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]  }
    if(row=="White_1977_retrorsine_rat_bile")          {YPRED         <- results_wide[[runid]]$Abil[is_timepoints]    }
    if(row=="Chu_1991_retrorsine_rat_urine")           {YPRED         <- results_wide[[runid]]$Auri[is_timepoints]    }
    if(row=="Yang_2018_retrorsine_mouse_plasma")       {YPRED         <- results_wide[[runid]]$Apla[is_timepoints]    }
    if(row=="Yang_2018_retrorsine_mouse_liver")        {YPRED         <- results_wide[[runid]]$Aliv[is_timepoints]    }
    if(row=="Yang_2018_GSHconjugates_mouse_liver")     {YPRED         <- results_wide[[runid]]$DHRGSH[is_timepoints]  }
    if(row=="Zhu_2017_DNAadducts_mouse_liver_10")      {YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]  }
    if(row=="Zhu_2017_DNAadducts_mouse_liver_20")      {YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]  }
    if(row=="Zhu_2017_DNAadducts_mouse_liver_40")      {YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]  }
    if(row=="Zhu_2017_DNAadducts_mouse_liver_60")      {YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]  }
    if(row=="Yang_2018_proteinadducts_mouse_liver")    {YPRED         <- results_wide[[runid]]$DHRPROT[is_timepoints] }
    if(row=="Li_2022_retrorsine_mouse_serum")          {YPRED         <- results_wide[[runid]]$Apla[is_timepoints]    }
    if(row=="Li_2022b_retrorsine_mouse_serum")         {YPRED         <- results_wide[[runid]]$Apla[is_timepoints]    }

    dataset$YPRED[which(dataset$ID==row & !(is.na(dataset$Y)))] <- YPRED
  }
  
  #add units
  dataset$YPREDUNIT[which(is.na(dataset$DOSE))] <- "nmol"
  
  #--- INCLUDE MORE TIMEPOINTS FOR PLOTTING OF YPRED in dataset
  
  if(add_timepoints==TRUE){
    
    #loop over ids
    for (row in ids) {
      # New rows to be inserted in dataframe
      new_rows             <- data.frame(matrix(NA,nrow = NROW(t), ncol = ncol(dataset)))
      colnames(new_rows)   <- colnames(dataset)
      new_rows$ID          <- row
      new_rows$X           <- t 
      new_rows$YTYPE       <- unique(dataset$YTYPE[which((dataset$ID==row)&(is.na(dataset[["DOSE"]])))])
      new_rows$SPECIESTYPE <- unique(dataset$SPECIESTYPE[which((dataset$ID==row)&(is.na(dataset[["DOSE"]])))])
      new_rows$XUNIT       <- unique(dataset$XUNIT[which((dataset$ID==row)&(is.na(dataset[["DOSE"]])))])
      new_rows$YPREDUNIT   <- unique(dataset$YPREDUNIT[which((dataset$ID==row)&(is.na(dataset[["DOSE"]])))])
      
      # Identify runid for id
      runid <-  unique(dataset$RUNID[which(dataset$ID==row)])
      
      # Identify rownumbers of timepoints from results_wide[[runid]]$time that match timepoints of t
      is_timepoints <- match(t,results_wide[[runid]]$time)
      
      # Transfer predictions for specific compartment at these timepoints to new_rows columns YPRED
      if(row=="Yang_2017_retrorsine_mouse_plasma")       {new_rows$YPRED         <- results_wide[[runid]]$Apla[is_timepoints]      }
      if(row=="Chu_1991_retrorsine_mouse_urine")         {new_rows$YPRED         <- results_wide[[runid]]$Auri[is_timepoints]      }
      if(row=="Yang_2017_proteinadducts_mouse_liver")    {new_rows$YPRED         <- results_wide[[runid]]$DHRPROT[is_timepoints]   }
      if(row=="Yang_2017_GSHconjugates_mouse_liver")     {new_rows$YPRED         <- results_wide[[runid]]$DHRGSH[is_timepoints]    }
      if(row=="Zhu_2017_DNAadducts_mouse_liver")         {new_rows$YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]    }
      if(row=="White_1977_retrorsine_rat_bile")          {new_rows$YPRED         <- results_wide[[runid]]$Abil[is_timepoints]      }
      if(row=="Chu_1991_retrorsine_rat_urine")           {new_rows$YPRED         <- results_wide[[runid]]$Auri[is_timepoints]      }
      if(row=="Yang_2018_retrorsine_mouse_plasma")       {new_rows$YPRED         <- results_wide[[runid]]$Apla[is_timepoints]      }
      if(row=="Yang_2018_retrorsine_mouse_liver")        {new_rows$YPRED         <- results_wide[[runid]]$Aliv[is_timepoints]      }
      if(row=="Yang_2018_GSHconjugates_mouse_liver")     {new_rows$YPRED         <- results_wide[[runid]]$DHRGSH[is_timepoints]    }
      if(row=="Zhu_2017_DNAadducts_mouse_liver_10")      {new_rows$YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]    }
      if(row=="Zhu_2017_DNAadducts_mouse_liver_20")      {new_rows$YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]    }
      if(row=="Zhu_2017_DNAadducts_mouse_liver_40")      {new_rows$YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]    }
      if(row=="Zhu_2017_DNAadducts_mouse_liver_60")      {new_rows$YPRED         <- results_wide[[runid]]$DHRDNA[is_timepoints]    }
      if(row=="Yang_2018_proteinadducts_mouse_liver")    {new_rows$YPRED         <- results_wide[[runid]]$DHRPROT[is_timepoints]   }
      if(row=="Li_2022_retrorsine_mouse_serum")          {new_rows$YPRED         <- results_wide[[runid]]$Apla[is_timepoints]   }
      if(row=="Li_2022b_retrorsine_mouse_serum")         {new_rows$YPRED         <- results_wide[[runid]]$Apla[is_timepoints]   }
      
      #add new rows to dataset
      dataset <- rbind(dataset, new_rows)
    }
    
    
  }
  
  return(dataset)
  
}




###############################################################################
# This function calculates the maximum likelihood              
# Authors: Anja Lehmann, Christoph Hethey                                                       
# Date: 2023-02-13                                            
###############################################################################

GetMinus2LL <- function(est,
                        t,
                        dataset, 
                        runtable) {
  
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
  
  
  
  #--- SOLVE ODE SYSTEM ---------------------------------------------------------
  
  dataset_pred <- PredRxODE(est,
                            t,
                            dataset=dataset_retrorsine, 
                            runtable=runtable_retrorsine, 
                            add_timepoints = FALSE,
                            add_error = FALSE)
  
  #--- CALCULATE SSR ------------------------------------------------------------
  
  #Residual sum of squares (RSS) = Sum of squared residuals (SSR) = Sum of squared errors of prediction (SSE)
  
  dataset_pred$RESIDUALS        <- NA
  dataset_pred$SQUAREDRESIDUALS <- NA
  dataset_pred$RESIDUALS        <- (dataset_pred$Y - dataset_pred$YPRED)
  dataset_pred$SQUAREDRESIDUALS <- (dataset_pred$Y - dataset_pred$YPRED)^2
  SSR                                      <- sum(dataset_pred$SQUAREDRESIDUALS[which(!is.na(dataset_pred$SQUAREDRESIDUALS))])

  
  
  # SSR for compound-species RTS in mouse plasma
  SSR_RTS_m_pla <- sum(dataset_pred$SQUAREDRESIDUALS[which(
                     (!is.na(dataset_pred$SQUAREDRESIDUALS)) & 
                        ((dataset_pred$ID=="Yang_2017_retrorsine_mouse_plasma")|
                         (dataset_pred$ID=="Li_2022_retrorsine_mouse_serum"))   
                     )])
  
  # SSR for compound-species RTS in rat bile
  SSR_RTS_r_bil <- sum(dataset_pred$SQUAREDRESIDUALS[which(
    (!is.na(dataset_pred$SQUAREDRESIDUALS)) & 
      ((dataset_pred$ID=="White_1977_retrorsine_rat_bile"))   
  )])
  
  # SSR for compound-species RTS in mouse urine
  SSR_RTS_m_uri <- sum(dataset_pred$SQUAREDRESIDUALS[which(
    (!is.na(dataset_pred$SQUAREDRESIDUALS)) & 
      ((dataset_pred$ID=="Chu_1991_retrorsine_mouse_urine"))   
  )])
  
  # SSR for compound-species RTS in rat urine
  SSR_RTS_r_uri <- sum(dataset_pred$SQUAREDRESIDUALS[which(
    (!is.na(dataset_pred$SQUAREDRESIDUALS)) & 
      ((dataset_pred$ID=="Chu_1991_retrorsine_rat_urine"))   
  )])
  

  # SSR for compound-species GSHconjugates
  SSR_GSHconjugates <- sum(dataset_pred$SQUAREDRESIDUALS[which(
                          (!is.na(dataset_pred$SQUAREDRESIDUALS)) &
                          (dataset_pred$ID=="Yang_2017_GSHconjugates_mouse_liver")
                                                                          )]
                           )
  
  # SSR for compound-species proteinadducts
  SSR_proteinadducts <- sum(dataset_pred$SQUAREDRESIDUALS[which(
                           (!is.na(dataset_pred$SQUAREDRESIDUALS)) &
                           (dataset_pred$ID=="Yang_2017_proteinadducts_mouse_liver")
                                                                          )]
                            )
  
  # SSR for compound-species DNAadducts
  SSR_DNAadducts <- sum(dataset_pred$SQUAREDRESIDUALS[which(
                       (!is.na(dataset_pred$SQUAREDRESIDUALS)) &
                       (dataset_pred$ID=="Zhu_2017_DNAadducts_mouse_liver")
                                                                       )]
                            )
  #SSR_check <- SSR_RTS+SSR_GSH+SSR_GSHconjugates+SSR_proteinadducts+SSR_DNAadducts
  # SSR_check <- SSR_RTS_m_pla+SSR_RTS_r_bil+SSR_RTS_m_uri+SSR_RTS_r_uri+
  #              SSR_GSH+SSR_GSHconjugates+SSR_proteinadducts+SSR_DNAadducts
  
  
  # if(SSR!=SSR_check){
  #   print("Check SSR calculation:")
  #   print(SSR)
  #   print(SSR_check)
  # }

  #--- DERIVE -2 log(L) for each compound-species with ------------------------
  #--- a_'compound-species'^2 as variance for iid normal distributed error ----
  

  # RTS_m_pla
  minus2LL_RTS_m_pla <- 1/est["a_RTS_m_pla"]^2 * SSR_RTS_m_pla +
    length(dataset_pred$X[which((!(is.na(dataset_pred$Y)))&
                                             ((dataset_pred$ID=="Yang_2017_retrorsine_mouse_plasma")|
                                              (dataset_pred$ID=="Li_2022_retrorsine_mouse_serum")))]) *
    log(2*pi*est["a_RTS_m_pla"]^2)
  
  # RTS_r_bil
  minus2LL_RTS_r_bil <- 1/est["a_RTS_r_bil"]^2 * SSR_RTS_r_bil +
    length(dataset_pred$X[which((!(is.na(dataset_pred$Y)))&
                                             (dataset_pred$ID=="White_1977_retrorsine_rat_bile"))]) *
    log(2*pi*est["a_RTS_r_bil"]^2)
  
  
  # RTS_m_uri
  minus2LL_RTS_m_uri <- 1/est["a_RTS_m_uri"]^2 * SSR_RTS_m_uri +
    length(dataset_pred$X[which((!(is.na(dataset_pred$Y)))&
                                             (dataset_pred$ID=="Chu_1991_retrorsine_mouse_urine"))]) *
    log(2*pi*est["a_RTS_m_uri"]^2)
  
  # RTS_r_uri
  minus2LL_RTS_r_uri <- 1/est["a_RTS_r_uri"]^2 * SSR_RTS_r_uri +
    length(dataset_pred$X[which((!(is.na(dataset_pred$Y)))&
                                             (dataset_pred$ID=="Chu_1991_retrorsine_rat_urine"))]) *
    log(2*pi*est["a_RTS_r_uri"]^2)
  
   # GSHconjugates
  minus2LL_GSHconjugates <- 1/est["a_DHRGSH_m_liv"]^2 * SSR_GSHconjugates +
                            length(dataset_pred$X[which((!(is.na(dataset_pred$Y)))&
                                  (dataset_pred$ID=="Yang_2017_GSHconjugates_mouse_liver"))]) * 
                            log(2*pi*est["a_DHRGSH_m_liv"]^2)
  # proteinadducts
  minus2LL_proteinadducts <- 1/est["a_DHRPROT_m_liv"]^2 * SSR_proteinadducts +
                             length(dataset_pred$X[which((!(is.na(dataset_pred$Y)))&
                                   (dataset_pred$ID=="Yang_2017_proteinadducts_mouse_liver"))]) *
                             log(2*pi*est["a_DHRPROT_m_liv"]^2)
  
  # DNAadducts
  minus2LL_DNAadducts <- 1/(10^-4*est["a_DHRDNA_m_liv"])^2 * SSR_DNAadducts +
                         length(dataset_pred$X[which((!(is.na(dataset_pred$Y)))&
                              (dataset_pred$ID=="Zhu_2017_DNAadducts_mouse_liver"))]) *
                         log(2*pi*10^-4*est["a_DHRDNA_m_liv"]^2)
  
  
  minus2LL <- sum(minus2LL_RTS_m_pla,minus2LL_RTS_r_bil,minus2LL_RTS_m_uri,minus2LL_RTS_r_uri,
                 minus2LL_GSHconjugates,minus2LL_proteinadducts,minus2LL_DNAadducts)
  

  #--- DERIVE -2 log(L) WITH a^2 AS VARIANCE FOR IID (INDEPENDENT IDENTICALLY DISTRUBUTED) NORMAL DISTRUBUTED ERROR

 
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
    print(dataset_pred)
    
    break
  }
  
  return(minus2LL)
}



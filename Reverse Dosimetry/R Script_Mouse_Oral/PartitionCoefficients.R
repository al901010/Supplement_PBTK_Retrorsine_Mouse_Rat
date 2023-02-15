###############################################################################
# This function calculates tissue-to-plasma partition coefficients
# Authors: Anja Lehmann, Christoph Hethey
# Date: 2023-02-13
#
# The function UseMethodRodgersAndRowland() is based on a MATBLAB script 
# originally written by W. Huisinga, A. Solms, L. Fronton, S. Pilari, Modeling 
# Interindividual Variability in Physiologically Based Pharmacokinetics and 
# Its Link to Mechanistic Covariate Modeling, CPT: Pharmacometrics & Systems 
# Pharmacology (2012) 1, e5; doi:10.1038/psp.2012.3                                        
###############################################################################


# We calculated the partition coefficients for steady-state conditions. 
# We described partitioning into tissue constituents (water, proteins, neutral
# lipids, neutral phospholipids, acidic phospholipids) as purely passive 
# process without active transport. The pH-dependent ionization of the molecule
# of interest is being taken into account as well, while only the neutral 
# species can passively diffuse across membranes.
#
# Rodgers and Rowlands - Assumptions:
#  Moderate to strong bases: 
#    intra-cellular distribution in neutral lipids, neutral phospholipids and 
#    acidic phospholipids;
#    no extra-cellular distribution 
#  Weak bases and acids:         
#    intra-cellular distribution in neutral lipids and neutral phospholipids;
#    extra-cellular distribution in albumin
#  Neutrals: 
#    intra-cellular distribution in neutral lipids and neutral phospholipids;
#    extra-cellular distribution in lipoproteins
#
# Schmitt - Assumptions:
#   Moderate to strong bases, weak bases, acids, neutrals: 
#     intra-cellular distribution in neutral lipids, neutral phopholipids, 
#     acidic phospholipids and proteins;
#     no extra-cellular distribution



UseMethodRodgersAndRowland <- function(DrugData,SpeciesData) {

  
  ###############################################################################
  #DATA: Kery_up (erythrocyte-to-unbound plasma partition coefficient)          
  #UNIT: [-]                                                                    
  ###############################################################################

  Kery_up <- (DrugData$BP-(1-SpeciesData$hct))/(SpeciesData$hct*DrugData$fuP)



  ###############################################################################
  #DATA: fn (fraction neutral)                                                  
  #UNIT: [-]                                                                    
  ###############################################################################

  fn <- rep(1,length.out=14)
  names(fn) <- c("blo","pla","ery","adi","bon","bra","gut",
                 "hea","kid","liv","lun","mus","ski","spl")
  
  if(DrugData$subclass=="neutral") {
    fn <- fn
    }
  if(DrugData$subclass=="acid") {
    fn <- 1/(1+10^(SpeciesData$RodgersAndRowland$pH-DrugData$pKa))
    }
  if((DrugData$subclass=="weak base")|
     (DrugData$subclass=="moderate to strong base")) {
    fn <- 1/(1+10^(DrugData$pKa-SpeciesData$RodgersAndRowland$pH))
    }

  
  ###############################################################################
  #DATA: Knlw (neutral lipids-to-water partition coefficient)                   
  #UNIT: [-]                                                                    
  ###############################################################################

  Knlw          <- rep(DrugData$Pow,length.out = 14)
  names(Knlw)   <- c("blo","pla","ery","adi","bon","bra","gut",
                     "hea","kid","liv","lun","mus","ski","spl")
  Knlw[["adi"]] <- DrugData$Pvow


  ##############################################################################
  #DATA: Knpw (neutral phospholipids-to-water partition coefficient)           
  #UNIT: [-]                                                                   
  ##############################################################################

  Knpw <- 0.3*Knlw+0.7


  ###############################################################################
  #DATA: KA_PR (Protein Affinity Constant)                                      
  #UNIT: [-]                                                                    
  ###############################################################################

  KA_PR <- rep(0,length.out=14)
  names(KA_PR) <- c("blo","pla","ery","adi","bon","bra","gut",
                    "hea","kid","liv","lun","mus","ski","spl") 
  
  if(DrugData$subclass=="neutral") {
    KA_PR <- (1/DrugData$fuP - 1 - fn[["pla"]]*Knlw[["ery"]]*SpeciesData$RodgersAndRowland$fV_nl[["pla"]] - 
              fn[["pla"]]*Knpw[["ery"]]*SpeciesData$RodgersAndRowland$fV_np[["pla"]])*SpeciesData$RodgersAndRowland$r_lpr
    } #neutrals binding to lipoproteins
  
  if((DrugData$subclass=="acid")|
     (DrugData$subclass=="weak base")) {
    KA_PR <- (1/DrugData$fuP - 1 - fn[["pla"]]*Knlw[["ery"]]*SpeciesData$RodgersAndRowland$fV_nl[["pla"]] - 
    fn[["pla"]]*Knpw[["ery"]]*SpeciesData$RodgersAndRowland$fV_np[["pla"]])*SpeciesData$RodgersAndRowland$r_alb
    } #acids and weak bases binding to albumin
  
  if(DrugData$subclass=="moderate to strong base") {
    KA_PR
    } #moderate to strong bases: no protein binding in extra-cellular space                                                                                                            
  

  ###############################################################################
  #DATA: KA_AP (Acidic Phospholipids Affinity Constant)                         
  #UNIT: [-]                                                                    
  ###############################################################################

  KA_AP <- 0
  
  if((DrugData$subclass=="neutral")|
     (DrugData$subclass=="acid")|
     (DrugData$subclass=="weak base")) {
    KA_AP
    }                                                                                                           
  if(DrugData$subclass=="moderate to strong base") {
    KA_AP <- (Kery_up - fn[["pla"]]/fn[["ery"]]*SpeciesData$RodgersAndRowland$fV_w_ic[["ery"]] - 
              fn[["pla"]]*Knlw[["ery"]]*SpeciesData$RodgersAndRowland$fV_nl[["ery"]] - 
              fn[["pla"]]*Knpw[["ery"]]*SpeciesData$RodgersAndRowland$fV_np[["ery"]])* 
              fn[["ery"]]/((1-fn[["ery"]])*fn[["pla"]])/SpeciesData$RodgersAndRowland$fV_ap[["ery"]]
    }  #moderate to strong bases: binding to intra-cellular acidic phospholipids                                               
  
  #############################################################################
  #Calculate tissue-to-unbound plasma partition coefficient Kpu               
  #############################################################################

  Kpu <- SpeciesData$RodgersAndRowland$fV_w_ex + 
         fn[["pla"]]/fn * SpeciesData$RodgersAndRowland$fV_w_ic + 
         fn[["pla"]] * Knlw * SpeciesData$RodgersAndRowland$fV_nl + 
         fn[["pla"]] * Knpw * SpeciesData$RodgersAndRowland$fV_np + 
         KA_PR + 
         fn[["pla"]] * (1-fn) / fn * KA_AP * SpeciesData$RodgersAndRowland$fV_ap 

  Kpu[["blo"]] <- DrugData$BP/DrugData$fuP     
  Kpu[["pla"]] <- 1/DrugData$fuP
  Kpu[["ery"]] <- Kery_up


  #############################################################################
  #Plausibility-CHECK                                                         
  #############################################################################

  if (DrugData$BP<(1-SpeciesData$hct)) {
    PlausabilityCheckBP <- "Blood-to-plasma ratio smaller than (1-hct)! Kery_up is negative!"
  
  } else {PlausibilityCheckBP <- "OK"}


  if((KA_PR[["adi"]]<0)|(KA_PR[["bon"]]<0)|(KA_PR[["bra"]]<0)|
     (KA_PR[["gut"]]<0)|(KA_PR[["hea"]]<0)|(KA_PR[["kid"]]<0)|
     (KA_PR[["liv"]]<0)|(KA_PR[["lun"]]<0)|(KA_PR[["mus"]]<0)|
     (KA_PR[["ski"]]<0)|(KA_PR[["spl"]]<0)<0) {
    PlausibilityCheckKA_PR <- "KA_PR is negative!"
  
  } else {PlausibilityCheckKA_PR <- "OK"}

  if (KA_AP<0) {
    PlausibilityCheckKA_AP <- "KA_AP is negative!"
  } else {PlausibilityCheckKA_AP <- "OK"}


  out <- list(Kery_up=Kery_up, fn=fn, Knlw=Knlw, Knpw=Knpw, KA_PR=KA_PR, 
              KA_AP=KA_AP,Kpu=Kpu,PlausibilityCheckBP=PlausibilityCheckBP,
              PlausibilityCheckKA_PR=PlausibilityCheckKA_PR,
              PlausibilityCheckKA_AP=PlausibilityCheckKA_AP)
  
  return(out)
}












UseMethodSchmitt <- function(DrugData,SpeciesData) {
  

  ###############################################################################
  #DATA: fn (fraction neutral)                                                  
  #UNIT: [-]                                                                    
  ###############################################################################
  
  fn <- rep(1,length.out=14)
  names(fn) <- c("blo","pla","ery","adi","bon","bra","gut",
                 "hea","kid","liv","lun","mus","ski","spl")
  
  if(DrugData$subclass=="neutral") {
    fn <- fn
  }
  if(DrugData$subclass=="acid") {
    fn <- 1/(1+10^(SpeciesData$Schmitt$pH-DrugData$pKa))
  }
  if((DrugData$subclass=="weak base")|
     (DrugData$subclass=="moderate to strong base")) {
    fn <- 1/(1+10^(DrugData$pKa-SpeciesData$Schmitt$pH))
  }
  

  
  #############################################################################
  #DATA:   alpha (ratio between distribution of the totally charged species   
  #        and the neutral form, D(totcharged)/D(neutral))                    
  #UNIT:   []                                                                 
  #SOURCE: Schmitt 2008                                                       
  #NOTE:   Distribution is usually found to be 3-4 orders of magnitude smaller
  #        for the totaly charged species than for the neutral species        
  #############################################################################  
  
  alpha <- 10^-3  
  
  #############################################################################
  #DATA: Dow  (distribution coefficient)                                      
  #UNIT: [-]                                                                  
  #NOTE: Partitioning between octanol and water strongly depends on the       
  #      charging state of a molecule (and therefore on the pH)               
  #############################################################################
  
  Dow <- rep(NA, length.out=14)
  names(Dow) <- c("blo","pla","ery","adi","bon","bra","gut",
                  "hea","kid","liv","lun","mus","ski","spl")
  
  if(DrugData$subclass=="neutral") {
    Dow <- DrugData$Pow
  }
  if(DrugData$subclass=="acid") {
    Dow <- DrugData$Pow* (fn*(1-alpha)+alpha)
  }
  if((DrugData$subclass=="weak base")|
     (DrugData$subclass=="moderate to strong base")) {
    Dow <- DrugData$Pow* (fn*(1-alpha)+alpha)
  }
  
  #############################################################################
  #DATA: Knlw (neutral lipids-to-water partition coefficient)                 
  #UNIT: [-]                                                                  
  #############################################################################
  
  Knlw <- Dow
  
  #############################################################################
  #DATA: Knpw (neutral phospholipids-to-water partition coefficient)          
  #UNIT: [-]                                                                  
  #############################################################################
  
  Knpw <- rep(DrugData$Pow, length.out=14)
  names(Knpw) <- c("blo","pla","ery","adi","bon","bra","gut",
                   "hea","kid","liv","lun","mus","ski","spl")
  
  
  #############################################################################
  #DATA: Kapw (acidic phospholipids-to-water partition coefficient)           
  #UNIT: [-]                                                                  
  #############################################################################
  
  Kapw <- rep(NA, length.out=14)
  names(Kapw) <- c("blo","pla","ery","adi","bon","bra","gut",
                   "hea","kid","liv","lun","mus","ski","spl")
  
  
  if(DrugData$subclass=="neutral") {
    Kapw <- Knpw
  }
  if(DrugData$subclass=="acid") {
    Kapw <- Knpw*(fn + 0.05*(1-fn))
  }
  if((DrugData$subclass=="weak base")|
     (DrugData$subclass=="moderate to strong base")) {
    Kapw <- Knpw*(fn + 20*(1-fn))
  }
  
  
  #############################################################################
  #DATA: Kprot (protein-to-water partition coefficient)                       
  #UNIT: [-]                                                                  
  #############################################################################
  
  Kprot <- (0.81+0.11*Knpw)/24.92*5
  
  #############################################################################
  #DATA:        Kip (partition between interstitial space and plasma)                
  #UNIT:        []                                                            
  #ASSUMPTIONS: 1) The ratio of interstitial fractional protein volume-to-    
  #                plasma fractional protein volume has typically             
  #                a value of 0.0252/0.068=0.37 (Sloop et al. 1987)           
  #             2) The water fraction in interstitial space is identical to   
  #                the water fraction in plasma, which has a value of 0.935   
  #                (Schmitt 2008)                                           
  #############################################################################  
  
   Kip <- (0.935 + 0.0252/0.068*(1/DrugData$fuP-0.935))*DrugData$fuP
  

  #############################################################################
  #DATA:        fuInt (unbound fraction in interstitial space)                
  #UNIT:        []                                                            
  #############################################################################  
  
   fuInt <- DrugData$fuP/Kip

  
  #############################################################################
  #DATA:    Kcp (partition between cellular space and plasma)                        
  #UNIT:    []                                                                
  #############################################################################  

  if(DrugData$subclass=="neutral") {
    Kcp <- (SpeciesData$Schmitt$fV_w_ic +
           Knlw*SpeciesData$Schmitt$fV_nl_cel +
           Knpw*SpeciesData$Schmitt$fV_np_cel +
           Kapw*SpeciesData$Schmitt$fV_ap_cel +
           Kprot*SpeciesData$Schmitt$fV_prot_cel) *
           DrugData$fuP
    
  }
  if(DrugData$subclass=="acid") {
    
    Kcp <- (SpeciesData$Schmitt$fV_w_ic +
           Knlw*SpeciesData$Schmitt$fV_nl_cel +
           Knpw*SpeciesData$Schmitt$fV_np_cel +
           Kapw*SpeciesData$Schmitt$fV_ap_cel +
           Kprot*SpeciesData$Schmitt$fV_prot_cel) *
           (  ( (1/(1+10^(7.4-DrugData$pKa)))*(1-alpha)+alpha ) / (fn*(1-alpha)+alpha) ) *
           DrugData$fuP
    
  }
  if((DrugData$subclass=="weak base")|
     (DrugData$subclass=="moderate to strong base")) {

    Kcp <- (SpeciesData$Schmitt$fV_w_ic +
            Knlw*SpeciesData$Schmitt$fV_nl_cel +
            Knpw*SpeciesData$Schmitt$fV_np_cel +
            Kapw*SpeciesData$Schmitt$fV_ap_cel +
            Kprot*SpeciesData$Schmitt$fV_prot_cel) *
            (  ( (1/(1+10^(DrugData$pKa-7.4)))*(1-alpha)+alpha ) / (fn*(1-alpha)+alpha) ) *
            DrugData$fuP
    
  }
  
  
  
  #############################################################################
  #DATA:    fuC (unbound fraction in cellular space)                        
  #UNIT:    []                                                                
  #############################################################################  
  
  fuC <-  1/(SpeciesData$Schmitt$fV_w_ic +
               Knlw*SpeciesData$Schmitt$fV_nl_cel +
               Knpw*SpeciesData$Schmitt$fV_np_cel +
               Kapw*SpeciesData$Schmitt$fV_ap_cel +
               Kprot*SpeciesData$Schmitt$fV_prot_cel)  # Huisinga: fV_w_ic
  
  #Set upper limit of fuC to 1
  for(i in seq_along(fuC)) {

    if(fuC[i]>1 | is.na(fuC[i]))  {
      fuC[i] <- 1
      }

    }

  #############################################################################
  #DATA: Kery_up (erythrocyte-to-unbound plasma partition coefficient)                                                                          
  #UNIT: [-]                                                                                                                                  
  #############################################################################
  
  Kery_up <- (DrugData$BP-(1-SpeciesData$hct))/(SpeciesData$hct*DrugData$fuP)
  
  
  #############################################################################
  #Calculate tissue-to-unbound plasma partition coefficient Kpu               
  #############################################################################
  
  Kpu <- Kip*SpeciesData$Schmitt$fV_int + Kcp*SpeciesData$Schmitt$fV_cel

  # blood:
  Kpu[["blo"]] <- DrugData$BP/DrugData$fuP     
  
  # plasma:
  Kpu[["pla"]] <- 1/DrugData$fuP

  # erythrocytes:
  Kpu[["ery"]] <- Kcp[["ery"]]*SpeciesData$Schmitt$fV_cel[["ery"]] 
  Kpu[["ery"]] <- Kery_up #can be used if BP is known
  

  
  #############################################################################
  #Plausibility-CHECK                                                                                                                           
  #############################################################################
  
  if (DrugData$BP<(1-SpeciesData$hct)) {
    PlausabilityCheckBP <- "Blood-to-plasma ratio smaller than (1-hct)! Kery_up is negative!"
    
  }  else {PlausibilityCheckBP <- "OK"}
  
  
  if((fuC[["adi"]]<0)|
     (fuC[["bon"]]<0)|
     (fuC[["bra"]]<0)|
     (fuC[["gut"]]<0)|
     (fuC[["hea"]]<0)|
     (fuC[["kid"]]<0)|
     (fuC[["liv"]]<0)|
     (fuC[["lun"]]<0)|
     (fuC[["mus"]]<0)|
     (fuC[["ski"]]<0)|
     (fuC[["spl"]]<0)) {
    PlausibilityCheckfuC <- "fuC is negative or bigger than 1!"
  } else {PlausibilityCheckfuC <- "OK"}
  
  
  out <- list(fuC=fuC, fuInt=fuInt,fn=fn, 
              Kip=Kip, Kcp=Kcp, Kpu=Kpu,
              PlausibilityCheckBP=PlausibilityCheckBP, 
              PlausibilityCheckfuC=PlausibilityCheckfuC)
  
  return(out)
}


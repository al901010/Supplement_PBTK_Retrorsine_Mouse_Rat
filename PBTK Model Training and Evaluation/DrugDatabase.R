
###############################################################################
# This function loads drug/compound-specific data               
# Authors: Anja Lehmann, Christoph Hethey                                                       
# Date: 2023-02-13                                            
###############################################################################

#--- DRUG-SPECIFIC DATA -------------------------------------------------------

# subclass:           "neutral", "acid", "weak base", or "moderate to strong base"
# pKa:                Acid dissociation constant
# Pow and logPow:     Octanol-to-water partition coefficient 
# Pvow and logPvow:   Vegetable oil-to-water partition coefficient
# fuP:                Fraction unbound in plasma
# BP:                 Blood-to-plasma partition coefficient  
# CLuri:              Urinary clearance [L/h]
# CLactin_perkgliver: Active transport into liver cell [L/h/kg liver]
# CLactef:            Active transport out of liver cell [L/h]
# PSdiff_perkgliver:  Passive diffusion into liver cell [L/h/kg liver]
# Vmax_perkgliver:    Maximum reaction velocity [mol/h/kg liver]
# Km:                 Michaelis Menten constant [mol/L]
# kip:                Peritoneal absorption rate constant [1/h]
# Papp:               Transcellular apparent permeability A-B (Caco-2) [cm/s]

GetDrugData <- function(drug,speciestype) {
 
  # Define species from speciestype
  if((speciestype=="rat10weeks")|(speciestype=="rat70weeks")) {
    species <- "rat"
    }  
  if((speciestype=="mouse9weeks")|(speciestype=="mouse70weeks")) {
    species <- "mouse"
    } 
  if((speciestype=="humannewborn")|(speciestype=="human1year")|
     (speciestype=="human5years")|(speciestype=="human10years")|
     (speciestype=="humanm15years")|(speciestype=="humanf15years")|
     (speciestype=="humanm35years")|(speciestype=="humanf35years")) {
    species <- "human"
    }  
  
  
  # RETRORSINE ################################################################
  if (drug=="Retrorsine") {
    
    subclass <- "weak base"                 # Lewis base
    pKa      <- 6.86                        # Predicted by SPARC at 37°C
    logPow   <- -1.26                       # Experimentally determined at 25°C
    Pow      <- 10^logPow                   # Calculate Pow
    logPvow  <- 1.115*logPow-1.35           # Calculate logPvow and Pvow according to Poulin & Theil
    Pvow     <- 10^logPvow
    
    switch(species,
      
           "mouse"={
             fuP               <- 0.600          
             BP                <- 1.08          
             CLactin_perkgliver<- 31.0           
             CLactef           <- 0            
             PSdiff_perkgliver <- 11.4           
             Vmax_perkgliver   <- 1.65*10^-3     
             Km                <- 5.54*10^-5     
             fsum              <- 1              
             kip               <- 166          
             Papp              <- 5.76e-06 
           },
                     
           "rat"={
             fuP               <- 0.600          
             BP                <- 1.08           
             CLactin_perkgliver<- 53.2           
             CLactef           <- 0              
             PSdiff_perkgliver <- 23.7           
             Vmax_perkgliver   <- 5.50*10^-3     
             Km                <- 4.43*10^-5     
             fsum              <- 1              
             kip               <- 166           
             Papp              <- 5.76e-06    
           }
           
           )

    DrugData <- list(subclass=subclass, 
                     pKa=pKa, 
                     logPow=logPow, 
                     Pow=Pow, 
                     logPvow=logPvow, 
                     Pvow=Pvow, 
                     fuP=fuP, 
                     BP=BP, 
                     CLactin_perkgliver=CLactin_perkgliver,
                     CLactef=CLactef,
                     PSdiff_perkgliver=PSdiff_perkgliver,
                     Vmax_perkgliver=Vmax_perkgliver, 
                     Km=Km, 
                     fsum=fsum,
                     kip=kip,
                     Papp=Papp
    )
                           
  }
  
  
  # RETURN DRUG PARAMETERS                      
  return(DrugData)
}






###############################################################################
# This function creates the parameter vector theta              
# Authors: Anja Lehmann, Christoph Hethey                                                       
# Date: 2023-02-13                                            
###############################################################################

GetParameterVector <- function (drug, 
                                speciestype, 
                                method_partcoeff) {
  

  #############################################################################
  #    SPECIES-SPECIFIC DATA                                                    
  #############################################################################
  
  SpeciesData <- GetSpeciesData(speciestype)
  
  # BODYWEIGHT ################################################################
  bw      <- SpeciesData$BW$bodyweight # body weight [kg]
  
  # Hematorcrit ###############################################################
  hct     <- SpeciesData$hct
  
  # VOLUMES ###################################################################
  Vvas    <- SpeciesData$fV$fV_vas_tot*SpeciesData$OW_and_V$V                   # Volume of tissues that is vascular space [L]
  SumVvas <- sum(Vvas["adi"],Vvas["bon"],Vvas["bra"],Vvas["gut"],Vvas["hea"], 
                 Vvas["kid"],Vvas["liv"],Vvas["lun"],Vvas["mus"],Vvas["ski"],
                 Vvas["spl"])                                                   # Sum of vascular volume across all tissues [L]
  Vblo    <- as.vector(SpeciesData$OW_and_V$V["blo"]-SumVvas)                   # Blood volume reduced by vascular volume of all tissues [L]
  
  Vven    <- as.vector(SpeciesData$fVvenart["ven"]*Vblo)   # Volume of venous blood [L]
  Vart    <- as.vector(SpeciesData$fVvenart["art"]*Vblo)   # Volume of arterial blood [L]
  
  Vpla    <- Vblo*(1-SpeciesData$hct)
  Very    <- Vblo - Vpla
  
  Vadi    <- as.vector(SpeciesData$OW_and_V$V["adi"]) # Tissue volumes V [L]
  Vbon    <- as.vector(SpeciesData$OW_and_V$V["bon"])
  Vbra    <- as.vector(SpeciesData$OW_and_V$V["bra"])
  Vgut    <- as.vector(SpeciesData$OW_and_V$V["gut"])
  Vhea    <- as.vector(SpeciesData$OW_and_V$V["hea"])
  Vkid    <- as.vector(SpeciesData$OW_and_V$V["kid"])
  Vliv    <- as.vector(SpeciesData$OW_and_V$V["liv"])
  Vlun    <- as.vector(SpeciesData$OW_and_V$V["lun"])
  Vmus    <- as.vector(SpeciesData$OW_and_V$V["mus"])
  Vski    <- as.vector(SpeciesData$OW_and_V$V["ski"])
  Vspl    <- as.vector(SpeciesData$OW_and_V$V["spl"]) 
  Vlivv   <- as.vector(SpeciesData$fV$fV_vas_tot["liv"]*SpeciesData$OW_and_V$V["liv"]) # Liver vascular volume [L]
  Vlivi   <- as.vector(SpeciesData$fV$fV_int_tot["liv"]*SpeciesData$OW_and_V$V["liv"]) # Liver interstitital volume  [L]
  Vlivvi  <- Vlivv + Vlivi                                                             # Liver volume of lumped vascular and interstitial space [L]
  Vlivc   <- as.vector(SpeciesData$fV$fV_cel_tot["liv"]*SpeciesData$OW_and_V$V["liv"]) # Liver intracellular volume [L]
  V       <- c(Vblo=Vblo,Vven=Vven, Vart=Vart, Vpla=Vpla,Very=Very,Vadi=Vadi,
               Vbon=Vbon,Vbra=Vbra,Vgut=Vgut,Vhea=Vhea,Vkid=Vkid,Vliv=Vliv,
               Vlun=Vlun,Vmus=Vmus,Vski=Vski,Vspl=Vspl,Vlivv=Vlivv,Vlivi=Vlivi,
               Vlivvi=Vlivvi,Vlivc=Vlivc)
  
  # BLOOD FLOWS ###############################################################
  Qadi    <- as.vector(SpeciesData$co_and_Q$Q["adi"]) # Tisse blood flows Q [L/min] and cardiac output co [L/min]
  Qbon    <- as.vector(SpeciesData$co_and_Q$Q["bon"])
  Qbra    <- as.vector(SpeciesData$co_and_Q$Q["bra"])
  Qgut    <- as.vector(SpeciesData$co_and_Q$Q["gut"])
  Qhea    <- as.vector(SpeciesData$co_and_Q$Q["hea"])
  Qkid    <- as.vector(SpeciesData$co_and_Q$Q["kid"])
  Qliv    <- as.vector(SpeciesData$co_and_Q$Q["liv"])
  Qmus    <- as.vector(SpeciesData$co_and_Q$Q["mus"])
  Qski    <- as.vector(SpeciesData$co_and_Q$Q["ski"])
  Qspl    <- as.vector(SpeciesData$co_and_Q$Q["spl"])
  Qc      <- sum(Qadi,Qbon,Qbra,Qhea,Qkid,Qliv,Qmus,Qski)     #update cardiac output: sum of tissue blood blows that go into venous blood
  Q       <- c(Qadi=Qadi, Qbon=Qbon, Qbra=Qbra, Qgut=Qgut, 
               Qhea=Qhea, Qkid=Qkid, Qliv=Qliv, Qmus=Qmus, 
               Qski=Qski, Qspl=Qspl, Qc=Qc)
  
  # GLOMERULAR FILTRATION RATE ################################################
  
  GFR     <- as.vector(SpeciesData$GFR) # Glomerular filtration rate [L/h]

  #############################################################################
  #    DRUG-SPECIFIC PARAMETERS                                                 
  #############################################################################
  DrugData <- GetDrugData(drug,speciestype)
  

  # Liver uptake and metabolism ###############################################
  CLactin      <- DrugData$CLactin_perkgliver*Vliv               # Active transport into liver cell [L/h]
  CLactef      <- DrugData$CLactef                               # Active transport out of liver cell [L/h]
  PSdiff       <- DrugData$PSdiff_perkgliver*Vliv                # Passive diffusion into liver cell [L/h]
  Km           <- DrugData$Km                                    # Michaelis Menten constant [M]
  Vmaxliv      <- DrugData$Vmax_perkgliver*Vliv                  # Maximum reaction velocity liver [M*L/h]

  CL           <- c(CLactin=CLactin,CLactef=CLactef,
                    PSdiff=PSdiff,Vmaxliv=Vmaxliv, Km=Km)
  
  # FRACTIONS OF METABOLITE FORMATION #########################################
  fsum         <- DrugData$fsum
  
  # PERITONEAL ABSORPTION RATE CONSTANT [1/h]
  kip <- DrugData$kip

  #############################################################################
  #    PARTITION COEFFICIENTS AND RELATED PARAMETERS                                                   
  #############################################################################
  
  RodgersAndRowland <- UseMethodRodgersAndRowland(DrugData,SpeciesData) 
  Schmitt <- UseMethodSchmitt(DrugData,SpeciesData)

  # Check in dataframe ########################################################
  KCheck         <- data.frame(RodgersAndRowland$Kpu, Schmitt$Kpu)
  KCheck$RoverS  <- KCheck$RodgersAndRowland.Kpu/KCheck$Schmitt.Kpu
  
  # SWITCH METHOD #############################################################
  switch(method_partcoeff,
         "RodgersAndRowland"={
           Kblo <- as.vector(RodgersAndRowland$Kpu["blo"])
           Kpla <- as.vector(RodgersAndRowland$Kpu["pla"])
           Kery <- as.vector(RodgersAndRowland$Kpu["ery"])
           Kadi <- as.vector(RodgersAndRowland$Kpu["adi"])
           Kbon <- as.vector(RodgersAndRowland$Kpu["bon"])
           Kbra <- as.vector(RodgersAndRowland$Kpu["bra"])
           Kgut <- as.vector(RodgersAndRowland$Kpu["gut"])
           Khea <- as.vector(RodgersAndRowland$Kpu["hea"])
           Kkid <- as.vector(RodgersAndRowland$Kpu["kid"])
           Kliv <- as.vector(RodgersAndRowland$Kpu["liv"])
           Klun <- as.vector(RodgersAndRowland$Kpu["lun"])
           Kmus <- as.vector(RodgersAndRowland$Kpu["mus"])
           Kski <- as.vector(RodgersAndRowland$Kpu["ski"])
           Kspl <- as.vector(RodgersAndRowland$Kpu["spl"])
           fnblo <- as.vector(RodgersAndRowland$fn["blo"])
           fnpla <- as.vector(RodgersAndRowland$fn["pla"])
           fnery <- as.vector(RodgersAndRowland$fn["ery"])
           fnadi <- as.vector(RodgersAndRowland$fn["adi"])
           fnbon <- as.vector(RodgersAndRowland$fn["bon"])
           fnbra <- as.vector(RodgersAndRowland$fn["bra"])
           fngut <- as.vector(RodgersAndRowland$fn["gut"])
           fnhea <- as.vector(RodgersAndRowland$fn["hea"])
           fnkid <- as.vector(RodgersAndRowland$fn["kid"])
           fnliv <- as.vector(RodgersAndRowland$fn["liv"])
           fnlun <- as.vector(RodgersAndRowland$fn["lun"])
           fnmus <- as.vector(RodgersAndRowland$fn["mus"])
           fnski <- as.vector(RodgersAndRowland$fn["ski"])
           fnspl <- as.vector(RodgersAndRowland$fn["spl"])
           },
         "Schmitt"          ={
           Kblo <- as.vector(Schmitt$Kpu["blo"])
           Kpla <- as.vector(Schmitt$Kpu["pla"])
           Kery <- as.vector(Schmitt$Kpu["ery"])
           Kadi <- as.vector(Schmitt$Kpu["adi"])
           Kbon <- as.vector(Schmitt$Kpu["bon"])
           Kbra <- as.vector(Schmitt$Kpu["bra"])
           Kgut <- as.vector(Schmitt$Kpu["gut"])
           Khea <- as.vector(Schmitt$Kpu["hea"])
           Kkid <- as.vector(Schmitt$Kpu["kid"])
           Kliv <- as.vector(Schmitt$Kpu["liv"])
           Klun <- as.vector(Schmitt$Kpu["lun"])
           Kmus <- as.vector(Schmitt$Kpu["mus"])
           Kski <- as.vector(Schmitt$Kpu["ski"])
           Kspl <- as.vector(Schmitt$Kpu["spl"])
           fnblo <- as.vector(Schmitt$fn["blo"])
           fnpla <- as.vector(Schmitt$fn["pla"])
           fnery <- as.vector(Schmitt$fn["ery"])
           fnadi <- as.vector(Schmitt$fn["adi"])
           fnbon <- as.vector(Schmitt$fn["bon"])
           fnbra <- as.vector(Schmitt$fn["bra"])
           fngut <- as.vector(Schmitt$fn["gut"])
           fnhea <- as.vector(Schmitt$fn["hea"])
           fnkid <- as.vector(Schmitt$fn["kid"])
           fnliv <- as.vector(Schmitt$fn["liv"])
           fnlun <- as.vector(Schmitt$fn["lun"])
           fnmus <- as.vector(Schmitt$fn["mus"])
           fnski <- as.vector(Schmitt$fn["ski"])
           fnspl <- as.vector(Schmitt$fn["spl"])
         }
  )
  
  
  # BP ratio [-] 
  BP           <- DrugData$BP  
  # Fraction unbound in plasma 
  fuP          <- DrugData$fuP 
  # Fraction unbound in interstitial space (assumed to be identical to fuP)  
  fuInt        <- fuP                    
  
  # Fraction unbound in liver cellular space (assumed to be identical to fuP) 
  fuliv <- fuP
 
  #############################################################################
  #DATA:   Kvasvi (liver vascular-to-lumped compartment partition coefficient)                    
  #UNIT:   []                                                                 
  #SOURCE: Schweinoch 2014                                                       
  #NOTE:   vi: lumped compartments vascular and interstitial space       
  #############################################################################  
  
  Kvasvi <- BP/fuP/((Vlivv/Vlivvi*BP/fuP)+(Vlivi/Vlivvi*1/fuInt))
  
  #############################################################################
  #DATA:   Kintuvi (liver unbound interstital-to-lumped compartment partition 
  #                 coefficient)                    
  #UNIT:   []                                                                 
  #SOURCE: Schweinoch 2014                                                       
  #NOTE:   vi: lumped compartments vascular and interstitial space       
  #############################################################################  
  
  Kintuvi <- fuInt/((Vlivv/Vlivvi*BP/fuP)+(Vlivi/Vlivvi*1/fuInt))
  
  
  
  partitioning  <- c(Kblo=Kblo,Kpla=Kpla,Kery=Kery,Kadi=Kadi,Kbon=Kbon,Kbra=Kbra,
                     Kgut=Kgut,Khea=Khea,Kkid=Kkid,Kliv=Kliv,Klun=Klun,Kmus=Kmus,
                     Kski=Kski,Kspl=Kspl,
                     Kvasvi=Kvasvi,Kintuvi=Kintuvi,
                     fnblo=fnblo,fnpla=fnpla,fnery=fnery,fnadi=fnadi,fnbon=fnbon,
                     fnbra=fnbra,fngut=fngut,fnhea=fnhea,fnkid=fnkid,fnliv=fnliv,
                     fnlun=fnlun,fnmus=fnmus,fnski=fnski,fnspl=fnspl,
                     BP=BP,fuP=fuP,fuInt=fuInt,fuliv=fuliv)
  
  
  #############################################################################
  # ORAL ABSOPRITON RELATED PARAMETERS
  #############################################################################
 
  #Fraction absorbed from gut tissue [-], according to Skolnik et al. 2010
  Fa <- (0.01+(1-0.01))/(1+(exp((-5.74-log10(DrugData$Papp))/0.39)))
  
  

  #############################################################################
  # THETA                                                                       
  #############################################################################
  
  theta <- c(bw=bw,
             hct=hct,
             V,
             Q,
             GFR=GFR,
             CL,
             fsum=fsum,
             kip=kip,
             partitioning,
             Fa=Fa)
  
  return(theta)
}


###############################################################################
# This function provides species-specific data
# Authors: Anja Lehmann, Christoph Hethey
# Date: 2023-02-13
#
# This function is based on a MATBLAB script originally written by            
# W. Huisinga, A. Solms, L. Fronton, S. Pilari, Modeling Interindividual      
# Variability in Physiologically Based Pharmacokinetics and Its Link to       
# Mechanistic Covariate Modeling, CPT: Pharmacometrics & Systems Pharmacology 
# (2012) 1, e5; doi:10.1038/psp.2012.3                                        
###############################################################################

GetSpeciesData <- function(speciestype) {
  
  GFR               <- GetGlomerularFiltrationRate(speciestype)             #vector       
  SA                <- GetIntestinalSurface(speciestype)                    #vector
  CYP0              <- GetCYPBaselineLiver(speciestype)                     #vector 
  GSH0liv           <- GetGSHBaselineLiver(speciestype)                     #vector
  kdegGSH           <- GetDegradationRateGSH(speciestype)                   #list
  fV                <- GetFractionalTissueVolumes(speciestype)              #list 
  fVvenart          <- GetFractionalVolumesVenousArterialBlood(speciestype) #vector
  dens              <- GetDensity(speciestype)                              #list
  BW                <- GetBodyWeight(speciestype)                           #matrix
  OW_and_V          <- GetOrganWeightAndVolume(speciestype)                 #list 
  co_and_Q          <- GetCardiacOutputAndBloodFlows(speciestype)           #list
  hct               <- GetHematocrit(speciestype)                           #1x1 vector
  RodgersAndRowland <- GetSpeciesDataRodgersAndRowland(speciestype)         #list 
  Schmitt           <- GetSpeciesDataSchmitt(speciestype)                   #list 
   
  out <- list(GFR=GFR,SA=SA,CYP0=CYP0, GSH0liv=GSH0liv,kdegGSH=kdegGSH, fV=fV, fVvenart =fVvenart , dens=dens, BW=BW, 
              OW_and_V=OW_and_V, co_and_Q=co_and_Q, hct=hct,RodgersAndRowland=RodgersAndRowland, Schmitt = Schmitt)
  
  return(out)
}

# Function that defines species from speciestype
GetSpecies <- function(speciestype) {
  if((speciestype=="rat10weeks")   |(speciestype=="rat70weeks"))   {species <- "rat"}  
  if((speciestype=="mouse9weeks")  |(speciestype=="mouse70weeks")) {species <- "mouse"} 
  if((speciestype=="humannewborn") |(speciestype=="human1year")   |
     (speciestype=="human5years")  |(speciestype=="human10years") |
     (speciestype=="humanm15years")|(speciestype=="humanf15years")|
     (speciestype=="humanm35years")|(speciestype=="humanf35years")){species <- "human"}    
  return(species)
}

#--- SPECIES-SPECIFIC DATA ----------------------------------------------------


###############################################################################
#DATA:        GFR (Glomerular filtration rate)                     
#UNIT:        [mL/min/g bw] converted to [L/h]                   
#SOURCE:      Benjamin et al. J Pharmacol Toxicol Methods 75: 101-110 (2015)                                            
#ASSUMPTIONS: - 
#INPUT:       speciestype; string defining species and age, e.g. "rat10weeks"
#OUTPUT:      see above: DATA
###############################################################################


GetGlomerularFiltrationRate <- function(speciestype) {
  
  GFRperBW <- matrix(NA,nrow=1,ncol=12)
  dimnames(GFRperBW) <- list("GFR",c("rat10weeks",   "rat70weeks",
                                     "mouse9weeks",  "mouse70weeks",
                                     "humannewborn", "human1year",
                                     "human5years",  "human10years",
                                     "humanm15years","humanf15years",
                                     "humanm35years","humanf35years"))
  
  #############   rat10weeks rat70weeks  mouse9weeks mouse70weeks humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years
  GFRperBW <- c(  0.8/100,   0.8/100,    0.8/100,    0.8/100,     NA,          NA,        NA,         NA,          NA,           NA,           0.15/100,     0.15/100       )
  
  #Convert GFRperBW from [mL/min/g bw] to GFR [L/h]
  #UNIT: *1000 to convert per g to per kg
  #      *0.06 to convert mL/min to L/h
  BW  <- GetBodyWeight(speciestype)
  GFR <- GFRperBW*1000*BW$BW*0.06

  if(speciestype=="rat10weeks")   {GFR <- GFR[ ,"rat10weeks"]} 
  if(speciestype=="rat70weeks")   {GFR <- GFR[ ,"rat70weeks"]} 
  if(speciestype=="mouse9weeks")  {GFR <- GFR[ ,"mouse9weeks"]} 
  if(speciestype=="mouse70weeks") {GFR <- GFR[ ,"mouse70weeks"]} 
  if(speciestype=="humannewborn") {GFR <- GFR[ ,"humannewborn"]} 
  if(speciestype=="human1year")   {GFR <- GFR[ ,"human1year"]} 
  if(speciestype=="human5years")  {GFR <- GFR[ ,"human5years"]} 
  if(speciestype=="human10years") {GFR <- GFR[ ,"human10years"]} 
  if(speciestype=="humanm15years"){GFR <- GFR[ ,"humanm15years"]} 
  if(speciestype=="humanf15years"){GFR <- GFR[ ,"humanf15years"]} 
  if(speciestype=="humanm35years"){GFR <- GFR[ ,"humanm35years"]} 
  if(speciestype=="humanf35years"){GFR <- GFR[ ,"humanf35years"]} 
  
  return(GFR)
}



###############################################################################
#DATA:        SA (Intestinal surface area considering the microvilli)                                                
#UNIT:        [m2] in [cm2]                                                            
#SOURCE:      rat:   Ferraris et al. Am J Physiol. 257:G689-97 (1989); Table 3         
#             mouse: Casteleyn et al. Lab. Anim. 44:176-183 (2010); Table 5  
#             human: ICRP, 1975. Report of the Task Group on Reference Man. 
#                    ICRP Publication 23. Pergamon Press, Oxford.
#STRAIN:      rat:   male Wood rat (Neotoma lepida), mean weight 100 g              
#             mouse: male CD-1 mouse (mean weight 35 g)             
#NOTE:        rat:   -                                 
#             mouse: sum of duodenum, jejunum, ileum                             
#             human: not specified if male or female, just given as adult                     
###############################################################################

GetIntestinalSurface <- function (speciestype) {
  
  SA <- matrix(NA,nrow=1,ncol=12)
  dimnames(SA) <- list("SA",c("rat10weeks",   "rat70weeks",
                              "mouse9weeks",  "mouse70weeks",
                              "humannewborn", "human1year",
                              "human5years",  "human10years",
                              "humanm15years","humanf15years",
                              "humanm35years","humanf35years"))
  
  #############   rat10weeks  rat70weeks  mouse9weeks mouse70weeks humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years
  SA[1, ] <- c(    4.10*10000,4.10*10000, 1.46*10000, 1.46*10000,  NA,          NA,        NA,         NA,          NA,           NA,           200*10000,    200*10000    )


  
  if(speciestype=="rat10weeks")   {SA <- SA[ ,"rat10weeks"]} 
  if(speciestype=="rat70weeks")   {SA <- SA[ ,"rat70weeks"]} 
  if(speciestype=="mouse9weeks")  {SA <- SA[ ,"mouse9weeks"]} 
  if(speciestype=="mouse70weeks") {SA <- SA[ ,"mouse70weeks"]} 
  if(speciestype=="humannewborn") {SA <- SA[ ,"humannewborn"]} 
  if(speciestype=="human1year")   {SA <- SA[ ,"human1year"]} 
  if(speciestype=="human5years")  {SA <- SA[ ,"human5years"]} 
  if(speciestype=="human10years") {SA <- SA[ ,"human10years"]} 
  if(speciestype=="humanm15years"){SA <- SA[ ,"humanm15years"]} 
  if(speciestype=="humanf15years"){SA <- SA[ ,"humanf15years"]} 
  if(speciestype=="humanm35years"){SA <- SA[ ,"humanm35years"]} 
  if(speciestype=="humanf35years"){SA <- SA[ ,"humanf35years"]} 
  
  return(SA)
}


###############################################################################
#DATA:        PRO0liv (absolute protein amount in liver) calculated from
#             rPRO0liv (relative protein amount in liver)                     
#UNIT:        [mg] converted from [mg/100 mg tissue weight]                        
#SOURCE:      mouse:   Verma and Chakraborty. Acta Pol. Pharm. 65: 3-9 (2008); Table 1                                   
#STRAIN:      mouse:   Mus musculus (30-33 g)                                 
#ASSUMPTIONS: rat data are identical to mouse data
###############################################################################

GetProteinAmountLiver <- function(speciestype) {

  OW_and_V <- GetOrganWeightAndVolume(speciestype)
  
  rPRO0liv <- matrix(NA,nrow=1,ncol=12)
  dimnames(rPRO0liv) <- list("liv", c("rat10weeks",   "rat70weeks",
                                      "mouse9weeks",  "mouse70weeks",
                                      "humannewborn", "human1year",
                                      "human5years",  "human10years",
                                      "humanm15years","humanf15years",
                                      "humanm35years","humanf35years"))
  ##################       rat10weeks rat70weeks  mouse9weeks mouse70weeks   humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years
  rPRO0liv["liv", ] <- c(  23.84,     23.84,      23.84,      23.84,         NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  
  #Convert liver weight from [kg] to [mg] by multiplication with 1000000
  #Convert relative protein amount in liver [mg/100 mg tissue] to absolute protein amount in liver[mg] 
  #by multiplication with liver weight and division by 100
  PRO0liv <- (rPRO0liv*OW_and_V$OW_experimental["liv",]*1000000)/100 
  
  return(PRO0liv)
}

###############################################################################
#DATA:        rCYP0liv (relative CYP amount in liver) converted to            
#             CYP0liv (absolute CYP amount in liver) converted to             
#             cCYP0livc (CYP concentration in liver cellular compartment)     
#UNIT:        [fmol/mg protein] converted to [mol] converted to [mol/L]         
#CYPTYPE:     mouse:   P3A11                                                  
#SOURCE:      mouse:   Moskaleva et al. PLoS ONE 10: 1-10 (2015); Table 1                               
#STRAIN:      mouse:   C57/BL6 line, 4-5 months old, male                 
###############################################################################  

GetCYPBaselineLiver <- function (speciestype) {
  
  PRO0liv <- GetProteinAmountLiver(speciestype)
  OW_and_V <- GetOrganWeightAndVolume(speciestype)
  fV       <- GetFractionalTissueVolumes(speciestype)   
  
  rCYP0liv <- matrix(NA,nrow=1,ncol=12)
  dimnames(rCYP0liv) <- list("liv", c("rat10weeks",   "rat70weeks",
                                      "mouse9weeks",  "mouse70weeks",
                                      "humannewborn", "human1year",
                                      "human5years",  "human10years",
                                      "humanm15years","humanf15years",
                                      "humanm35years","humanf35years"))
  ##################   rat10weeks rat70weeks  mouse9weeks   mouse70weeks  humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years
  rCYP0liv["liv", ] <- c(   NA,        NA,        7.97,         7.97,         NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  

  #Convert relative CYP amount in liver [fmol/mg protein] to absolute CYP amount 
  #in liver [mol] by multiplication with total liver protein amount [mg] and by
  #division by 10^15
  CYP0liv <- rCYP0liv*PRO0liv/10^15
 
  #Calculate CYP concentration in liver cellular compartment cCYP0livc[mol/L]
  #by division with volume of liver cellular compartment [L]
  cCYP0livc <- CYP0liv/OW_and_V$OW_experimental["liv",]*fV$fV_cel_tot["liv"] # note that for human and mouse fV is assumed to be identical to rat 
  
  
  if(speciestype=="rat10weeks")   {CYP0liv <- CYP0liv[ ,"rat10weeks"]} 
  if(speciestype=="rat70weeks")   {CYP0liv <- CYP0liv[ ,"rat70weeks"]} 
  if(speciestype=="mouse9weeks")  {CYP0liv <- CYP0liv[ ,"mouse9weeks"]} 
  if(speciestype=="mouse70weeks") {CYP0liv <- CYP0liv[ ,"mouse70weeks"]} 
  if(speciestype=="humannewborn") {CYP0liv <- CYP0liv[ ,"humannewborn"]} 
  if(speciestype=="human1year")   {CYP0liv <- CYP0liv[ ,"human1year"]} 
  if(speciestype=="human5years")  {CYP0liv <- CYP0liv[ ,"human5years"]} 
  if(speciestype=="human10years") {CYP0liv <- CYP0liv[ ,"human10years"]} 
  if(speciestype=="humanm15years"){CYP0liv <- CYP0liv[ ,"humanm15years"]} 
  if(speciestype=="humanf15years"){CYP0liv <- CYP0liv[ ,"humanf15years"]} 
  if(speciestype=="humanm35years"){CYP0liv <- CYP0liv[ ,"humanm35years"]} 
  if(speciestype=="humanf35years"){CYP0liv <- CYP0liv[ ,"humanf35years"]} 
  
  if(speciestype=="rat10weeks")   {cCYP0livc <- cCYP0livc[ ,"rat10weeks"]} 
  if(speciestype=="rat70weeks")   {cCYP0livc <- cCYP0livc[ ,"rat70weeks"]} 
  if(speciestype=="mouse9weeks")  {cCYP0livc <- cCYP0livc[ ,"mouse9weeks"]} 
  if(speciestype=="mouse70weeks") {cCYP0livc <- cCYP0livc[ ,"mouse70weeks"]} 
  if(speciestype=="humannewborn") {cCYP0livc <- cCYP0livc[ ,"humannewborn"]} 
  if(speciestype=="human1year")   {cCYP0livc <- cCYP0livc[ ,"human1year"]} 
  if(speciestype=="human5years")  {cCYP0livc <- cCYP0livc[ ,"human5years"]} 
  if(speciestype=="human10years") {cCYP0livc <- cCYP0livc[ ,"human10years"]} 
  if(speciestype=="humanm15years"){cCYP0livc <- cCYP0livc[ ,"humanm15years"]} 
  if(speciestype=="humanf15years"){cCYP0livc <- cCYP0livc[ ,"humanf15years"]} 
  if(speciestype=="humanm35years"){cCYP0livc <- cCYP0livc[ ,"humanm35years"]} 
  if(speciestype=="humanf35years"){cCYP0livc <- cCYP0livc[ ,"humanf35years"]} 
  
  
  out <- list(CYP0liv=CYP0liv, cCYP0livc=cCYP0livc)
  
  return(out) 
}


###############################################################################
#DATA:        mouse:                                                          
#             rGSH0liv (relative glutathione amount in liver) converted to    
#              GSH0liv (absolute GSH amount in liver) and            
#             cGSH0liv (GSH concentration in liver)               
#            
#             rat:                                                            
#             cGSH0liv (GSH concentation in liver) converted to               
#             cGSH0livc (GSH concentration in liver cellular !!! compartment) 
#             human:                                                          
#             cGSH0liv (GSH concentration in liver)                           
#UNIT:        mouse:                                                          
#             [ug/mg protein] converted to [mol] and [mol/L]         
#             rat:                                                            
#             [mol/g] converted to [mol/L]                                   
#             human: [mmol/L] converted to [mol/L]                            
#SOURCE:      rat:     Griffith and Meister. Proc. Nati. Acad. Sci. 76: 5606-5610 (1979); Table 2
#             mouse:   Stohs et al. AGE. 3: 11-14 (1980); Table 1  
#             human:   Skamarauskas et al. Hepatology 59: 2321-2330 (2014); p.2328                             
#STRAIN:      rat:     rats 90 - 100 g (no additional information)            
#             mouse:   male mice of the CBF-1  strain                         
#ASSUMPTIONS: rat: rats of weight 90 - 100 g are assumed to be of age 10 weeks
#                  rat70weeks is assumed to be identical to rat10weeks        
#             human: Since there were no information provided on age of study 
#                    participants, it is assumed that cGSH0liv was determined 
#                    in adults. It is further assumed, that cGSH0liv of       
#                    children and newborn is identical to adult values.       
###############################################################################  
GetGSHBaselineLiver <- function(speciestype) {
  
  PRO0liv  <- GetProteinAmountLiver(speciestype)
  OW_and_V <- GetOrganWeightAndVolume(speciestype)
  fV       <- GetFractionalTissueVolumes(speciestype)   
  
  rGSH0liv <- matrix(NA,nrow=1,ncol=12)
  dimnames(rGSH0liv) <- list("liv", c("rat10weeks",   "rat70weeks",
                                      "mouse9weeks",  "mouse70weeks",
                                      "humannewborn", "human1year",
                                      "human5years",  "human10years",
                                      "humanm15years","humanf15years",
                                      "humanm35years","humanf35years"))
  ##################   rat10weeks rat70weeks  mouse9weeks   mouse70weeks  humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years
  rGSH0liv["liv", ] <- c(   NA,        NA,       5.93,         6.77,         NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  
  #Convert relative GSH [ug/mg protein] to absolute GSH [mol]
  #by multiplication with absolute protein amount [mg]
  #by division by 1000000 [ug] in [g]
  #by division with GSH molar mass 307.32 [g/mol]
  GSH0liv  <- rGSH0liv*PRO0liv/1000000/307.32
  
  #Calculate GSH concentration in liver cellular compartment GSH0liv[mol/L]
  #by division with volume of liver cellular compartment [L]
  cGSH0liv <- GSH0liv/OW_and_V$OW_experimental["liv",]   
  
  ######################################################################
  # Rat data cGSHliv [umol/g] to [mol/kg] by division with 1000        
  ######################################################################
  cGSH0liv_rat <- cGSH0liv[,1:2]
  
  ###########        rat10weeks  rat70weeks 
  cGSH0liv_rat <- c( 4.51/1000,    4.51/1000 )
  
  cGSH0liv[,1:2] <- cGSH0liv_rat
  
  ######################################################################
  # Human data cGSHliv [mmol/L] converted to [mol/L]                   
  ######################################################################
  cGSH0liv_human <- cGSH0liv[,5:12]
  
  ###############      humannewborn  human1year  human5years  human10years  humanm15years  humanf15years  humanm35years  humanf35years
  cGSH0liv_human <- c( 2.64/1000,    2.64/1000,  2.64/1000,   2.64/1000,    2.64/1000,     2.64/1000,     2.64/1000,     2.64/1000      )
  
  cGSH0liv[,5:12] <- cGSH0liv_human
  
  #Again GSH0liv to include rat and human data
  GSH0liv <- cGSH0liv*OW_and_V$OW_experimental["liv",] 
  
  if(speciestype=="rat10weeks")   {GSH0liv <- GSH0liv[ ,"rat10weeks"]} 
  if(speciestype=="rat70weeks")   {GSH0liv <- GSH0liv[ ,"rat70weeks"]} 
  if(speciestype=="mouse9weeks")  {GSH0liv <- GSH0liv[ ,"mouse9weeks"]} 
  if(speciestype=="mouse70weeks") {GSH0liv <- GSH0liv[ ,"mouse70weeks"]} 
  if(speciestype=="humannewborn") {GSH0liv <- GSH0liv[ ,"humannewborn"]} 
  if(speciestype=="human1year")   {GSH0liv <- GSH0liv[ ,"human1year"]} 
  if(speciestype=="human5years")  {GSH0liv <- GSH0liv[ ,"human5years"]} 
  if(speciestype=="human10years") {GSH0liv <- GSH0liv[ ,"human10years"]} 
  if(speciestype=="humanm15years"){GSH0liv <- GSH0liv[ ,"humanm15years"]} 
  if(speciestype=="humanf15years"){GSH0liv <- GSH0liv[ ,"humanf15years"]} 
  if(speciestype=="humanm35years"){GSH0liv <- GSH0liv[ ,"humanm35years"]} 
  if(speciestype=="humanf35years"){GSH0liv <- GSH0liv[ ,"humanf35years"]} 
  
  return(GSH0liv)
}


###############################################################################
#DATA:        kdegGSH (1st-order degradation rate constant of GSH)            
#UNIT:        [1/h]                                                           
#SOURCE:      Griffith and Meister. Proc. Nati. Acad. Sci. 76: 5606-5610 (1979)                                                 
#STRAIN:      rats 90-100 g, mice 28-35 g (no additional information)         
#ASSUMPTIONS: 1) Rats are 10 weeks old and data for rat70weeks are identical  
#                to rat10weeks.                                               
#             2) Mice are 9 weeks old and data for mouse 70weeks are identical
#                to mouse9weeks.                                              
###############################################################################

GetDegradationRateGSH <- function(speciestype) {
  
  kdegGSH <- matrix(NA,nrow=1,ncol=12)
  dimnames(kdegGSH) <- list("liv", c("rat10weeks",   "rat70weeks",
                                     "mouse9weeks",  "mouse70weeks",
                                     "humannewborn", "human1year",
                                     "human5years",  "human10years",
                                     "humanm15years","humanf15years",
                                     "humanm35years","humanf35years"))
  ##################       rat10weeks rat70weeks  mouse9weeks mouse70weeks   humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years
  kdegGSH["liv", ] <- c(   0.596,     0.596,      0.528,      0.528,         NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )

  if(speciestype=="rat10weeks")   {kdegGSH <- kdegGSH[ ,"rat10weeks"]}
  if(speciestype=="rat70weeks")   {kdegGSH <- kdegGSH[ ,"rat70weeks"]}
  if(speciestype=="mouse9weeks")  {kdegGSH <- kdegGSH[ ,"mouse9weeks"]}
  if(speciestype=="mouse70weeks") {kdegGSH <- kdegGSH[ ,"mouse70weeks"]}
  if(speciestype=="humannewborn") {kdegGSH <- kdegGSH[ ,"humannewborn"]}
  if(speciestype=="human1year")   {kdegGSH <- kdegGSH[ ,"human1year"]}
  if(speciestype=="human5years")  {kdegGSH <- kdegGSH[ ,"human5years"]}
  if(speciestype=="human10years") {kdegGSH <- kdegGSH[ ,"human10years"]}
  if(speciestype=="humanm15years"){kdegGSH <- kdegGSH[ ,"humanm15years"]}
  if(speciestype=="humanf15years"){kdegGSH <- kdegGSH[ ,"humanf15years"]}
  if(speciestype=="humanm35years"){kdegGSH <- kdegGSH[ ,"humanm35years"]}
  if(speciestype=="humanf35years"){kdegGSH <- kdegGSH[ ,"humanf35years"]}

  return(kdegGSH)

}

GetFractionalTissueVolumes <- function (speciestype) {
  
  ###############################################################################
  #DATA:        fOW_vas and fOW_int (Fraction of experimental organ weight that 
  #             is vascular volume and interstitial space volume in non-bled    
  #             rats)                                                           
  #UNIT:        [fraction]                                                         
  #SOURCE:      rat:   Kawai et al. J Pharmacokinet Biopharm. 22: 327-364 (1994); Table B-I
  #                    based on measurements in nonbled rats, see Appendix B, p. 362                                        
  #NOTE:        The vascular volume fractions in Kawai et al. are in good        
  #             agreement with the mean residual blood data in 
  #             Brown et al. Toxicol. Ind. Health 13: 407-484 (1997); Table 30 (rat). 
  #             There, however, it is mentioned (on p.457, Section 'Blood volume data') 
  #             that the values in Table 30 are not representations of the 
  #             fraction of the total blood volume that resides in the tissue.                  
  #ASSUMPTIONS: 1) We assume that fractions with respect to total organ weight  
  #                are approximately identical to the fractions with respect to 
  #                experimental organ weight, since fractions were determined  
  #                from non-bled rats according to the footnote of table B-I,   
  #                Kawai et al.                                       
  #             2) Moreover, we assume that fractions with respect to volume    
  #                are identical to those with respect to weight.                
  #             3) mouse: data identical to rat data. This is supported by      
  #                       Brown et al., Table 30, comparing the data for mouse   
  #                       and rat, when taking the reported range of values     
  #                       (n=3) into account. Note that the mean has been       
  #                       reported, although it is very sensitive to outliers.       
  #             3) human: data identical to rat values                          
  ###############################################################################
 
  #fraction of experimental organ weight that is vascular space
  fOW_vas_experimental <- matrix(NA,nrow=14,ncol=3)
  dimnames(fOW_vas_experimental) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                                           "hea","kid","liv","lun","mus","ski","spl"),
                                         c("rat","human","mouse"))
  
  ###################################   rat     mouse   human 
  fOW_vas_experimental["blo", ] <- c(   NA,     NA,     NA    )
  fOW_vas_experimental["pla", ] <- c(   NA,     NA,     NA    )
  fOW_vas_experimental["ery", ] <- c(   NA,     NA,     NA    )
  fOW_vas_experimental["adi", ] <- c(   0.010,  0.010,  0.010 )
  fOW_vas_experimental["bon", ] <- c(   0.041,  0.041,  0.041 )
  fOW_vas_experimental["bra", ] <- c(   0.037,  0.037,  0.037 )
  fOW_vas_experimental["gut", ] <- c(   0.024,  0.024,  0.024 )
  fOW_vas_experimental["hea", ] <- c(   0.262,  0.262,  0.262 )
  fOW_vas_experimental["kid", ] <- c(   0.105,  0.105,  0.105 )
  fOW_vas_experimental["liv", ] <- c(   0.115,  0.115,  0.115 )
  fOW_vas_experimental["lun", ] <- c(   0.262,  0.262,  0.262 )
  fOW_vas_experimental["mus", ] <- c(   0.026,  0.026,  0.026 )
  fOW_vas_experimental["ski", ] <- c(   0.019,  0.019,  0.019 )
  fOW_vas_experimental["spl", ] <- c(   0.282,  0.282,  0.282 )
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {fOW_vas_experimental <- fOW_vas_experimental[ ,"rat"]} 
  if(species=="mouse") {fOW_vas_experimental <- fOW_vas_experimental[ ,"mouse"]} 
  if(species=="human") {fOW_vas_experimental <- fOW_vas_experimental[ ,"human"]} 
  
  #fraction of experimental organ weight that is interstitial space
  fOW_int_experimental <- matrix(NA,nrow=14,ncol=3)
  dimnames(fOW_int_experimental) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                                           "hea","kid","liv","lun","mus","ski","spl"),
                                         c("rat","human","mouse"))

  ###################################   rat     mouse   human 
  fOW_int_experimental["blo", ] <- c(   NA,     NA,     NA    )
  fOW_int_experimental["pla", ] <- c(   NA,     NA,     NA    )
  fOW_int_experimental["ery", ] <- c(   NA,     NA,     NA    )
  fOW_int_experimental["adi", ] <- c(   0.135,  0.135,  0.135 )
  fOW_int_experimental["bon", ] <- c(   0.100,  0.100,  0.100 )
  fOW_int_experimental["bra", ] <- c(   0.004,  0.004,  0.004 )
  fOW_int_experimental["gut", ] <- c(   0.094,  0.094,  0.094 )
  fOW_int_experimental["hea", ] <- c(   0.100,  0.100,  0.100 )
  fOW_int_experimental["kid", ] <- c(   0.200,  0.200,  0.200 )
  fOW_int_experimental["liv", ] <- c(   0.163,  0.163,  0.163 )
  fOW_int_experimental["lun", ] <- c(   0.188,  0.188,  0.188 )
  fOW_int_experimental["mus", ] <- c(   0.120,  0.120,  0.120 )
  fOW_int_experimental["ski", ] <- c(   0.302,  0.302,  0.302 )
  fOW_int_experimental["spl", ] <- c(   0.150,  0.150,  0.150 )
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {fOW_int_experimental <- fOW_int_experimental[ ,"rat"]} 
  if(species=="mouse") {fOW_int_experimental <- fOW_int_experimental[ ,"mouse"]} 
  if(species=="human") {fOW_int_experimental <- fOW_int_experimental[ ,"human"]} 
  
  # Calculate fraction of experimental organ weight that is intra-cellular space
  fOW_cel_experimental <- 1- (fOW_vas_experimental + fOW_int_experimental)
  fOW_cel_experimental
  
  #ASSUMPTION: We assume that fractions with respect to total organ weight are 
  #            approximately identical to the fractions with respect to experimental 
  #            organ weight, since fractions were determined from non-bled rats 
  #            according to the footnote of table B-I, Kawai et al. (2004).
  fOW_vas_tot <- fOW_vas_experimental
  fOW_int_tot <- fOW_int_experimental
  fOW_cel_tot <- fOW_cel_experimental
  
  #ASSUMPTION: Moreover, we assume that fractions with respect to volume are 
  #            identical to those with respect to weight.
  fV_vas_tot <- fOW_vas_tot
  fV_int_tot <- fOW_int_tot
  fV_cel_tot <- fOW_cel_tot
  
  #RESCALE
  #Determining fraction of interstitial and intra-cellular space with respect 
  #to tissue weight NOT INCLUDING regional vascular blood so that fV_int + fV_cel = 1
  fV_int <- fV_int_tot/(fV_int_tot + fV_cel_tot)
  fV_cel <- 1-fV_int
  
  fV <- list(fV_vas_tot=fV_vas_tot,fV_int_tot=fV_int_tot, fV_cel_tot=fV_cel_tot, 
             fV_int=fV_int, fV_cel=fV_cel)
  return(fV)
}


###############################################################################
#DATA:        fVvenart (fractional volume venous and arterial blood)          
#UNIT:        [fraction of blood]                                               
#ASSUMPTIONS: rat, mouse: Assume that venous and arterial blood are 2/3 and   
#                         1/3 of total blood                                  
###############################################################################

GetFractionalVolumesVenousArterialBlood <- function(speciestype) {
  
  fVvenart <- matrix(NA,nrow=2,ncol=12)
  dimnames(fVvenart) <- list(c("ven","art"), c("rat10weeks",   "rat70weeks",
                                               "mouse9weeks",  "mouse70weeks",
                                               "humannewborn", "human1year",
                                               "human5years",  "human10years",
                                               "humanm15years","humanf15years",
                                               "humanm35years","humanf35years")) 
  
  
  #######################   rat10weeks rat70weeks  mouse9weeks mouse70weeks humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years 
  fVvenart["ven", ] <- c(   2/3,       2/3,        2/3,        2/3,         NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA            )
  fVvenart["art", ] <- c(   1/3,       1/3,        1/3,        1/3,         NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA            )
  
  
  #############################################################################
  # For humanm35years: Regional vascular blood volumes and fractions of       
  #                    total blood                                            
  # UNIT:        L and fraction of total blood                                
  # SOURCE:      ICRP, 2002. Basic Anatomical and Physiological Data for Use 
  #              in Radiological Protection Reference Values. ICRP Publication 89. 
  #              Ann. ICRP 32 (3-4); Table 2.13                              
  # ASSUMPTIONS: 1) For female adults, male data for fraction of venous and   
  #                 arterial blood were adopted                                
  #              2) For children of age 1,5, 10 and 15f, adult female data    
  #                 were adopted for all fVblood entries. For children of age 
  #                 15m, corresponding adult male data were adopted. This is  
  #                 in line with the NHANES study (for age 5 and older)       
  #############################################################################
  
  ICRP_fVblo_heart_chambers              <-  9.00/100
  ICRP_fVblo_pulmonary_arteries          <-  3.00/100
  ICRP_fVblo_pulmonary_veins             <-  5.50/100
  ICRP_fVblo_pulmonary_capillaries       <-  2.00/100
  ICRP_fVblo_systemic_aortalargearteries <-  6.00/100
  ICRP_fVblo_systemic_capillaries        <-  5.00/100
  ICRP_fVblo_systemic_smallveins         <- 41.50/100
  ICRP_fVblo_systemic_smallarteries      <- 10.00/100
  ICRP_fVblo_systemic_largeveins         <- 18.00/100
  
  
  fVvenart["art",5:12] <- 0.5*ICRP_fVblo_heart_chambers               + 
                              ICRP_fVblo_pulmonary_veins              + 
                          0.5*ICRP_fVblo_pulmonary_capillaries        +
                              ICRP_fVblo_systemic_aortalargearteries  +
                          0.5*ICRP_fVblo_systemic_capillaries         +
                              ICRP_fVblo_systemic_smallarteries
  
  fVvenart["ven",5:12] <- 0.5*ICRP_fVblo_heart_chambers               + 
                              ICRP_fVblo_pulmonary_arteries           + 
                          0.5*ICRP_fVblo_pulmonary_capillaries        +
                          0.5*ICRP_fVblo_systemic_capillaries         +
                              ICRP_fVblo_systemic_smallveins          +
                              ICRP_fVblo_systemic_largeveins
  
  if(speciestype=="rat10weeks")   {fVvenart <- fVvenart[ ,"rat10weeks"]} 
  if(speciestype=="rat70weeks")   {fVvenart <- fVvenart[ ,"rat70weeks"]} 
  if(speciestype=="mouse9weeks")  {fVvenart <- fVvenart[ ,"mouse9weeks"]} 
  if(speciestype=="mouse70weeks") {fVvenart <- fVvenart[ ,"mouse70weeks"]} 
  if(speciestype=="humannewborn") {fVvenart <- fVvenart[ ,"humannewborn"]} 
  if(speciestype=="human1year")   {fVvenart <- fVvenart[ ,"human1year"]} 
  if(speciestype=="human5years")  {fVvenart <- fVvenart[ ,"human5years"]} 
  if(speciestype=="human10years") {fVvenart <- fVvenart[ ,"human10years"]} 
  if(speciestype=="humanm15years"){fVvenart <- fVvenart[ ,"humanm15years"]} 
  if(speciestype=="humanf15years"){fVvenart <- fVvenart[ ,"humanf15years"]} 
  if(speciestype=="humanm35years"){fVvenart <- fVvenart[ ,"humanm35years"]} 
  if(speciestype=="humanf35years"){fVvenart <- fVvenart[ ,"humanf35years"]} 
  
 return(fVvenart)

}



###############################################################################
#DATA:        dens (densisty of tissue of humans)                             
#UNIT:        [kg/L = g/cm^3]                                                 
#SOURCE:      human: all organs except for adi and bon: Brown et al. Toxicol. 
#                     Ind. Health 13: 407-484 (1997); Table 19                           
#                    adipose: ICRP, 1975. Report of the Task Group on Reference Man. 
#                     ICRP Publication 23. Pergamon Press, Oxford; p.44                          
#                    bone: ICRP, 2002. Basic Anatomical and Physiological Data 
#                     for Use in Radiological Protection Reference Values. ICRP 
#                     Publication 89. Ann. ICRP 32 (3-4); Table 2.20 (whole 
#                     skeleton, adults)                                                  
#ASSUMPTIONS: 1) rat:   data identical to human data                          
#             2) mouse: data identical to human data                          
###############################################################################
GetDensity <- function (speciestype) {
  
  dens <- matrix(1,nrow=14,ncol=3)
  dimnames(dens) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                           "hea","kid","liv","lun","mus","ski","spl"), 
                         c("rat","mouse","human"))     
  
  dens["adi", ] <- 0.916
  dens["bon", ] <- 1.3
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {dens <- dens[ ,"rat"]} 
  if(species=="mouse") {dens <- dens[ ,"mouse"]} 
  if(species=="human") {dens <- dens[ ,"human"]} 
  
  return(dens)
}


###############################################################################
#DATA:        BW (body weight)                                                
#UNIT:        [kg]                                                            
#SOURCE:      rat:   Brown et al. Toxicol. Ind. Health 13: 407-484 (1997); Table 3 
#                    (relation between weight and age)                                    
#             mouse: Brown et al. Toxicol. Ind. Health 13: 407-484 (1997); Table 1               
#             human: ICRP, 2002. Basic Anatomical and Physiological Data for Use
#                    in Radiological Protection Reference Values. ICRP Publication 
#                    89. Ann. ICRP 32 (3-4); Table 2.9                             
#STRAIN:      rat:   most data were reported for the F344 strain              
#             mouse: most data were reported for the B6C3F1 strain            
#NOTE:        rat:   Male rats of 250 g are only 10 weeks of age and in rapid  
#                    growth phase. Growth is much slower between age 20-55 weeks 
#                    (ca. 380-470g) and reaches a plateau of 475g at age 56-91 weeks.                                 
#             mouse: Male mice of 25 g are only 9 weeks of age and in rapid    
#                    growth phase. Growth is much slower between age 20-67 weeks 
#                    (ca. 31.5-40 g) and reaches a plateau of 40g at age 68-91 weeks.                               
#             human: As in the source, we associate an average age of 35 with 
#                    the adult (age 20 to age 50).                       
###############################################################################

GetBodyWeight <- function (speciestype) {
  
  BW <- matrix(NA,nrow=1,ncol=12)
  dimnames(BW) <- list("BW",c("rat10weeks",   "rat70weeks",
                              "mouse9weeks",  "mouse70weeks",
                              "humannewborn", "human1year",
                              "human5years",  "human10years",
                              "humanm15years","humanf15years",
                              "humanm35years","humanf35years"))
  
  #############   rat10weeks rat70weeks  mouse9weeks mouse70weeks humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years
  BW[1, ] <- c(   0.25,      0.475,      0.025,      0.040,       3.50,        10.00,     19.00,      32.00,       56.00,        53.00,        73.00,        60.00        )
  #                     between 56-91w
  
  if(speciestype=="rat10weeks")   {bodyweight <- BW[ ,"rat10weeks"]} 
  if(speciestype=="rat70weeks")   {bodyweight <- BW[ ,"rat70weeks"]} 
  if(speciestype=="mouse9weeks")  {bodyweight <- BW[ ,"mouse9weeks"]} 
  if(speciestype=="mouse70weeks") {bodyweight <- BW[ ,"mouse70weeks"]} 
  if(speciestype=="humannewborn") {bodyweight <- BW[ ,"humannewborn"]} 
  if(speciestype=="human1year")   {bodyweight <- BW[ ,"human1year"]} 
  if(speciestype=="human5years")  {bodyweight <- BW[ ,"human5years"]} 
  if(speciestype=="human10years") {bodyweight <- BW[ ,"human10years"]} 
  if(speciestype=="humanm15years"){bodyweight <- BW[ ,"humanm15years"]} 
  if(speciestype=="humanf15years"){bodyweight <- BW[ ,"humanf15years"]} 
  if(speciestype=="humanm35years"){bodyweight <- BW[ ,"humanm35years"]} 
  if(speciestype=="humanf35years"){bodyweight <- BW[ ,"humanf35years"]} 
  
  out <- list(bodyweight=bodyweight, BW=BW)
  
  return(out)
}


GetOrganWeightAndVolume <- function (speciestype) {
  
  #############################################################################
  #DATA:     fBWOW (fraction of total body weight that is experimental organ  
  #          weight or total blood volume) According to Brown et al. p.411,   
  #          in most cases, the values provided reflect the weight of organs  
  #          that are drained of blood.                       
  #UNIT:     [fraction] converted from [percentage]         
  #SOURCE:   rat:   Brown et al. Toxicol. Ind. Health 13: 407-484 (1997); 
  #                 most tissues: Table 5 
  #                 adi: Table 13 
  #                 bon: p.425  
  #                 blo: Diehl et al. J. Appl. Toxicol. 21: 15-23 (2001)          
  #          mouse: Brown et al. Toxicol. Ind. Health 13: 407-484 (1997); 
  #                 most tissues: Table 4
  #                 adi: Table 10
  #                 bon: mean of values reported on p.424 
  #                 blo: Diehl et al. J. Appl. Toxicol. 21: 15-23 (2001)                                                            
  #NOTE:     gut is sum of stomach, small and large intestine                 
  #############################################################################
  
  fOWBW_experimental <- matrix(NA,nrow=14,ncol=12)
  dimnames(fOWBW_experimental) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                                         "hea","kid","liv","lun","mus","ski","spl"),
                                       c("rat10weeks",   "rat70weeks",
                                         "mouse9weeks",  "mouse70weeks",
                                         "humannewborn", "human1year",
                                         "human5years",  "human10years",
                                         "humanm15years","humanf15years",
                                         "humanm35years","humanf35years"))
  
  #################################   rat10weeks  rat70weeks  mouse9weeks mouse70weeks humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years
  fOWBW_experimental["blo", ] <- c(    6.4 /100,   6.4 /100,   7.2 /100,   7.2 /100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["pla", ] <- c(          NA,         NA,         NA,         NA,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["ery", ] <- c(          NA,         NA,         NA,         NA,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["adi", ] <- c(   10.6 /100,  16.0 /100,   7   /100,   7   /100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["bon", ] <- c(    7.3 /100,   7.3 /100,  10.73/100,  10.73/100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["bra", ] <- c(    0.57/100,   0.57/100,   1.65/100,   1.65/100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["gut", ] <- c(    2.7 /100,   2.7 /100,   4.22/100,   4.22/100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["hea", ] <- c(    0.33/100,   0.33/100,   0.50/100,   0.50/100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["kid", ] <- c(    0.73/100,   0.73/100,   1.67/100,   1.67/100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["liv", ] <- c(    3.66/100,   3.66/100,   5.49/100,   5.49/100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["lun", ] <- c(    0.50/100,   0.50/100,   0.73/100,   0.73/100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["mus", ] <- c(   40.43/100,  40.43/100,  38.40/100,  38.40/100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["ski", ] <- c(   19.03/100,  19.03/100,  16.53/100,  16.53/100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fOWBW_experimental["spl", ] <- c(    0.2 /100,   0.2 /100,   0.35/100,   0.35/100,   NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  
  #Calculate experimental organ weights
  #UNIT: [kg]=[fraction*kg]
  BW       <- GetBodyWeight(speciestype)
  OW_experimental <- sweep(fOWBW_experimental,MARGIN=2,BW$BW ,'*')
  
  #############################################################################
  #DATA:     OW (experimental organ weights for human)                        
  #UNIT:     [kg] converted from [g]                       
  #SOURCE:   ICRP, 2002. Basic Anatomical and Physiological Data for Use in 
  #          Radiological Protection Reference Values. ICRP Publication 89. 
  #          Ann. ICRP 32 (3-4); Table 2.8                                     
  #NOTE:     Tissue only (NOT INCLUDING VASCULAR SPACE)                       
  #          separable adipose, total skeleton, gut is sum of small intestine 
  #          wall, large intestine right colon wall, large intestine left     
  #          colon wall and large intestine rectosigmoid wall, heart tissue   
  #          only, lung tissue only                                           
  #############################################################################
  
  OWhuman_experimental <- OW_experimental[,5:12]
  
  ##############################        humannewborn human1year  human5years human10years humanm15years humanf15years humanm35years humanf35years
  OWhuman_experimental["blo", ] <- c(   270  /1000,   500/1000,  1400/1000,   2400/1000,   4500/1000,    3300/1000,    5300/1000,    3900/1000   )
  OWhuman_experimental["pla", ] <- c(   NA,          NA,         NA,          NA,          NA,           NA,           NA,           NA          )
  OWhuman_experimental["ery", ] <- c(   NA,          NA,         NA,          NA,          NA,           NA,           NA,           NA          )
  OWhuman_experimental["adi", ] <- c(   890  /1000,  3600/1000,  5000/1000,   7500/1000,   9500/1000,   16000/1000,   14500/1000,   19000/1000   )
  OWhuman_experimental["bon", ] <- c(   370  /1000,  1170/1000,  2430/1000,   4500/1000,   7950/1000,    7180/1000,   10500/1000,    7800/1000   )
  OWhuman_experimental["bra", ] <- c(   380  /1000,   950/1000,  1245/1000,   1310/1000,   1420/1000,    1300/1000,    1450/1000,    1300/1000   )
  OWhuman_experimental["gut", ] <- c(    47  /1000,   135/1000,   340/1000,    580/1000,    820/1000,     820/1000,    1020/1000,     960/1000   )
  OWhuman_experimental["hea", ] <- c(    20  /1000,    50/1000,    85/1000,    140/1000,    230/1000,     220/1000,     330/1000,     250/1000   )
  OWhuman_experimental["kid", ] <- c(    25  /1000,    70/1000,   110/1000,    180/1000,    250/1000,     240/1000,     310/1000,     275/1000   )
  OWhuman_experimental["liv", ] <- c(   130  /1000,   330/1000,   570/1000,    830/1000,   1300/1000,    1300/1000,    1800/1000,    1400/1000   )
  OWhuman_experimental["lun", ] <- c(    30  /1000,    80/1000,   125/1000,    210/1000,    330/1000,     290/1000,     500/1000,     420/1000   )
  OWhuman_experimental["mus", ] <- c(   800  /1000,  1900/1000,  5600/1000,  11000/1000,  24000/1000,   17000/1000,   29000/1000,   17500/1000   )
  OWhuman_experimental["ski", ] <- c(   175  /1000,   350/1000,   570/1000,    820/1000,   2000/1000,    1700/1000,    3300/1000,    2300/1000   )
  OWhuman_experimental["spl", ] <- c(     9.5/1000,    29/1000,    50/1000,     80/1000,    130/1000,     130/1000,     150/1000,     130/1000   )
  
  #Merge OWhuman_experimental and OW_experimental
  OW_experimental[,5:12] <- OWhuman_experimental
  
  #ASSUMPTION: We assume that the experimental organ weights (including residual
  #            blood to some varying degree) are approximately equal to the
  #            tissue organ weight (not including residual blood), since
  #            according to Brown et al. p.411, in most cases, the values
  #            provided reflect the weight of organs that are drained of blood.

  OW <- OW_experimental
  
  if(speciestype=="rat10weeks")   {OW <- OW[ ,"rat10weeks"]} 
  if(speciestype=="rat70weeks")   {OW <- OW[ ,"rat70weeks"]} 
  if(speciestype=="mouse9weeks")  {OW <- OW[ ,"mouse9weeks"]} 
  if(speciestype=="mouse70weeks") {OW <- OW[ ,"mouse70weeks"]} 
  if(speciestype=="humannewborn") {OW <- OW[ ,"humannewborn"]} 
  if(speciestype=="human1year")   {OW <- OW[ ,"human1year"]} 
  if(speciestype=="human5years")  {OW <- OW[ ,"human5years"]} 
  if(speciestype=="human10years") {OW <- OW[ ,"human10years"]} 
  if(speciestype=="humanm15years"){OW <- OW[ ,"humanm15years"]} 
  if(speciestype=="humanf15years"){OW <- OW[ ,"humanf15years"]} 
  if(speciestype=="humanm35years"){OW <- OW[ ,"humanm35years"]} 
  if(speciestype=="humanf35years"){OW <- OW[ ,"humanf35years"]} 
  
  #Organ volumes
  #UNIT: [L]=[kg/kg/L]
  dens <- GetDensity(speciestype)  
  V    <- OW/dens
  
  
  #Get volumes for plasma and erythrocytes
  #############################################################################
  #DATA:    V (volumes of plasma and erythrocytes)                                                                                             
  #UNIT:    [L] converted from [mL]                                                                                           
  #SOURCE:  ICRP, 2002. Basic Anatomical and Physiological Data for Use in 
  #         Radiological Protection Reference Values. ICRP Publication 89. 
  #         Ann. ICRP 32 (3-4); Chapter 7.3                                            
  #############################################################################
  
  GetOrganVolumesPlasmaErythrocytes <- function(speciestype){
    switch(speciestype,
           "humanm35years" ={
             V[["pla"]] <- 3000/1000
             V[["ery"]] <- 2300/1000},
           "humanf35years" ={
             V[["pla"]] <- 2400/1000
             V[["ery"]] <- 1500/1000}
    )
    return(V)
  }
  
  V <- GetOrganVolumesPlasmaErythrocytes(speciestype)
  
  out <- list(OW=OW, V=V, OW_experimental=OW_experimental)
  
  return (out)
  
}





###############################################################################
#DATA:        co and fQco (cardiac output and fraction of co that is regional 
#                          blood flows)                                       
#UNIT:        co [L/h] converted from [L/min] and fQco [fraction]        
#SOURCE (co): rat:   Brown et al. Toxicol. Ind. Health 13: 407-484 (1997); p.441 
#                    (based on Arms and Travis, Reference Physiological Parameters 
#                    in Pharmacokinetic Modeling. U.S. EPA Final Report, 
#                    EPA 600/6-88/004 (1988), involving different scaling steps)     
#             mouse: Brown et al. Toxicol. Ind. Health 13: 407-484 (1997); p.440                  
#             human: ICRP, 2002. Basic Anatomical and Physiological Data for Use 
#                    in Radiological Protection Reference Values. ICRP 
#                    Publication 89. Ann. ICRP 32 (3-4); Table 2.39. 
#                    For newborn and children age 1, the values of Alverson et al.
#                    J Ultrasound Med 6: 519-524 (1987) (cited in Abraham et al., 
#                    Arch Toxicol. 79: 63-73 (2005)) have been taken, since the
#                    values of the ICRP report appear to be too low 
#                    (age1: ICRP.co=1.2, Alverson.co = 1.8. For newborn the 
#                    difference is 0.6 vs. 0.7. The ICRP value for age 1 is a
#                    marked outlier from the expectations of allometric scaling.                                      
#SOURCE (fQco): rat:   Brown et al. Toxicol. Ind. Health 13: 407-484 (1997); Table 25            
#                      spl: Poulin and Theil, J Pharm Sci. 91: 129-156 (2002)                  
#               mouse: Brown et al. Toxicol. Ind. Health 13: 407-484 (1997); Table 24            
#                      adi, bon, bra, gut, spl: El-Masri and Portier. Drug Metab 
#                      Dispos.26:585-594 (1998)                 
#               human: ICRP, 2002. Basic Anatomical and Physiological Data for 
#                      Use in Radiological Protection Reference Values. ICRP 
#                      Publication 89. Ann. ICRP 32 (3-4); Table 2.40                          
###############################################################################

GetCardiacOutputAndBloodFlows <- function(speciestype) {

  #--- CARDIAC OUTPUT ---------------------------------------------------------
  BW <- GetBodyWeight(speciestype)
  
  co_experimental <- matrix(NA,nrow=1,ncol=12)
  dimnames(co_experimental) <- list("co_experimental",c("rat10weeks",   "rat70weeks",
                                                        "mouse9weeks",  "mouse70weeks",
                                                        "humannewborn", "human1year",
                                                        "human5years",  "human10years",
                                                        "humanm15years","humanf15years",
                                                        "humanm35years","humanf35years"))
  
  ##########################   rat10weeks                            rat70weeks                            mouse9weeks                            mouse70weeks                            humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years
  co_experimental[1, ] <- c(   0.235*(BW$BW[,"rat10weeks"]^0.75)*60, 0.235*(BW$BW[,"rat70weeks"]^0.75)*60, 0.275*(BW$BW[,"mouse9weeks"]^0.75)*60, 0.275*(BW$BW[,"mouse70weeks"]^0.75)*60, 0.6*60,      1.2*60,    3.4*60,     5.0*60,      6.1*60,       6.1*60,       6.5*60,       5.9*60          )
  
  if(speciestype=="rat10weeks")   {co <- co_experimental[ ,"rat10weeks"]} 
  if(speciestype=="rat70weeks")   {co <- co_experimental[ ,"rat70weeks"]} 
  if(speciestype=="mouse9weeks")  {co <- co_experimental[ ,"mouse9weeks"]} 
  if(speciestype=="mouse70weeks") {co <- co_experimental[ ,"mouse70weeks"]} 
  if(speciestype=="humannewborn") {co <- co_experimental[ ,"humannewborn"]} 
  if(speciestype=="human1year")   {co <- co_experimental[ ,"human1year"]} 
  if(speciestype=="human5years")  {co <- co_experimental[ ,"human5years"]} 
  if(speciestype=="human10years") {co <- co_experimental[ ,"human10years"]} 
  if(speciestype=="humanm15years"){co <- co_experimental[ ,"humanm15years"]} 
  if(speciestype=="humanf15years"){co <- co_experimental[ ,"humanf15years"]} 
  if(speciestype=="humanm35years"){co <- co_experimental[ ,"humanm35years"]} 
  if(speciestype=="humanf35years"){co <- co_experimental[ ,"humanf35years"]} 
  
  #--- FRACTION OF CARDIAC OUTPUT THAT IS REGIONL BLOOD FLOW ------------------
  fQco_experimental <- matrix(NA,nrow=14,ncol=12)
  dimnames(fQco_experimental) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                                        "hea","kid","liv","lun","mus","ski","spl"),
                                      c("rat10weeks",   "rat70weeks",
                                        "mouse9weeks",  "mouse70weeks",
                                        "humannewborn", "human1year",
                                        "human5years",  "human10years",
                                        "humanm15years","humanf15years",
                                        "humanm35years","humanf35years"))
  
  ################################   rat10weeks  rat70weeks  mouse9weeks  mouse70weeks humannewborn human1year human5years human10years humanm15years humanf15years humanm35years humanf35years
  fQco_experimental["blo", ] <- c(   NA,         NA,         NA,          NA,          NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fQco_experimental["pla", ] <- c(   NA,         NA,         NA,          NA,          NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fQco_experimental["ery", ] <- c(   NA,         NA,         NA,          NA,          NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fQco_experimental["adi", ] <- c(    7.0/100,   7.0/100,     7  /100,     7  /100,    NA,          NA,        NA,         NA,          NA,           NA,            5.00/100,     8.50/100    )
  fQco_experimental["bon", ] <- c(   12.2/100,   12.2/100,   11  /100,    11  /100,    NA,          NA,        NA,         NA,          NA,           NA,            5.00/100,     5.00/100    )
  fQco_experimental["bra", ] <- c(    2.0/100,   2.0/100,     3.3/100,     3.3/100,    NA,          NA,        NA,         NA,          NA,           NA,           12.00/100,    12.00/100    )
  fQco_experimental["gut", ] <- c(   13.1/100,   13.1/100,   14.1/100,    14.1/100,    NA,          NA,        NA,         NA,          NA,           NA,           14.00/100,    16.00/100    ) #rat, gut: Reference?
  fQco_experimental["hea", ] <- c(    4.9/100,   4.9/100,     6.6/100,     6.6/100,    NA,          NA,        NA,         NA,          NA,           NA,            4.00/100,     5.00/100    )
  fQco_experimental["kid", ] <- c(   14.1/100,   14.1/100,    9.1/100,     9.1/100,    NA,          NA,        NA,         NA,          NA,           NA,           19.00/100,    17.00/100    )
  fQco_experimental["liv", ] <- c(   17.4/100,   17.4/100,   16.2/100,    16.2/100,    NA,          NA,        NA,         NA,          NA,           NA,           25.50/100,    27.00/100    )
  fQco_experimental["lun", ] <- c(   NA,         NA,         NA,          NA,          NA,          NA,        NA,         NA,          NA,           NA,           NA,           NA           )
  fQco_experimental["mus", ] <- c(   27.8/100,   27.8/100,   15.9/100,    15.9/100,    NA,          NA,        NA,         NA,          NA,           NA,           17.00/100,    12.00/100    )
  fQco_experimental["ski", ] <- c(    5.8/100,   5.8/100,     5.8/100,     5.8/100,    NA,          NA,        NA,         NA,          NA,           NA,            5.00/100,     5.00/100    )
  fQco_experimental["spl", ] <- c(    2.0/100,   2.0/100,     1  /100,     1  /100,    NA,          NA,        NA,         NA,          NA,           NA,            3.00/100,     3.00/100    )
  
  #Regional blood flow [L/h] by multiplication of fractional blood flow
  #with cardiac output
  Q <- sweep(fQco_experimental,MARGIN=2,co_experimental,'*')
   
  
  #--- CHILDREN ---------------------------------------------------------------
  
  if((speciestype=="humannewborn")|(speciestype=="human1year")|(speciestype=="human5years")|
     (speciestype=="human10years")|(speciestype=="humanm15years")|(speciestype=="humanf15years")) {
     #############################################################################
     # CHILDREN: For all children age classes: fQco values were estimated based  
     #           on the approach presented in Abraham et al. Arch Toxicol. 79:
     #           63-73 (2005), except when experimental data are available, 
     #           i.e. brain and kidney                                
     #############################################################################  
     OW_and_V <- GetOrganWeightAndVolume(speciestype)
     dens     <- GetDensity(speciestype)
  
     Qchildren <- Q[,5:10]
     # Intermediate target regional blood flow (Qinter) scaled solely according to 
     # ratio of target-to-reference organ weights (SF)
     # Reference individuals: humanf35years or humanm35years
  
       # Ratio of target-to-reference organ weights (SF)
        if(speciestype=="humanm15years") {
          OWref <- OW_and_V$OW_experimental[,"humanm35years"]
          SF    <- OW_and_V$OW/OWref
        }
        if((speciestype=="humannewborn")|(speciestype=="human1year")|(speciestype=="human5years")|(speciestype=="human10years")|(speciestype=="humanf15years")) {
          OWref <- OW_and_V$OW_experimental[,"humanf35years"]
          SF    <- OW_and_V$OW/OWref
        }
  
        # Intermediate target regional blood flow (Qinter)
        if(speciestype=="humanm15years") {
          Qref   <- fQco_experimental[,"humanm35years"]*co_experimental[,"humanm35years"]
          Qinter <- SF*Qref 
        }
        if((speciestype=="humannewborn")|(speciestype=="human1year")|(speciestype=="human5years")|(speciestype=="human10years")|(speciestype=="humanf15years")) {
          Qref   <- fQco_experimental[,"humanf35years"]*co_experimental[,"humanf35years"]
          Qinter <- SF*Qref 
        }
  
     #--- CHILDREN BRAIN ------------------------------------------------------
     # Blood flow data in [L/h] converted from [mL/min/100g] from Chiron et al. 
     # J. Nucl. Med. 33:696-703 (1992)
  
     #######################   humannewborn                      human1year                         human5years                        human10years                       humanm15years                      humanf15years
     Qchildren["bra",]  <- c(  50/1000/10*OW_and_V$V[["bra"]]*60,59/1000/10*OW_and_V$V[["bra"]]*60, 71/1000/10*OW_and_V$V[["bra"]]*60, 68/1000/10*OW_and_V$V[["bra"]]*60, 57/1000/10*OW_and_V$V[["bra"]]*60, 57/1000/10*OW_and_V$V[["bra"]]*60  )


     #--- CHILDREN KIDNEY -----------------------------------------------------
     # Assume that blood flow per kg kidney tissue is independent of age according
     # to Grunert et al. Eur J Pediatr. 149:287-292 (1990)
  
     #######################   humannewborn   human1year     human5years    human10years   humanm15years  humanf15years
     Qchildren["kid",]  <- c(  Qinter["kid"], Qinter["kid"], Qinter["kid"], Qinter["kid"], Qinter["kid"], Qinter["kid"]  )
  
  
     #--- CHILDREN ADIPOSE, BONE, GUT, HEART, LIVER, MUSCLE, SKIN, SPLEEN -----
   
     #Scale regional tissue blood flows to match cardiac output. 
  
     #1) Scale all blood flows that flow into the vein, excluding those tissues 
     #   where experimental data were used (brain and kidney), 
     #   i.e. adipose, bone, heart, liver, muscle, ski
     SF_co <- (co_experimental[,5:10] - Qchildren["kid",] - Qchildren["bra",])/
              (Qinter["adi"]+Qinter["bon"]+Qinter["bra"]+Qinter["gut"]+Qinter["hea"]+
               Qinter["kid"]+Qinter["liv"]+Qinter["mus"]+Qinter["ski"]+Qinter["spl"])
  
     Qchildren["adi",]  <- SF_co*Qinter["adi"]
     Qchildren["bon",]  <- SF_co*Qinter["bon"]
     Qchildren["hea",]  <- SF_co*Qinter["hea"]
     Qchildren["liv",]  <- SF_co*Qinter["liv"]
     Qchildren["mus",]  <- SF_co*Qinter["mus"]
     Qchildren["ski",]  <- SF_co*Qinter["ski"]

     #2) Scale regional tissue blood flows that flow into the liver,
     #   i.e. gut and spleen
     #   To this end, determine hepartic artery blood flow (hepart)
  
     # Intermediate target hepatic artery blood flow (Qinter_hepart)
     if(speciestype=="humanm15years") {
       OVref         <- OW_and_V$OW_experimental[,"humanm35years"]/dens
       OVchildren_liv<- OW_and_V$OW_experimental["liv",5:10]/dens["liv"]
       Qhepart       <- Q["liv","humanm35years"]-Q["gut","humanm35years"]-Q["spl","humanm35years"]
       Qhepart_perkg <- Qhepart/OVref["liv"]
       Qinter_hepart <- Qhepart_perkg*OVchildren_liv
     }
  
     if((speciestype=="humannewborn")|(speciestype=="human1year")|(speciestype=="human5years")|(speciestype=="human10years")|(speciestype=="humanf15years")) {
       OVref         <- OW_and_V$OW_experimental[,"humanf35years"]/dens
       OVchildren_liv<- OW_and_V$OW_experimental["liv",5:10]/dens["liv"]
       Qhepart       <- Q["liv","humanf35years"]-Q["gut","humanf35years"]-Q["spl","humanf35years"]
       Qhepart_perkg <- Qhepart/OVref["liv"]
       Qinter_hepart <- Qhepart_perkg*OVchildren_liv
     }
  
     SF_hepart <- Qchildren["liv",]/(Qinter_hepart+Qinter["gut"]+Qinter["spl"])
  
     Qchildren["gut",] <- SF_hepart*Qinter["gut"]
     Qchildren["spl",] <- SF_hepart*Qinter["spl"]

     # --- MERGE CHILDREN DATA AND Q ------------------------------------------
     Q[,5:10] <- Qchildren
  }
  
  
  
  if(speciestype=="rat10weeks")   {Q <- Q[ ,"rat10weeks"]} 
  if(speciestype=="rat70weeks")   {Q <- Q[ ,"rat70weeks"]} 
  if(speciestype=="mouse9weeks")  {Q <- Q[ ,"mouse9weeks"]} 
  if(speciestype=="mouse70weeks") {Q <- Q[ ,"mouse70weeks"]} 
  if(speciestype=="humannewborn") {Q <- Q[ ,"humannewborn"]} 
  if(speciestype=="human1year")   {Q <- Q[ ,"human1year"]} 
  if(speciestype=="human5years")  {Q <- Q[ ,"human5years"]} 
  if(speciestype=="human10years") {Q <- Q[ ,"human10years"]} 
  if(speciestype=="humanm15years"){Q <- Q[ ,"humanm15years"]} 
  if(speciestype=="humanf15years"){Q <- Q[ ,"humanf15years"]} 
  if(speciestype=="humanm35years"){Q <- Q[ ,"humanm35years"]} 
  if(speciestype=="humanf35years"){Q <- Q[ ,"humanf35years"]} 
  
  out <- list(co=co, Q=Q)
  return(out)
}


###############################################################################
#DATA:        hct (hematocrit as fraction of red blood cells-to-blood volume) 
#UNIT:        [-]                                                             
#SOURCE:      Windberger and Baskurt, Handbook of hemorheology and hemodynamics, 
#             1st ed. Amsterdam: IOS Press (2007)                                               
#ASSUMPTION:  hct value of is function of age (only defined for age           
#             "m35years" and "f35years")                                      
###############################################################################

GetHematocrit <- function(speciestype) {
  species <- GetSpecies(speciestype) 
  
  hct <- NA
  switch(species,
         "rat"  ={hct <- 0.40},
         "mouse"={hct <- 0.40},
         "human"={hct <- 0.40}
  )
  return(hct)
}


#--- SPECIES-SPECIFIC DATA RODGERS AND ROWLAND --------------------------------

GetSpeciesDataRodgersAndRowland <- function (speciestype) {

  #############################################################################
  #DATA:        pH (pH values in plasma (7.4), interstitital water (7.4), intra-        
  #             cellular water (7.0) and ery (7.22)  in rat)                                          
  #UNIT:        [-]                                                           
  #SOURCE:      Rodgers et al. J Pharm Sci. 94: 1259-1276 (2005); p. 1263      
  #ASSUMPTION:  pH values assumed to be tissue independent, as in Rodgers et al.                                                        
  #ASSUMPTION:  pH values of human and mouse assumed to be identical to rat values                                                    
  #############################################################################

  pH <- matrix(7,nrow=14,ncol=3)
  dimnames(pH) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                         "hea","kid","liv","lun","mus","ski","spl"),
                       c("rat","mouse","human"))   

  pH["pla", ] <- 7.4
  pH["ery", ] <- 7.22

  species <- GetSpecies(speciestype)
  if(species=="rat")   {pH <- pH[ ,"rat"]} 
  if(species=="mouse") {pH <- pH[ ,"mouse"]} 
  if(species=="human") {pH <- pH[ ,"human"]} 

  #############################################################################
  #DATA:        fV_w_tot_experimental (Fraction of tissue volume that is      
  #             total tissue water)                                           
  #UNIT:        [fraction]                                                    
  #SOURCE:      rat:    Rodgers and Rowland, J Pharm Sci. 95: 1238-1257 (2006); 
  #                     Table 1, total tissue water is reported, corrected for         
  #                     residual (see eq. (A2)).                              
  #                     According to email correspondance of W.Huisinga with  
  #                     T. Rodgers, f_residual was taken from ref. 15         
  #                     (Kawai et al. 1994), not from ref. 16 as mentioned    
  #                     in the article)                                                                           
  #             human:  Poulin and Theil, J Pharm Sci. 91: 129-156 (2002)
  #                     ery: Poulin and Theil, J Pharm Sci. 98: 4941-4961 (2009)           
  #ASSUMPTIONS: 1) mouse values are identical to rat values                   
  #             2) fractions of experimental tissue volumes are identical to  
  #                the fractions of tissue volume                             
  #############################################################################

  fV_w_tot_experimental <- matrix(NA,nrow=14,ncol=3)
  dimnames(fV_w_tot_experimental) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                                            "hea","kid","liv","lun","mus","ski","spl"), 
                                          c("rat","mouse","human"))

  #####################################  rat     mouse   human 
  fV_w_tot_experimental["blo", ] <- c(   NA,     NA,     NA    )
  fV_w_tot_experimental["pla", ] <- c(   NA,     NA,     0.945 )
  fV_w_tot_experimental["ery", ] <- c(   NA,     NA,     0.63  )
  fV_w_tot_experimental["adi", ] <- c(   0.144,  0.144,  0.18  )
  fV_w_tot_experimental["bon", ] <- c(   0.417,  0.417,  0.439 )
  fV_w_tot_experimental["bra", ] <- c(   0.753,  0.753,  0.77  )
  fV_w_tot_experimental["gut", ] <- c(   0.738,  0.738,  0.718 )
  fV_w_tot_experimental["hea", ] <- c(   0.568,  0.568,  0.758 )
  fV_w_tot_experimental["kid", ] <- c(   0.672,  0.672,  0.783 )
  fV_w_tot_experimental["liv", ] <- c(   0.642,  0.642,  0.751 )
  fV_w_tot_experimental["lun", ] <- c(   0.574,  0.574,  0.811 )
  fV_w_tot_experimental["mus", ] <- c(   0.726,  0.726,  0.76  )
  fV_w_tot_experimental["ski", ] <- c(   0.658,  0.658,  0.718 )
  fV_w_tot_experimental["spl", ] <- c(   0.562,  0.562,  0.778 )

  fV_w_tot_experimental_human <- fV_w_tot_experimental[ ,"human"]
 
  species <- GetSpecies(speciestype)
  if(species=="rat")   {fV_w_tot_experimental <- fV_w_tot_experimental[ ,"rat"]} 
  if(species=="mouse") {fV_w_tot_experimental <- fV_w_tot_experimental[ ,"mouse"]} 
  if(species=="human") {fV_w_tot_experimental <- fV_w_tot_experimental[ ,"human"]} 

  #ASSUMPTION: fractions of experimental tissue volumes are identical to the 
  #            fractions of tissue volume (since Vtis is approximately equal 
  #            to Vexperimental)
  fV_w_tot <- fV_w_tot_experimental

  #############################################################################
  #DATA:        fV_w_ex_experimental (Fraction of experimental tissue volumes 
  #             (based on the tissue wet weight) that is extracellular tissue 
  #             water)                                                        
  #UNIT:        [fraction]                                                    
  #SOURCE:      rat:  Rodgers et al. J Pharm Sci. 94: 1259-1276 (2005)     
  #ASSUMPTIONS: 1) mouse:  values are identical to rat values               
  #             2) human:  Assume that ratio w_ex-to-w_tot is the same as in rat                                                
  #############################################################################

  fV_w_ex_experimental <- matrix(NA,nrow=14,ncol=3)
  dimnames(fV_w_ex_experimental) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                                           "hea","kid","liv","lun","mus","ski","spl"), 
                                         c("rat","mouse","human"))

  ####################################  rat     mouse   human 
  fV_w_ex_experimental["blo", ] <- c(   NA,     NA,     NA    )
  fV_w_ex_experimental["pla", ] <- c(   NA,     NA,     NA    )
  fV_w_ex_experimental["ery", ] <- c(   0,      0,      NA    )
  fV_w_ex_experimental["adi", ] <- c(   0.135,  0.135,  NA    )
  fV_w_ex_experimental["bon", ] <- c(   0.100,  0.100,  NA    )
  fV_w_ex_experimental["bra", ] <- c(   0.162,  0.162,  NA    )
  fV_w_ex_experimental["gut", ] <- c(   0.282,  0.282,  NA    )
  fV_w_ex_experimental["hea", ] <- c(   0.320,  0.320,  NA    )
  fV_w_ex_experimental["kid", ] <- c(   0.273,  0.273,  NA    )
  fV_w_ex_experimental["liv", ] <- c(   0.161,  0.161,  NA    )
  fV_w_ex_experimental["lun", ] <- c(   0.336,  0.336,  NA    )
  fV_w_ex_experimental["mus", ] <- c(   0.118,  0.118,  NA    )
  fV_w_ex_experimental["ski", ] <- c(   0.382,  0.382,  NA    )
  fV_w_ex_experimental["spl", ] <- c(   0.207,  0.207,  NA    )

  fV_w_ex_experimental_rat <- fV_w_ex_experimental[ ,"rat"]

  species <- GetSpecies(speciestype)
  if(species=="rat")   {fV_w_ex_experimental <- fV_w_ex_experimental[ ,"rat"]} 
  if(species=="mouse") {fV_w_ex_experimental <- fV_w_ex_experimental[ ,"mouse"]} 
  if(species=="human") {fV_w_ex_experimental <- fV_w_ex_experimental[ ,"human"]} 

  #############################################################################
  #DATA:        fV_w_ic_experimental (Fraction of experimental organ weight   
  #             that is intracellular tissue water)                           
  #UNIT:        [fraction]                                                    
  #SOURCE:      rat:    Rodgers et al. J Pharm Sci. 94: 1259-1276 (2005),      
  #                     fraction of intra-cellular water was determined as    
  #                     difference of extra-cellular water to total           
  #                     experimental tissue water                             
  #                     i.e. fV_w_ic_experimental = fV_w_tot_experimental -   
  #                          fV_w_ex_experimental                             
  #ASSUMPTIONS: 1) mouse:  values are identical to rat values                 
  #             2) human:  assume that ratio wex-to-wtot is the same as in rat
  #             3) fractions of experimental tissue volumes are identical to
  #                the fractions of tissue volume                             
  #############################################################################

  fV_w_ic_experimental <- matrix(NA,nrow=14,ncol=3)
  dimnames(fV_w_ic_experimental) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                                           "hea","kid","liv","lun","mus","ski","spl"), 
                                         c("rat","mouse","human"))

  ####################################  rat     mouse   human 
  fV_w_ic_experimental["blo", ] <- c(   NA,     NA,     NA    )
  fV_w_ic_experimental["pla", ] <- c(   NA,     NA,     NA    )
  fV_w_ic_experimental["ery", ] <- c(   0.603,  0.603,  NA    )
  fV_w_ic_experimental["adi", ] <- c(   0.017,  0.017,  NA    )
  fV_w_ic_experimental["bon", ] <- c(   0.346,  0.346,  NA    )
  fV_w_ic_experimental["bra", ] <- c(   0.620,  0.620,  NA    )
  fV_w_ic_experimental["gut", ] <- c(   0.475,  0.475,  NA    )
  fV_w_ic_experimental["hea", ] <- c(   0.456,  0.456,  NA    )
  fV_w_ic_experimental["kid", ] <- c(   0.483,  0.483,  NA    )
  fV_w_ic_experimental["liv", ] <- c(   0.573,  0.573,  NA    )
  fV_w_ic_experimental["lun", ] <- c(   0.446,  0.446,  NA    )
  fV_w_ic_experimental["mus", ] <- c(   0.630,  0.630,  NA    )
  fV_w_ic_experimental["ski", ] <- c(   0.291,  0.291,  NA    )
  fV_w_ic_experimental["spl", ] <- c(   0.579,  0.579,  NA    )

  fV_w_ic_experimental_rat <- fV_w_ic_experimental[ ,"rat"]

  species <- GetSpecies(speciestype)
  if(species=="rat")   {fV_w_ic_experimental <- fV_w_ic_experimental[ ,"rat"]} 
  if(species=="mouse") {fV_w_ic_experimental <- fV_w_ic_experimental[ ,"mouse"]} 
  if(species=="human") {fV_w_ic_experimental <- fV_w_ic_experimental[ ,"human"]} 

  #ASSUMPTION: fractions of experimental tissue volumes are identical to the 
  #            fractions of tissue volume (since Vtis is approximately equal 
  #            to Vexperimental)
  fV_w_ic <- fV_w_ic_experimental

  #As stated above, fV_w_ic_experimental = fV_w_tot_experimental - fV_w_ex_experimental. 
  #Assuming that both, fV_w_tot_experimental and fV_w_ex_experimental are perturbed by 
  #water in residual blood and to the same extent, residual blood perturbations 
  #cancel out in fV_w_ic_experimental.
 
  #We then determine again fV_w_ex by the difference of fV_w_ic to fV_w_tot.
  fV_w_ex <- fV_w_tot - fV_w_ic

  #############################################################################
  #Adapt fV_w_ex and fV_w_ic for "human"                                      
  #############################################################################
  species <- GetSpecies(speciestype)
  if (species=="human") {
  
    #assume that ratio wex-to-wtot is the same as in rat
    ratio <- fV_w_ex_experimental_rat/ (fV_w_ex_experimental_rat + fV_w_ic_experimental_rat)
    ratio
    #fraction of extra-celluar water 
    fV_w_ex_human <- ratio*fV_w_tot_experimental_human
    fV_w_ex_human
    #replace fV_w_ex by values calculated for human
    fV_w_ex <- fV_w_ex_human
    #fraction of intra-celluar water
    fV_w_ic_human <- fV_w_tot_experimental_human-fV_w_ex_human
    #replace fV_w_ic by values calculated for human
    fV_w_ic <- fV_w_ic_human
  }  

  #############################################################################
  #DATA:        fV_nl (fraction of tissue volume that is neutral lipids)                                                                         
  #UNIT:        [fraction]                                                                                                                       
  #SOURCE:      rat:   Rodgers and Rowland, J Pharm Sci. 95: 1238-1257 (2006)                                                                                  
  #                    Values have been corrected for residual blood          
  #                    contributions, i.e., reported values do not contain    
  #                    resdiual blood contributions (see Rodgers and Rowland, 
  #                    p.1252 top left paragraph)                                             
  #                    adi: in Rodgers and Rowland (2006) incorrectly         
  #                         reported under neutral phospholipids            
  #                    pla: p. 1241, paragarph "Tissue Specific Input          
  #                         Parameters"                                       
  #                    ery: Rodgers et al. J Pharm Sci. 94: 1259-1276 (2005); Table 1                                                                              
  #             human: Poulin and Theil, J Pharm Sci. 91: 129-156 (2002)                                                                              
  #                    ery: Rodgers and Rowland, Pharm Res. 24: 918-933 (2007); Table VII                         
  #ASSUMPTIONS: 1) mouse values are identical to rat values                                                                                       
  #             2) fractions of experimental tissue volumes are identical to  
  #                the fractions of tissue volume                             
  #############################################################################

  fV_nl <- matrix(NA,nrow=14,ncol=3)
  dimnames(fV_nl) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                            "hea","kid","liv","lun","mus","ski","spl"), 
                          c("rat","mouse","human"))

  ####################   rat     mouse   human 
  fV_nl["blo", ] <- c(   NA,     NA,     NA      )
  fV_nl["pla", ] <- c(   0.0023, 0.0023, 0.0035  )
  fV_nl["ery", ] <- c(   0.0017, 0.0017, 0.0033  )
  fV_nl["adi", ] <- c(   0.8530, 0.8530, 0.79    )
  fV_nl["bon", ] <- c(   0.0174, 0.0174, 0.074   )
  fV_nl["bra", ] <- c(   0.0391, 0.0391, 0.051   )
  fV_nl["gut", ] <- c(   0.0375, 0.0375, 0.0487  )
  fV_nl["hea", ] <- c(   0.0135, 0.0135, 0.0115  )
  fV_nl["kid", ] <- c(   0.0121, 0.0121, 0.0207  )
  fV_nl["liv", ] <- c(   0.0135, 0.0135, 0.0348  )
  fV_nl["lun", ] <- c(   0.0215, 0.0215, 0.003   )
  fV_nl["mus", ] <- c(   0.0100, 0.0100, 0.0238  )
  fV_nl["ski", ] <- c(   0.0603, 0.0603, 0.0284  )
  fV_nl["spl", ] <- c(   0.0071, 0.0071, 0.0201  )

  species <- GetSpecies(speciestype)
  if(species=="rat")   {fV_nl <- fV_nl[ ,"rat"]} 
  if(species=="mouse") {fV_nl <- fV_nl[ ,"mouse"]} 
  if(species=="human") {fV_nl <- fV_nl[ ,"human"]} 

  ##############################################################################
  #DATA:        fV_np (fraction of tissue volume that is neutral phospholipids)                                                                 
  #UNIT:        [fraction]                                                     
  #SOURCE:      rat:   Rodgers and Rowland, J Pharm Sci. 95: 1238-1257 (2006)                                                                                 
  #                    Values have been corrected for residual blood           
  #                    contributions, i.e., reported values do not contain     
  #                    resdiual blood contributions (see Rodgers and Rowland,    
  #                    p.1252 top left paragraph)                              
  #                    adi: in Rodgers and Rowland (2006) incorrectly reported 
  #                         under neutral lipids                               
  #                    pla: p. 1241, paragarph "Tissue Specific Input          
  #                         Parameters"                                        
  #                    ery: Rodgers et al. J Pharm Sci. 94: 1259-1276 (2005); Table 1                                                                             
  #             human: Poulin and Theil, J Pharm Sci. 91: 129-156 (2002)                                                                             
  #                    ery: Rodgers and Rowland, Pharm Res. 24: 918-933 (2007); Table VII                         
  #ASSUMPTIONS: 1) mouse: values are identical to rat values                     
  #             2) fractions of experimental tissue volumes are identical to   
  #                the fractions of tissue volume                              
  ##############################################################################

  fV_np <- matrix(NA,nrow=14,ncol=3)
  dimnames(fV_np) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                            "hea","kid","liv","lun","mus","ski","spl"), 
                          c("rat","mouse","human"))

  ####################   rat     mouse   human 
  fV_np["blo", ] <- c(   NA,     NA,     NA      )
  fV_np["pla", ] <- c(   0.0013, 0.0013, 0.00225 )
  fV_np["ery", ] <- c(   0.0029, 0.0029, 0.0012  )
  fV_np["adi", ] <- c(   0.0016, 0.0016, 0.002   )
  fV_np["bon", ] <- c(   0.0016, 0.0016, 0.0011  )
  fV_np["bra", ] <- c(   0.0015, 0.0015, 0.0565  )
  fV_np["gut", ] <- c(   0.0124, 0.0124, 0.0163  )
  fV_np["hea", ] <- c(   0.0106, 0.0106, 0.0166  )
  fV_np["kid", ] <- c(   0.0240, 0.0240, 0.0162  )
  fV_np["liv", ] <- c(   0.0238, 0.0238, 0.0252  )
  fV_np["lun", ] <- c(   0.0123, 0.0123, 0.009   )
  fV_np["mus", ] <- c(   0.0072, 0.0072, 0.0072  )
  fV_np["ski", ] <- c(   0.0044, 0.0044, 0.0111  )
  fV_np["spl", ] <- c(   0.0107, 0.0107, 0.0198  ) 

  species <- GetSpecies(speciestype)
  if(species=="rat")   {fV_np <- fV_np[ ,"rat"]} 
  if(species=="mouse") {fV_np <- fV_np[ ,"mouse"]} 
  if(species=="human") {fV_np <- fV_np[ ,"human"]} 

  #############################################################################
  #DATA:        fOW_ap (intra-cellular acidic phospholipids in rat)           
  #UNIT:        [mg/g tissue] (original) scaled to fraction                   
  #SOURCE:      rat:   Rodgers et al. J Pharm Sci. 94: 1259-1276 (2005)       
  #                    Not so clear, whether fractions are with respect to    
  #                    residual blood corrected or contaminated tissue weight.
  #                    We assume that values have been corrected for residual 
  #                    blood (as fV_nl or fV_np in the same article)          
  #             human: Poulin and Theil, J Pharm Sci, Vol. 91 (2002)           
  #                    Rodgers and Rowland, J Pharm Res, Vol. 24 (2007)      
  #                    (erythrocyte values only), (Table VII)                       
  #ASSUMPTIONS: 1) mouse: values are identical to rat values                
  #             2) human: values are identical to rat values                  
  #############################################################################

  fOW_ap <- matrix(NA,nrow=14,ncol=3)
  dimnames(fOW_ap) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                             "hea","kid","liv","lun","mus","ski","spl"), 
                           c("rat","mouse","human"))

  ####################   rat     mouse   human 
  fOW_ap["blo", ] <- c(   NA,     NA,     NA    )
  fOW_ap["pla", ] <- c(   NA,     NA,     NA    )
  fOW_ap["ery", ] <- c(   0.5,    0.5,    0.5   )
  fOW_ap["adi", ] <- c(   0.40,   0.40,   0.40  )
  fOW_ap["bon", ] <- c(   0.67,   0.67,   0.67  )
  fOW_ap["bra", ] <- c(   0.40,   0.40,   0.40  )
  fOW_ap["gut", ] <- c(   2.41,   2.41,   2.41  )
  fOW_ap["hea", ] <- c(   2.25,   2.25,   2.25  )
  fOW_ap["kid", ] <- c(   5.03,   5.03,   5.03  )
  fOW_ap["liv", ] <- c(   4.56,   4.56,   4.56  )
  fOW_ap["lun", ] <- c(   3.91,   3.91,   3.91  )
  fOW_ap["mus", ] <- c(   1.53,   1.53,   1.53  )
  fOW_ap["ski", ] <- c(   1.32,   1.32,   1.32  )
  fOW_ap["spl", ] <- c(   3.18,   3.18,   3.18  )

  species <- GetSpecies(speciestype)
  if(species=="rat")   {fOW_ap <- fOW_ap[ ,"rat"]} 
  if(species=="mouse") {fOW_ap <- fOW_ap[ ,"mouse"]} 
  if(species=="human") {fOW_ap <- fOW_ap[ ,"human"]} 
 
  fOW_ap

  #scaling from weight fraction mg/g to weight fraction 10e-3 g/g tissue that 
  #are identified with the volume fraction
  fV_ap <- fOW_ap/1000

  #############################################################################
  #DATA:        r_alb (albumin as tissue-to-plasma ratio in rat)                                                                                
  #UNIT:        []                                                                                                                              
  #SOURCE:      rat:   Rodgers and Rowland, J Pharm Sci. 95: 1238-1257 (2006)     
  #                    Not so clear, whether fractions are with respect to    
  #                    residual blood corrected or contaminated tissue weight.
  #                    We assume that values have been corrected for residual 
  #                    blood (as fV_nl or fV_np in the same article)          
  #             human: Poulin and Theil, J Pharm Sci. 91: 129-156 (2002)              
  #                    ery: Rodgers and Rowland, Pharm Res. 24: 918-933 (2007); Table VII                       
  #ASSUMPTIONS: 1) mouse: values are identical to rat values                   
  #             2) human: values are identical to rat values                  
  #############################################################################

  r_alb <- matrix(NA,nrow=14,ncol=3)
  dimnames(r_alb) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                            "hea","kid","liv","lun","mus","ski","spl"), 
                          c("rat","mouse","human"))

  ####################   rat     mouse   human 
  r_alb["blo", ] <- c(   NA,     NA,     NA    )
  r_alb["pla", ] <- c(   NA,     NA,     NA    )
  r_alb["ery", ] <- c(   NA,     NA,     NA    )
  r_alb["adi", ] <- c(   0.049,  0.049,  0.049 )
  r_alb["bon", ] <- c(   0.100,  0.100,  0.100 )
  r_alb["bra", ] <- c(   0.048,  0.048,  0.048 )
  r_alb["gut", ] <- c(   0.158,  0.158,  0.158 ) 
  r_alb["hea", ] <- c(   0.157,  0.157,  0.157 )
  r_alb["kid", ] <- c(   0.130,  0.130,  0.130 ) 
  r_alb["liv", ] <- c(   0.086,  0.086,  0.086 )
  r_alb["lun", ] <- c(   0.212,  0.212,  0.212 )
  r_alb["mus", ] <- c(   0.064,  0.064,  0.064 ) 
  r_alb["ski", ] <- c(   0.277,  0.277,  0.277 )
  r_alb["spl", ] <- c(   0.097,  0.097,  0.097 )

  species <- GetSpecies(speciestype)
  if(species=="rat")   {r_alb <- r_alb[ ,"rat"]} 
  if(species=="mouse") {r_alb <- r_alb[ ,"mouse"]} 
  if(species=="human") {r_alb <- r_alb[ ,"human"]} 

  #############################################################################
  #DATA:        r_lpr (lipoprotein as tissue-to-plasma ratio in rat)                                                                             
  #UNIT:        []                                                                                                                              
  #SOURCE:      rat:   Rodgers and Rowland, J Pharm Sci. 95: 1238-1257 (2006)              
  #                    Not so clear, whether fractions are with respect to    
  #                    residual blood corrected or contaminated tissue weight.                
  #                    We assume that values have been corrected for residual 
  #                    blood (as fV_nl or fV_np in the same article)          
  #             human: Poulin and Theil, J Pharm Sci. 91: 129-156 (2002)            
  #                    ery: Rodgers and Rowland, Pharm Res. 24: 918-933 (2007); Table VII                                        
  #ASSUMPTIONS: 1) mouse: values are identical to rat values                                                                                     
  #             2) human: values are identical to rat values                                                                                    
  #############################################################################

  r_lpr <- matrix(NA,nrow=14,ncol=3)
  dimnames(r_lpr) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                            "hea","kid","liv","lun","mus","ski","spl"), 
                          c("rat","mouse","human"))

  ####################   rat     mouse   human 
  r_lpr["blo", ] <- c(   NA,     NA,     NA    )
  r_lpr["pla", ] <- c(   NA,     NA,     NA    )
  r_lpr["ery", ] <- c(   NA,     NA,     NA    )
  r_lpr["adi", ] <- c(   0.068,  0.068,  0.068 )
  r_lpr["bon", ] <- c(   0.050,  0.050,  0.050 )
  r_lpr["bra", ] <- c(   0.041,  0.041,  0.041 )
  r_lpr["gut", ] <- c(   0.141,  0.141,  0.141 )
  r_lpr["hea", ] <- c(   0.160,  0.160,  0.160 )
  r_lpr["kid", ] <- c(   0.137,  0.137,  0.137 )
  r_lpr["liv", ] <- c(   0.161,  0.161,  0.161 )
  r_lpr["lun", ] <- c(   0.168,  0.168,  0.168 )
  r_lpr["mus", ] <- c(   0.059,  0.059,  0.059 )
  r_lpr["ski", ] <- c(   0.096,  0.096,  0.096 )
  r_lpr["spl", ] <- c(   0.207,  0.207,  0.207 )

  species <- GetSpecies(speciestype)
  if(species=="rat")   {r_lpr <- r_lpr[ ,"rat"]} 
  if(species=="mouse") {r_lpr <- r_lpr[ ,"mouse"]} 
  if(species=="human") {r_lpr <- r_lpr[ ,"human"]} 

  out <- list(pH=pH,fV_w_tot=fV_w_tot, fV_w_ex=fV_w_ex, fV_w_ic=fV_w_ic, 
              fV_nl=fV_nl, fV_np=fV_np, fV_ap=fV_ap, r_alb=r_alb, r_lpr=r_lpr)

  return(out)
}


#--- SPECIES-SPECIFIC DATA SCHMITT --------------------------------------------

GetSpeciesDataSchmitt <- function (speciestype) {
  
  #############################################################################
  #DATA:        pH (pH values in plasma/interstitital water and intra-        
  #             cellular water)                                                 
  #UNIT:        [-]                                                                                                                            
  #SOURCE:      Schmitt Corrigendum Toxicol. Vitro 22: 1666 (2008)
  #             pla: Rodgers et al. J Pharm Sci. 94: 1259-1276 (2005); p. 1263 
  #NOTE:        Values taken from (Waddell and Bates, 1969; Malan et al.,     
  #             1985; Wood and Schaefer, 1978; Schanker and Less, 1977;       
  #             Harrison and Walker, 1979 and Civelelek et al., 1996).        
  #             Mean values were calculated when more than one value was      
  #             found for the same tissue.                                    
  #ASSUMPTION:  pH values of human, mouse and rat are assumed to be identical 
  #############################################################################
  
  pH <- matrix(NA,nrow=14,ncol=3)
  dimnames(pH) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                         "hea","kid","liv","lun","mus","ski","spl"), 
                       c("rat","mouse","human"))     
  
  #################   rat     mouse   human 
  pH["blo", ] <- c(   NA,     NA,     NA    )
  pH["pla", ] <- c(   7.4,    7.4,    7.4   )
  pH["ery", ] <- c(   7.20,   7.20,   7.20  )
  pH["adi", ] <- c(   7.10,   7.10,   7.10  )
  pH["bon", ] <- c(   7.00,   7.00,   7.00  )
  pH["bra", ] <- c(   7.10,   7.10,   7.10  )
  pH["gut", ] <- c(   7.00,   7.00,   7.00  )
  pH["hea", ] <- c(   7.10,   7.10,   7.10  )
  pH["kid", ] <- c(   7.22,   7.22,   7.22  )
  pH["liv", ] <- c(   7.23,   7.23,   7.23  )
  pH["lun", ] <- c(   6.60,   6.60,   6.60  )
  pH["mus", ] <- c(   6.81,   6.81,   6.81  )
  pH["ski", ] <- c(   7.00,   7.00,   7.00  )
  pH["spl", ] <- c(   7.00,   7.00,   7.00  )
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {pH <- pH[ ,"rat"]} 
  if(species=="mouse") {pH <- pH[ ,"mouse"]} 
  if(species=="human") {pH <- pH[ ,"human"]} 
  
  ##################################################################################
  #DATA:        fL_nl (fraction of total lipid that is neutral lipid in rat)       
  #             fL_np (fraction of total lipid that is neutral phospholipid in rat)
  #             fL_ap (fraction of total lipid that is acidic phospholipid in rat) 
  #             fV_lip_cel (fraction of cell volume that is lipids in human)       
  #UNIT:        [fraction]                                                         
  #SOURCE:      Schmitt Corrigendum Toxicol. Vitro 22: 1666 (2008)                                          
  #ASSUMPTIONS: 1) fL_nl, fL_np, fL_ap: mouse values and human values are          
  #                identical to reported rat values                                  
  #             2) fV_lip_cel: rat values and mouse values are identical to        
  #                reported human values                                           
  ##################################################################################

  fL_nl <- matrix(NA,nrow=14,ncol=3)
  dimnames(fL_nl) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                            "hea","kid","liv","lun","mus","ski","spl"), 
                          c("rat","mouse","human"))
  
  ####################   rat     mouse   human 
  fL_nl["blo", ] <- c(   NA,     NA,     NA      )
  fL_nl["pla", ] <- c(   NA,     NA,     NA      )
  fL_nl["ery", ] <- c(   0.30,   0.30,   0.30    )
  fL_nl["adi", ] <- c(   1.00,   1.00,   1.00    )
  fL_nl["bon", ] <- c(   0.85,   0.85,   0.85    )
  fL_nl["bra", ] <- c(   0.39,   0.39,   0.39    )
  fL_nl["gut", ] <- c(   0.69,   0.69,   0.69    )
  fL_nl["hea", ] <- c(   0.48,   0.48,   0.48    )
  fL_nl["kid", ] <- c(   0.26,   0.26,   0.26    )
  fL_nl["liv", ] <- c(   0.29,   0.29,   0.29    )
  fL_nl["lun", ] <- c(   0.51,   0.51,   0.51    )
  fL_nl["mus", ] <- c(   0.49,   0.49,   0.49    )
  fL_nl["ski", ] <- c(   0.90,   0.90,   0.90    )
  fL_nl["spl", ] <- c(   0.30,   0.30,   0.30    )
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {fL_nl <- fL_nl[ ,"rat"]} 
  if(species=="mouse") {fL_nl <- fL_nl[ ,"mouse"]} 
  if(species=="human") {fL_nl <- fL_nl[ ,"human"]} 
  
  fL_np <- matrix(NA,nrow=14,ncol=3)
  dimnames(fL_np) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                            "hea","kid","liv","lun","mus","ski","spl"), 
                          c("rat","mouse","human"))
  
  ####################   rat     mouse   human 
  fL_np["blo", ] <- c(   NA,     NA,     NA      )
  fL_np["pla", ] <- c(   NA,     NA,     NA      )
  fL_np["ery", ] <- c(   0.59,   0.59,   0.59    )
  fL_np["adi", ] <- c(   0.0022, 0.0022, 0.0022  )
  fL_np["bon", ] <- c(   0.11,   0.11,   0.11    )
  fL_np["bra", ] <- c(   0.48,   0.48,   0.48    )
  fL_np["gut", ] <- c(   0.26,   0.26,   0.26    )
  fL_np["hea", ] <- c(   0.43,   0.43,   0.43    )
  fL_np["kid", ] <- c(   0.61,   0.61,   0.61    )
  fL_np["liv", ] <- c(   0.59,   0.59,   0.59    )
  fL_np["lun", ] <- c(   0.38,   0.38,   0.38    )
  fL_np["mus", ] <- c(   0.42,   0.42,   0.42    )
  fL_np["ski", ] <- c(   0.08,   0.08,   0.08    )
  fL_np["spl", ] <- c(   0.54,   0.54,   0.54    )

  species <- GetSpecies(speciestype)
  if(species=="rat")   {fL_np <- fL_np[ ,"rat"]} 
  if(species=="mouse") {fL_np <- fL_np[ ,"mouse"]} 
  if(species=="human") {fL_np <- fL_np[ ,"human"]} 
  
  fL_ap <- matrix(NA,nrow=14,ncol=3)
  dimnames(fL_ap) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                            "hea","kid","liv","lun","mus","ski","spl"), 
                          c("rat","mouse","human"))
  
  ####################   rat     mouse   human 
  fL_ap["blo", ] <- c(   NA,     NA,     NA      )
  fL_ap["pla", ] <- c(   NA,     NA,     NA      )
  fL_ap["ery", ] <- c(   0.10,   0.10,   0.10    )
  fL_ap["adi", ] <- c(   0.0006, 0.0006, 0.0006  )
  fL_ap["bon", ] <- c(   0.04,   0.04,   0.04    )
  fL_ap["bra", ] <- c(   0.13,   0.13,   0.13    )
  fL_ap["gut", ] <- c(   0.05,   0.05,   0.05    )
  fL_ap["hea", ] <- c(   0.09,   0.09,   0.09    )
  fL_ap["kid", ] <- c(   0.13,   0.13,   0.13    )
  fL_ap["liv", ] <- c(   0.11,   0.11,   0.11    )
  fL_ap["lun", ] <- c(   0.11,   0.11,   0.11    )
  fL_ap["mus", ] <- c(   0.09,   0.09,   0.09    )
  fL_ap["ski", ] <- c(   0.02,   0.02,   0.02    )
  fL_ap["spl", ] <- c(   0.15,   0.15,   0.15    )
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {fL_ap <- fL_ap[ ,"rat"]} 
  if(species=="mouse") {fL_ap <- fL_ap[ ,"mouse"]} 
  if(species=="human") {fL_ap <- fL_ap[ ,"human"]} 
  
  #############################################################################
  #DATA:        fV_tot_lip_cel (fraction of cell volume that is total lipids  
  #             in human)                                                     
  #UNIT:        [fraction]                                                    
  #SOURCE:      Schmitt Corrigendum Toxicol. Vitro 22: 1666 (2008)                                   
  #NOTE:        Values from taken from (ICRP, 1975). Original values given as 
  #             fraction of total tissue mass were rescaled to cellular       
  #             volume as follows: Water fraction of total tissue reduced by  
  #             interstitial volume and subsequently                          
  #             all values normalized by cellular fraction.                   
  #ASSUMPTIONS: 1) rat values and mouse values are identical to reported       
  #                human values                                               
  #############################################################################
  
  fV_tot_lip_cel <- matrix(NA,nrow=14,ncol=3)
  dimnames(fV_tot_lip_cel) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                                     "hea","kid","liv","lun","mus","ski","spl"), 
                                   c("rat","mouse","human"))
  
  #############################   rat     mouse   human 
  fV_tot_lip_cel["blo", ] <- c(   NA,     NA,     NA      )
  fV_tot_lip_cel["pla", ] <- c(   NA,     NA,     NA      )
  fV_tot_lip_cel["ery", ] <- c(   0.01,   0.01,   0.01    )
  fV_tot_lip_cel["adi", ] <- c(   0.92,   0.92,   0.92    )
  fV_tot_lip_cel["bon", ] <- c(   0.02,   0.02,   0.02    )
  fV_tot_lip_cel["bra", ] <- c(   0.11,   0.11,   0.11    )
  fV_tot_lip_cel["gut", ] <- c(   0.07,   0.07,   0.07    )
  fV_tot_lip_cel["hea", ] <- c(   0.11,   0.11,   0.11    )
  fV_tot_lip_cel["kid", ] <- c(   0.06,   0.06,   0.06    )
  fV_tot_lip_cel["liv", ] <- c(   0.08,   0.08,   0.08    )
  fV_tot_lip_cel["lun", ] <- c(   0.04,   0.04,   0.04    )
  fV_tot_lip_cel["mus", ] <- c(   0.01,   0.01,   0.01    )
  fV_tot_lip_cel["ski", ] <- c(   0.14,   0.14,   0.14    )
  fV_tot_lip_cel["spl", ] <- c(   0.02,   0.02,   0.02    )
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {fV_tot_lip_cel <- fV_tot_lip_cel[ ,"rat"]} 
  if(species=="mouse") {fV_tot_lip_cel <- fV_tot_lip_cel[ ,"mouse"]} 
  if(species=="human") {fV_tot_lip_cel <- fV_tot_lip_cel[ ,"human"]} 
  
  #calculate fraction of neutral lipids in cell
  fV_nl_cel <- fV_tot_lip_cel*fL_nl
  
  #calculate fraction of neutral phopspholipids in cell
  fV_np_cel <- fV_tot_lip_cel*fL_np
  
  #calculate fraction of acidic phopspholipids in cell
  fV_ap_cel <- fV_tot_lip_cel*fL_ap
  
  #############################################################################
  #DATA:        fV_prot_cel (fraction of cell volume that is protein in human)                                                                  
  #UNIT:        [fraction]                                                    
  #SOURCE:      Schmitt Corrigendum Toxicol. Vitro 22: 1666 (2008)                                     
  #NOTE:        Values from taken from (ICRP, 1975). Original values given as 
  #             fraction of total tissue mass were rescaled to                
  #             cellular volume as follows: Water fraction of total tissue    
  #             reduced by interstitial volume and subsequently               
  #             all values normalized by cellular fraction.                   
  #ASSUMPTIONS: 1) rat values and mouse values are identical to reported      
  #                human values                                               
  #############################################################################

  fV_prot_cel <- matrix(NA,nrow=14,ncol=3)
  dimnames(fV_prot_cel) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                                  "hea","kid","liv","lun","mus","ski","spl"), 
                                c("rat","mouse","human"))
  
  ##########################   rat     mouse   human 
  fV_prot_cel["blo", ] <- c(   NA,     NA,     NA      )
  fV_prot_cel["pla", ] <- c(   NA,     NA,     NA      )
  fV_prot_cel["ery", ] <- c(   0.33,   0.33,   0.33    )
  fV_prot_cel["adi", ] <- c(   0.06,   0.06,   0.06    )
  fV_prot_cel["bon", ] <- c(   0.21,   0.21,   0.21    )
  fV_prot_cel["bra", ] <- c(   0.08,   0.08,   0.08    )
  fV_prot_cel["gut", ] <- c(   0.15,   0.15,   0.15    )
  fV_prot_cel["hea", ] <- c(   0.19,   0.19,   0.19    )
  fV_prot_cel["kid", ] <- c(   0.21,   0.21,   0.21    )
  fV_prot_cel["liv", ] <- c(   0.21,   0.21,   0.21    )
  fV_prot_cel["lun", ] <- c(   0.11,   0.11,   0.11    )
  fV_prot_cel["mus", ] <- c(   0.19,   0.19,   0.19    )
  fV_prot_cel["ski", ] <- c(   0.23,   0.23,   0.23    )
  fV_prot_cel["spl", ] <- c(   0.30,   0.30,   0.30    )
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {fV_prot_cel <- fV_prot_cel[ ,"rat"]} 
  if(species=="mouse") {fV_prot_cel <- fV_prot_cel[ ,"mouse"]} 
  if(species=="human") {fV_prot_cel <- fV_prot_cel[ ,"human"]} 
  
  #############################################################################
  #DATA:        fV_w_ic (fraction of tissue volume that is cellular water in  
  #             human)                                                        
  #UNIT:        [fraction]                                                    
  #SOURCE:      Schmitt Corrigendum Toxicol. Vitro 22: 1666 (2008)                                     
  #NOTE:        Values from taken from (ICRP, 1975). Original values given as 
  #             fraction of total tissue mass were rescaled to cellular       
  #             volume as follows: Water fraction of total tissue reduced by  
  #             interstitial volume and subsequently all values normalized    
  #             by cellular fraction.                                         
  #ASSUMPTIONS: 1) rat values and mouse values are identical to reported      
  #                human values                                               
  #############################################################################

  fV_w_ic <- matrix(NA,nrow=14,ncol=3)
  dimnames(fV_w_ic) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                              "hea","kid","liv","lun","mus","ski","spl"), 
                            c("rat","mouse","human"))
  
  ######################   rat     mouse   human 
  fV_w_ic["blo", ] <- c(   NA,     NA,     NA      )
  fV_w_ic["pla", ] <- c(   NA,     NA,     NA      )
  fV_w_ic["ery", ] <- c(   0.63,   0.63,   0.63    )
  fV_w_ic["adi", ] <- c(   0.03,   0.03,   0.03    )
  fV_w_ic["bon", ] <- c(   0.26,   0.26,   0.26    )
  fV_w_ic["bra", ] <- c(   0.79,   0.79,   0.79    )
  fV_w_ic["gut", ] <- c(   0.78,   0.78,   0.78    )
  fV_w_ic["hea", ] <- c(   0.70,   0.70,   0.70    )
  fV_w_ic["kid", ] <- c(   0.73,   0.73,   0.73    )
  fV_w_ic["liv", ] <- c(   0.68,   0.68,   0.68    )
  fV_w_ic["lun", ] <- c(   0.74,   0.74,   0.74    )
  fV_w_ic["mus", ] <- c(   0.76,   0.76,   0.76    )
  fV_w_ic["ski", ] <- c(   0.47,   0.47,   0.47    )
  fV_w_ic["spl", ] <- c(   0.75,   0.75,   0.75    )
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {fV_w_ic <- fV_w_ic[ ,"rat"]} 
  if(species=="mouse") {fV_w_ic <- fV_w_ic[ ,"mouse"]} 
  if(species=="human") {fV_w_ic <- fV_w_ic[ ,"human"]} 
  
  #############################################################################
  #DATA:        fV_cel (fraction of tissue volume that is cell volume in rat) 
  #UNIT:        [fraction]                                                    
  #SOURCE:      Schmitt Corrigendum Toxicol. Vitro 22: 1666 (2008)                                     
  #NOTE:        Values from taken from (Kawai et al., 1994). Original values  
  #             given as fraction of total organ volume were rescaled to      
  #             tissue volume by subtracting vascular volume.                 
  #ASSUMPTIONS: 1) mouse values and human values are identical to reported    
  #                rat values                                                 
  #############################################################################
  
  fV_cel <- matrix(NA,nrow=14,ncol=3)
  dimnames(fV_cel) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                             "hea","kid","liv","lun","mus","ski","spl"), 
                           c("rat","mouse","human"))
  
  #####################   rat     mouse   human 
  fV_cel["blo", ] <- c(   NA,     NA,     NA      )
  fV_cel["pla", ] <- c(   NA,     NA,     NA      )
  fV_cel["ery", ] <- c(   1.00,   1.00,   1.00    )
  fV_cel["adi", ] <- c(   0.86,   0.86,   0.86    )
  fV_cel["bon", ] <- c(   0.90,   0.90,   0.90    )
  fV_cel["bra", ] <- c(   1.00,   1.00,   1.00    )
  fV_cel["gut", ] <- c(   0.90,   0.90,   0.90    )
  fV_cel["hea", ] <- c(   0.86,   0.86,   0.86    )
  fV_cel["kid", ] <- c(   0.78,   0.78,   0.78    )
  fV_cel["liv", ] <- c(   0.82,   0.82,   0.82    )
  fV_cel["lun", ] <- c(   0.50,   0.50,   0.50    )
  fV_cel["mus", ] <- c(   0.88,   0.88,   0.88    )
  fV_cel["ski", ] <- c(   0.69,   0.69,   0.69    )
  fV_cel["spl", ] <- c(   0.79,   0.79,   0.79    )
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {fV_cel <- fV_cel[ ,"rat"]} 
  if(species=="mouse") {fV_cel <- fV_cel[ ,"mouse"]} 
  if(species=="human") {fV_cel <- fV_cel[ ,"human"]} 
  
  #############################################################################
  #DATA:        fV_int (fraction of tissue volume that is interstitial volume 
  #             in rat)                                                       
  #UNIT:        [fraction]                                                    
  #SOURCE:      Schmitt Corrigendum Toxicol. Vitro 22: 1666 (2008)                                    
  #NOTE:        Values from taken from (Kawai et al., 1994). Original values  
  #             given as fraction of total organ volume were rescaled to      
  #             tissue volume by subtracting vascular volume.                 
  #ASSUMPTIONS: 1) mouse values and human values are identical to reported    
  #                rat values                                                 
  #############################################################################
  
  fV_int <- matrix(NA,nrow=14,ncol=3)
  dimnames(fV_int) <- list(c("blo","pla","ery","adi","bon","bra","gut",
                             "hea","kid","liv","lun","mus","ski","spl"), 
                           c("rat","mouse","human"))
  
  #####################   rat     mouse   human 
  fV_int["blo", ] <- c(   NA,     NA,     NA      )
  fV_int["pla", ] <- c(   NA,     NA,     NA      )
  fV_int["ery", ] <- c(   NA,     NA,     NA      )
  fV_int["adi", ] <- c(   0.14,   0.14,   0.14    )
  fV_int["bon", ] <- c(   0.10,   0.10,   0.10    )
  fV_int["bra", ] <- c(   0.004,  0.004,  0.004   )
  fV_int["gut", ] <- c(   0.096,  0.096,  0.096   )
  fV_int["hea", ] <- c(   0.14,   0.14,   0.14    )
  fV_int["kid", ] <- c(   0.22,   0.22,   0.22    )
  fV_int["liv", ] <- c(   0.18,   0.18,   0.18    )
  fV_int["lun", ] <- c(   0.50,   0.50,   0.50    )
  fV_int["mus", ] <- c(   0.12,   0.12,   0.12    )
  fV_int["ski", ] <- c(   0.31,   0.31,   0.31    )
  fV_int["spl", ] <- c(   0.21,   0.21,   0.21    )
  
  species <- GetSpecies(speciestype)
  if(species=="rat")   {fV_int <- fV_int[ ,"rat"]} 
  if(species=="mouse") {fV_int <- fV_int[ ,"mouse"]} 
  if(species=="human") {fV_int <- fV_int[ ,"human"]} 
  
  out <- list(pH=pH, fV_tot_lip_cel=fV_tot_lip_cel, fV_nl_cel=fV_nl_cel, fV_np_cel=fV_np_cel, 
              fV_ap_cel=fV_ap_cel, fV_prot_cel=fV_prot_cel, 
              fV_w_ic=fV_w_ic, fV_cel=fV_cel, fV_int=fV_int)
  
  return (out)
}



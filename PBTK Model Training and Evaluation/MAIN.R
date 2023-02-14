###############################################################################
# This main script 
# - processes observed data of retrorsine                  
# - initializes and performs MCMC parameter estimation 
# - performs model prediction with optimized parameters
# - plots training and evaluation data and simulations
# Authors:  Anja Lehmann, Christoph Hethey                                                       
# Date: 2023-02-11 
#
# The functions GetSpeciesData() and UseMethodRodgersAndRowland() 
# are based on a MATBLAB script originally written by                             
# W. Huisinga, A. Solms, L. Fronton, S. Pilari, Modeling Interindividual       
# Variability in Physiologically Based Pharmacokinetics and Its Link to       
# Mechanistic Covariate Modeling, CPT: Pharmacometrics & Systems Pharmacology 
# (2012) 1, e5; doi:10.1038/psp.2012.3                                        
###############################################################################

#--- INITIALIZATION 

# Remove all objects from current workspace
rm(list=ls())

# Reset graphics
graphics.off()
cat("\014")

# Set working directory to current path
PATH <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(PATH)

# Load dependencies
source("MyFunctions.R")
source("DataProcessing.R")
source("DrugDatabase.R")
source("SpeciesDatabase.R")
source("ParameterVector.R")
source("PartitionCoefficients.R")
source("Pred.R")
source("Minus2LL.R")

# Load packages
library(readr)
library(crayon)
library(RxODE)
library(deSolve)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(dfoptim)
library(neldermead)
library(cowplot)
library(profvis)
library(FME)
library(ucminf)
library(nloptr)
library(matlib)
library(bayestestR)

# Define color-blind friendly palette (Okate and Ito 2002)
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", 
                      "#CC79A7", "#999999", "#F0E442") 

###############################################################################
# PART: OBSERVED DATA                                                       
###############################################################################

#read dataset
dataset_retrorsine <- read_delim("ObservedData.csv",";", 
                                 escape_double = FALSE, 
                                 col_types = cols(DOSEORIGINAL = col_number(), 
                                                  XORIGINAL = col_number(),
                                                  WEIGHT = col_number(), 
                                                  YORIGINAL = col_number()),
                                 trim_ws = TRUE)

# Process data
dataset_retrorsine <- GetDataProcessing(dataset=dataset_retrorsine,deagg=TRUE) 
View(dataset_retrorsine)

#--- LOG TRANSFORM OBSERVED DATA 
dataset_retrorsine$Y <- TransData(dataset_retrorsine$Y)

###############################################################################
# PART: INITIALIZATION OF PREDICTION                                              
###############################################################################

#--- DEFINE SPECIESTYPES

# speciestype is a string defining species, sex (only human 15 and 35 years), 
# and age, e.g. "rat10weeks"
# Speciestypes, for which data is available in SpeciesData.R: 
# "rat10weeks",   "rat70weeks", 
# "mouse9weeks",  "mouse70weeks",
# "humannewborn", "human1year", "human5years",  "human10years", 
# "humanm15years","humanf15years","humanm35years","humanf35years"

#Extract speciestypes
speciestypes <- as.matrix(unique(subset(dataset_retrorsine,select = SPECIESTYPE)))

#--- LOAD THETA FOR EACH SPECIESTYPE and PRINT PLAUSIBILITY CHECKS 

for(i in seq_along(speciestypes)) {
  speciestype       <- speciestypes[i,]
  drug <- "Retrorsine"
  
  theta <- GetParameterVector(drug,speciestype,
                              method_partcoeff="RodgersAndRowland")
  # "Retrorsine" has to be included in DrugDatabase.R
  # method_partcoeff choose between "RodgersAndRowland" / 
  # "Schmitt", NOTE: also change in Pred.R
  
  #Assign variable names xxx_speciestype
  assign(paste("theta",speciestype,sep="_"), theta)
}

#--- CREATE RUNTABLE AND RUNID 

# Create runtable and give runid for each run
runtable_retrorsine <- data.frame(unique(subset(dataset_retrorsine,!(is.na(dataset_retrorsine$DOSE)),
                                     select = c("DOSE","DOSEUNIT","SPECIESTYPE","ADM"))))
runtable_retrorsine$RUNID <- seq.int(nrow(runtable_retrorsine))
View(runtable_retrorsine)
  
# Create subdataset without NAs in DOSE column
subdataset_retrorsine <- subset(dataset_retrorsine, !(is.na(dataset_retrorsine$DOSE)))
  
# Initialize new column RUNID in subdataset and dataset
subdataset_retrorsine$RUNID <- as.numeric(NA)
dataset_retrorsine$RUNID    <- as.numeric(NA)
  
# Loop through all rows of subdataset to identify which rows (of columns DOSE and SPECIESTYPE)
# are identical to runtable columns DOSE AND SPECIESTYPE and allocate runid in column RUNID
for(row in 1:nrow(runtable_retrorsine)) {
    speciestype <- runtable_retrorsine[[row,"SPECIESTYPE"]]
    dose        <- runtable_retrorsine[[row,"DOSE"]]
    adm         <- runtable_retrorsine[[row,"ADM"]]
    runid       <- runtable_retrorsine[[row,"RUNID"]]
    for(row in 1:nrow(subdataset_retrorsine)) {
      if((subdataset_retrorsine[[row,"DOSE"]]==dose) & 
         (subdataset_retrorsine[[row,"SPECIESTYPE"]]==speciestype) &
         (subdataset_retrorsine[[row,"ADM"]]==adm) ) 
      {subdataset_retrorsine[[row,"RUNID"]] <- runid}
    }
  }  
  
# Transfer runid from subdataset to original dataset
for(row in 1:nrow(subdataset_retrorsine)) {
    id    <- subdataset_retrorsine[[row, "ID"]]
    runid <- subdataset_retrorsine[[row,"RUNID"]] 
    for(row in 1:nrow(dataset_retrorsine)) {
      if(dataset_retrorsine[[row,"ID"]] == id) 
      {dataset_retrorsine[[row,"RUNID"]] <- runid}
    }
  }

#--- DEFINE ODE SYSTEM FOR RxODE PACKAGE 

mRetrorsine <-RxODE({
d/dt (Alum)   = -ka*Fa*Alum;
d/dt (Aper)    = -kip*Aper;
d/dt (Aven)   = Qadi*Aadi/Vadi/Kadi+Qbon*Abon/Vbon/Kbon+Qbra*Abra/Vbra/Kbra+Qhea*Ahea/Vhea/Khea+Qkid*Akid/Vkid/Kkid+Qliv*Kvasvi*Alivvi/Vlivvi+Qmus*Amus/Vmus/Kmus+Qski*Aski/Vski/Kski-Qc*Aven/Vven;
d/dt (Aart)   = Qc*Alun/Vlun/Klun-Qc*Aart/Vart;
d/dt (Aadi)   = Qadi*Aart/Vart-Qadi*Aadi/Vadi/Kadi;
d/dt (Abon)   = Qbon*Aart/Vart-Qbon*Abon/Vbon/Kbon;
d/dt (Abra)   = Qbra*Aart/Vart-Qbra*Abra/Vbra/Kbra;
d/dt (Agut)   = Qgut*Aart/Vart-Qgut*Agut/Vgut/Kgut-((Vmaxliv/10)/(Km+fuP*(Agut/Vgut)))*fuP*Agut/Vgut+ka*Fa*Alum;
d/dt (Ahea)   = Qhea*Aart/Vart-Qhea*Ahea/Vhea/Khea;
d/dt (Akid)   = Qkid*Aart/Vart-Qkid*Akid/Vkid/Kkid-fuP*GFR*Akid/Vkid;
d/dt (Alivvi) = kip*Aper+(Qliv-Qspl-Qgut)*Aart/Vart+Qspl*Aspl/Vspl/Kspl+Qgut*Agut/Vgut/Kgut-Qliv*Kvasvi*Alivvi/Vlivvi+CLactef*fuliv*Alivc/Vlivc+PSdiff*fnliv/fnpla*fuliv*Alivc/Vlivc-CLactin*Alivvi/Vlivvi*Kintuvi-PSdiff*Alivvi/Vlivvi*Kintuvi;
d/dt (Alivc)  = CLactin*Kintuvi*Alivvi/Vlivvi+PSdiff*Kintuvi*Alivvi/Vlivvi-CLactef*fuliv*Alivc/Vlivc-PSdiff*fnliv/fnpla*fuliv*Alivc/Vlivc-CLcanef*fuliv*Alivc/Vlivc-(fDHRPROT+fDHRDNA*10^-6+fDHRGSH*10^-4+(fsum-fDHRGSH*10^-4-fDHRPROT-fDHRDNA*10^-6))*(Vmaxliv/(Km+fuliv*(Alivc/Vlivc)))*fuliv*Alivc/Vlivc;
d/dt (Alun)   = Qc*Aven/Vven-Qc*Alun/Vlun/Klun;
d/dt (Amus)   = Qmus*Aart/Vart-Qmus*Amus/Vmus/Kmus;
d/dt (Aski)   = Qski*Aart/Vart-Qski*Aski/Vski/Kski;
d/dt (Aspl)   = Qspl*Aart/Vart-Qspl*Aspl/Vspl/Kspl;
d/dt (DHRGSH) = fDHRGSH*10^-4*(Vmaxliv/(Km+fuliv*(Alivc/Vlivc)))*fuliv*Alivc/Vlivc-lambda1DHRGSH*DHRGSH;
d/dt (DHRPROT) = fDHRPROT*(Vmaxliv/(Km+fuliv*(Alivc/Vlivc)))*fuliv*Alivc/Vlivc-lambda1DHRPROT*DHRPROT;
d/dt (DHRDNA) = fDHRDNA*10^-6*(Vmaxliv/(Km+fuliv*(Alivc/Vlivc)))*fuliv*Alivc/Vlivc-(lambda2DHRDNA*10^-3+(lambda1DHRDNA*10^-2-lambda2DHRDNA*10^-3)*exp(-kDHRDNA*t))*DHRDNA;
d/dt (Auri)   = fuP*GFR*Akid/Vkid;
d/dt (Abil)   = CLcanef*fuliv*Alivc/Vlivc;
d/dt (bag)    = (fsum-fDHRGSH*10^-4-fDHRPROT-fDHRDNA*10^-6)*(Vmaxliv/(Km+fuliv*(Alivc/Vlivc)))*fuliv*Alivc/Vlivc+lambda1DHRGSH*DHRGSH+lambda1DHRPROT*DHRPROT +(lambda2DHRDNA*10^-3+(lambda1DHRDNA*10^-2-lambda2DHRDNA*10^-3)*exp(-kDHRDNA*t))*DHRDNA +((Vmaxliv/10)/(Km+fuP*(Agut/Vgut)))*fuP*Agut/Vgut;

fsum = fDHRGSH*10^-4+fDHRPROT+fDHRDNA*10^-6+(fsum-fDHRGSH*10^-4-fDHRPROT-fDHRDNA*10^-6);

})


#--- ODE to derive fraction of metabolites formed; i.e. no depletion functions for metabolites
mRetrorsine2 <-RxODE({
  d/dt (Alum)   = -ka*Fa*Alum;
  d/dt (Aper)    = -kip*Aper;
  d/dt (Aven)   = Qadi*Aadi/Vadi/Kadi+Qbon*Abon/Vbon/Kbon+Qbra*Abra/Vbra/Kbra+Qhea*Ahea/Vhea/Khea+Qkid*Akid/Vkid/Kkid+Qliv*Kvasvi*Alivvi/Vlivvi+Qmus*Amus/Vmus/Kmus+Qski*Aski/Vski/Kski-Qc*Aven/Vven;
  d/dt (Aart)   = Qc*Alun/Vlun/Klun-Qc*Aart/Vart;
  d/dt (Aadi)   = Qadi*Aart/Vart-Qadi*Aadi/Vadi/Kadi;
  d/dt (Abon)   = Qbon*Aart/Vart-Qbon*Abon/Vbon/Kbon;
  d/dt (Abra)   = Qbra*Aart/Vart-Qbra*Abra/Vbra/Kbra;
  d/dt (Agut)   = Qgut*Aart/Vart-Qgut*Agut/Vgut/Kgut-((Vmaxliv/10)/(Km+fuP*(Agut/Vgut)))*fuP*Agut/Vgut+ka*Fa*Alum;
  d/dt (Ahea)   = Qhea*Aart/Vart-Qhea*Ahea/Vhea/Khea;
  d/dt (Akid)   = Qkid*Aart/Vart-Qkid*Akid/Vkid/Kkid-fuP*GFR*Akid/Vkid;
  d/dt (Alivvi) = kip*Aper+(Qliv-Qspl-Qgut)*Aart/Vart+Qspl*Aspl/Vspl/Kspl+Qgut*Agut/Vgut/Kgut-Qliv*Kvasvi*Alivvi/Vlivvi+CLactef*fuliv*Alivc/Vlivc+PSdiff*fnliv/fnpla*fuliv*Alivc/Vlivc-CLactin*Alivvi/Vlivvi*Kintuvi-PSdiff*Alivvi/Vlivvi*Kintuvi;
  d/dt (Alivc)  = CLactin*Kintuvi*Alivvi/Vlivvi+PSdiff*Kintuvi*Alivvi/Vlivvi-CLactef*fuliv*Alivc/Vlivc-PSdiff*fnliv/fnpla*fuliv*Alivc/Vlivc-CLcanef*fuliv*Alivc/Vlivc-(fDHRPROT+fDHRDNA*10^-6+fDHRGSH*10^-4+(fsum-fDHRGSH*10^-4-fDHRPROT-fDHRDNA*10^-6))*(Vmaxliv/(Km+fuliv*(Alivc/Vlivc)))*fuliv*Alivc/Vlivc;
  d/dt (Alun)   = Qc*Aven/Vven-Qc*Alun/Vlun/Klun;
  d/dt (Amus)   = Qmus*Aart/Vart-Qmus*Amus/Vmus/Kmus;
  d/dt (Aski)   = Qski*Aart/Vart-Qski*Aski/Vski/Kski;
  d/dt (Aspl)   = Qspl*Aart/Vart-Qspl*Aspl/Vspl/Kspl;
  d/dt (DHRGSH) = fDHRGSH*10^-4*(Vmaxliv/(Km+fuliv*(Alivc/Vlivc)))*fuliv*Alivc/Vlivc;
  d/dt (DHRPROT) = fDHRPROT*(Vmaxliv/(Km+fuliv*(Alivc/Vlivc)))*fuliv*Alivc/Vlivc;
  d/dt (DHRDNA) = fDHRDNA*10^-6*(Vmaxliv/(Km+fuliv*(Alivc/Vlivc)))*fuliv*Alivc/Vlivc;
  d/dt (Auri)   = fuP*GFR*Akid/Vkid;
  d/dt (Abil)   = CLcanef*fuliv*Alivc/Vlivc;
  d/dt (AMetgut) = ((Vmaxliv/10)/(Km+fuP*(Agut/Vgut)))*fuP*Agut/Vgut;
  d/dt (AMetlivother) = (fsum-fDHRGSH*10^-4-fDHRPROT-fDHRDNA*10^-6)*(Vmaxliv/(Km+fuliv*(Alivc/Vlivc)))*fuliv*Alivc/Vlivc;

  fsum = fDHRGSH*10^-4+fDHRPROT+fDHRDNA*10^-6+(fsum-fDHRGSH*10^-4-fDHRPROT-fDHRDNA*10^-6);

})


#################################################################################
# PART: MAXIMUM LIKELIHOOD ESTIMATION                                       
#################################################################################

################################################################################# 
# CAUTION: There is a bug in function nmkb of package dfoptim. To solve the bug: 
# the position of the elements in vector 'est' have to be assigned to 
# the names of the elements of 'est' in functions GetMinus2LL and GetPred 
#################################################################################

# Extract all timepoints of observed data from dataset
timepoints <- unique(subset(subset(dataset_retrorsine,
                                   is.na(dataset_retrorsine[["DOSE"]])), 
                            select = X))
timepoints <- sort(as.vector(timepoints$X))

#--- DEFINE START VALUES FOR PARAMETERS TO BE ESTIMATED 

  est0 <- c(a_RTS_m_pla=      0.1,
            a_RTS_r_bil=      0.1,
            a_RTS_m_uri=      0.1,
            a_RTS_r_uri=      0.1,
            a_DHRGSH_m_liv=   0.1,
            a_DHRPROT_m_liv=  0.1,
            a_DHRDNA_m_liv=   1, 
            lambda1DHRGSH=    5,
            lambda1DHRPROT=   0.1,
            lambda1DHRDNA=    5,     
            lambda2DHRDNA=    2,
            kDHRDNA=          0.06,
            fDHRGSH=          0.5, 
            fDHRPROT=         0.002,
            fDHRDNA=          1.3, 
            CLcanefr=         0.5,
            ka =              0.8
            )
            
#--- DEFINE BOUNDARIES FOR PARAMETERS TO BE ESTIMATED 
lboundary <- c(1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,
               0, 
               0, 
               0, 
               0, 
               0, 
               0, 
               0, 
               0, 
               0, 
               0) 
uboundary <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf, 
               Inf, 
               Inf, 
               Inf, 
               Inf, 
               Inf, 
               Inf, 
               Inf, 
               Inf, 
               Inf, 
               Inf) 

#--- Hooke Jeeves direct search from Package dfoptim

# set.seed(10)
# 
# out_hjkb <- hjkb(par=est0,
#                  t=timepoints,
#                  fn=GetMinus2LL,
#                  lower = lboundary,
#                  upper = uboundary,
#                  control = list(tol=1e-10), # default 1e-6
#                  dataset=dataset_retrorsine,
#                  runtable=runtable_retrorsine)
# out_hjkb
# 
# # Save vector as .rda file
# save(out_hjkb, file="out_hjkb.rda")

# load file
load(file = "out_hjkb.rda")

# Assign optimised parameters to vector 
est_hat_hjkb <- out_hjkb$par

#--- AKAIKE INFORMATION CRITERION 
minus2LL <- out_hjkb$value
minus2LL
AIC <- 2*length(est0) + minus2LL
AIC

###############################################################################
# PART: MCMC
###############################################################################

#--- PERFORM MCMC (PACKAGE FME) TO ESTIMATE PARAMETER DISTRIBUTION
numberiter <- 200000 
burnin     <- 0.1*numberiter
update     <- 0.1*numberiter
dr         <- 6
dispersion <- sample(seq(0.7,1.3,by=0.1),length(est_hat_hjkb),replace=TRUE)
dispersion2 <- sample(seq(0.7,1.3,by=0.1),length(est_hat_hjkb),replace=TRUE)

# # Generate first chain
# out_MCMC1 <- modMCMC(f=GetMinus2LL,
#                      p=est_hat_hjkb,   # initial values for parameters
#                      jump=NULL,          # SD of proposal normal distribution; if NULL 10% of p
#                      # it can be efficient to use covar from model fit
#                      lower=lboundary,
#                      upper=uboundary,
#                      prior=NULL,         # -2log(parameter prior probability); NULL: non-informative prior, all parameters are equally likely
#                      var0=NULL,          # NULL: it is assumned that model variance is 1 and the return element from f is -2logL
#                      wvar0=NULL,         # "weight" for initial model variance, NULL: error variances are assumed to be fixed
#                      n0=NULL,            # parameter used for weihghing initial model variance, NULL: n0=wvar0*n
#                      niter=numberiter,   # no of iterations for the MCMC
#                      updatecov=update,   # setting updatecov smaller than niter will trigger adaptive MH,
#                      # proposal distribution is only updated during burnin when burninlenth is positive
#                      burninlength=burnin,# no of initial iterations to be removed, about 10% of niter
#                      ntrydr=dr,          # max no of tries for delayed rejection procedure
#                      dataset=dataset_retrorsine,
#                      runtable=runtable_retrorsine,
#                      t=timepoints)
# 
# save(out_MCMC1, file="out_MCMC1.rda")
# 
# # Generate second chain
# out_MCMC2 <- modMCMC(f=GetMinus2LL,
#                      p=dispersion*est_hat_hjkb,
#                      jump=NULL,
#                      lower=lboundary,
#                      upper=uboundary,
#                      prior=NULL,
#                      var0=NULL,
#                      wvar0=NULL,
#                      n0=NULL,
#                      niter=numberiter,
#                      updatecov=update,
#                      burninlength=burnin,
#                      ntrydr=dr,
#                      dataset=dataset_retrorsine,
#                      runtable=runtable_retrorsine,
#                      t=timepoints)
# 
# save(out_MCMC2, file="out_MCMC2.rda")
# 
# 
# # Generate third chain
# out_MCMC3 <- modMCMC(f=GetMinus2LL,
#                      p=dispersion2*est_hat_hjkb,
#                      jump=NULL,
#                      lower=lboundary,
#                      upper=uboundary,
#                      prior=NULL,
#                      var0=NULL,
#                      wvar0=NULL,
#                      n0=NULL,
#                      niter=numberiter,
#                      updatecov=update,
#                      burninlength=burnin,
#                      ntrydr=dr,
#                      dataset=dataset_retrorsine,
#                      runtable=runtable_retrorsine,
#                      t=timepoints)
# 
# save(out_MCMC3, file="out_MCMC3.rda")

# load files
load(file = "out_MCMC1.rda")
load(file = "out_MCMC2.rda")
load(file = "out_MCMC3.rda")

# Join  chains
out_MCMC <- rbind(out_MCMC1$pars,
                  out_MCMC2$pars,
                  out_MCMC3$pars)

# MCMC output can be used as functions from the coda package
MC1 <- as.mcmc(out_MCMC1$pars)
MC2 <- as.mcmc(out_MCMC2$pars)
MC3 <- as.mcmc(out_MCMC3$pars)

# Gelman-Rubin convergence dignostc
combinedchains <- mcmc.list(MC1,MC2,MC3) 

gelmandiag     <- gelman.diag(combinedchains, confidence = 0.95,
                              autoburnin = TRUE, multivariate = TRUE)
gelmandiag

# pairs plot
sample_params <- out_MCMC[sample(nrow(out_MCMC),size=1000,replace = TRUE),]

sample_params_plot <- subset(sample_params, select = c(ka,CLcanefr,
                                                       fDHRDNA,lambda1DHRDNA,
                                                       lambda2DHRDNA,kDHRDNA,
                                                       fDHRGSH,lambda1DHRGSH,
                                                       fDHRPROT,lambda1DHRPROT))

sample_params_plot   <- sample_params_plot[seq(1, NROW(sample_params_plot),1),]
pairs(sample_params_plot, diag.panel = panel.hist,upper.panel = panel.cor,
      panel=panel.smooth,cex.labels = 1.5, font.labels = 2, cex = 1.5,pch = 21)


# summary of posterior distribution

# Initialize
summary<- data.frame(matrix(NA,nrow = 1, ncol = 1))

# loop through posterior 
for(i in colnames(out_MCMC)) {
  print(i)
  
  #summary <- summary(out_MCMC[,i])
  mode_i    <- map_estimate(out_MCMC[,i])
  # mode_i    <- estimate_mode(out_MCMC[,i])
  mean_i    <- mean(out_MCMC[,i])
  median_i  <- median(out_MCMC[,i])
  ci_i      <- ci(out_MCMC[,i], ci=0.95,method = "HDI")
  summary_i <- rbind(median=median_i,mean=mean_i,mode=mode_i,ci=ci_i$CI,
                     ci_low=ci_i$CI_low,ci_high=ci_i$CI_high)
  colnames(summary_i) <- i
  
  print(summary_i)
  summary <- cbind(summary, summary_i)
}

# remove initialisation column
summary <- summary[colSums(!is.na(summary)) > 0]

View(summary)

# minus 2LL out_MCMC1
minus2LL_out_MCMC1 <- out_MCMC1$bestfunp
minus2LL_out_MCMC1
AIC <- 2*length(est0) + minus2LL_out_MCMC1
AIC


###############################################################################
# PREDICTION USING ESTIMATED PARAMETERS FROM MCMC ANALYSIS
###############################################################################

#--- SOLVE ODE SYSTEM USING OPTIMIZED PARAMETERS 

#mode
est_mode <- unlist(summary[3,])

dataset_retrorsine_pred <- PredRxODE(est=est_mode,
                                     t=c(timepoints,seq(0,672,by=0.1)),
                                     dataset=dataset_retrorsine, 
                                     runtable=runtable_retrorsine, 
                                     add_timepoints = TRUE,
                                     add_error = FALSE)

View(dataset_retrorsine_pred)


#--- SOLVE  ODE SYSTEM TO PREDICT PO DOSE

#--- ADD RUNS TO RUNTABLE
id_max <- max(runtable_retrorsine$RUNID)

add_run  <- data.frame(1/10^3/351.399*10^9, "nmolperkg", "mouse9weeks","po",
                       id_max+1)
add_run_next  <- data.frame(1/10^3/351.399*10^9, "nmolperkg", "rat10weeks","po",
                            id_max+2)
names(add_run)  <- c("DOSE","DOSEUNIT","SPECIESTYPE","ADM","RUNID")
names(add_run_next)  <- c("DOSE","DOSEUNIT","SPECIESTYPE","ADM","RUNID")


runtable_retrorsine <- rbind(runtable_retrorsine,add_run,add_run_next)

dataset_retrorsine_pred2 <- PredRxODE2(est=est_mode,
                                     t=seq(0,24,by=0.01),
                                     dataset=dataset_retrorsine,
                                     runtable=runtable_retrorsine,
                                     add_timepoints = FALSE,
                                     add_error = FALSE)

#################################################################################
# PART: SSR                                   
#################################################################################
  
dataset_retrorsine_pred$RESIDUALS        <- NA
dataset_retrorsine_pred$SQUAREDRESIDUALS <- NA
dataset_retrorsine_pred$RESIDUALS        <- (dataset_retrorsine_pred$Y - 
                                               dataset_retrorsine_pred$YPRED)
dataset_retrorsine_pred$SQUAREDRESIDUALS <- (dataset_retrorsine_pred$Y - 
                                               dataset_retrorsine_pred$YPRED)^2

# SSR for compound-species RTS
SSR_RTS <- sum(dataset_retrorsine_pred$SQUAREDRESIDUALS[which(
  (!is.na(dataset_retrorsine_pred$SQUAREDRESIDUALS)) & 
    ((dataset_retrorsine_pred$ID=="Yang_2017_retrorsine_mouse_plasma") |
       (dataset_retrorsine_pred$ID=="Li_2022_retrorsine_mouse_serum") |
       (dataset_retrorsine_pred$ID=="White_1977_retrorsine_rat_bile")    |
       (dataset_retrorsine_pred$ID=="Chu_1991_retrorsine_mouse_urine") |
       (dataset_retrorsine_pred$ID=="Chu_1991_retrorsine_rat_urine"))
)]
)
  

# SSR for compound-species GSHconjugates
SSR_GSHconjugates <- sum(dataset_retrorsine_pred$SQUAREDRESIDUALS[which(
  (!is.na(dataset_retrorsine_pred$SQUAREDRESIDUALS)) &
    (dataset_retrorsine_pred$ID=="Yang_2017_GSHconjugates_mouse_liver")
)]
)
  
# SSR for compound-species proteinadducts
SSR_proteinadducts <- sum(dataset_retrorsine_pred$SQUAREDRESIDUALS[which(
  (!is.na(dataset_retrorsine_pred$SQUAREDRESIDUALS)) &
    (dataset_retrorsine_pred$ID=="Yang_2017_proteinadducts_mouse_liver")
)]
)

# SSR for compound-species DNAadducts
SSR_DNAadducts <- sum(dataset_retrorsine_pred$SQUAREDRESIDUALS[which(
  (!is.na(dataset_retrorsine_pred$SQUAREDRESIDUALS)) &
    (dataset_retrorsine_pred$ID=="Zhu_2017_DNAadducts_mouse_liver")
)]
)

SSR <- SSR_RTS+SSR_GSHconjugates+SSR_proteinadducts+SSR_DNAadducts


###############################################################################
# PART: POSTPROCESSING
###############################################################################

#--- DETRANSFORM PREDICTIONS 

for(row in 1:nrow(runtable_retrorsine)) {
  results_long[[row]]$value <- DeTransData(results_long[[row]]$value)
}  

dataset_retrorsine_pred$Y     <- DeTransData(dataset_retrorsine_pred$Y)
dataset_retrorsine_pred$YPRED <- DeTransData(dataset_retrorsine_pred$YPRED)

#--- CONVERT Y AND YPRED INTO Y_PERCENT AND YPRED_PERCENT IN DATASET 
#--- CONVERT LOQ AND LOD INTO LOQ_PERCENT AND LOD_PERCENT IN DATASET 

# Initialize new columns
dataset_retrorsine_pred$Y_PERCENT <- NA
dataset_retrorsine_pred$Y_PERCENTUNIT <- "%"
dataset_retrorsine_pred$YPRED_PERCENT <- NA
dataset_retrorsine_pred$YPRED_PERCENTUNIT <- "%"

dataset_retrorsine_pred$LOQ_PERCENT <- NA
dataset_retrorsine_pred$LOQ_PERCENTUNIT <- "%"
dataset_retrorsine_pred$LOD_PERCENT <- NA
dataset_retrorsine_pred$LOD_PERCENTUNIT <- "%"

#--- CONVERT value INTO value_PERCENT AND CONVERT LONG FORMAT TO WIDE FORMAT
# FOR LOOP THROUGH ROWS OF RUNTABLE
for(row in 1:nrow(runtable_retrorsine)) {
  
  # Identify dose from runtable (nmol/kg)
  dose <- runtable_retrorsine$DOSE[which(runtable_retrorsine$RUNID==row)]
  
  # Identify speciestype
  speciestype <- runtable_retrorsine$SPECIESTYPE[which(runtable_retrorsine$RUNID==row)]
  
  # Identify bodyweight (kg) of speciestype
  bodyweight <- get(paste("theta",speciestype,sep="_"))["bw"]
  
  
  results_long[[row]]$value_PERCENT <- results_long[[row]]$value*100/(dose*bodyweight)
}

# Extract ids
ids <- unique(dataset_retrorsine$ID)


#--- DATASET

# Loop through ids
for (ID in ids) {

  # Identify runid for id
  runid <-  unique(dataset_retrorsine$RUNID[which(dataset_retrorsine$ID==ID)])
  
  # Identify dose from runtable (nmol/kg)
  dose <- runtable_retrorsine$DOSE[which(runtable_retrorsine$RUNID==runid)]
  
  # Identify speciestype
  speciestype <- runtable_retrorsine$SPECIESTYPE[which(runtable_retrorsine$RUNID==runid)]
  
  # Identify bodyweight (kg) of speciestype
  bodyweight <- get(paste("theta",speciestype,sep="_"))["bw"]
  
  # Identify rownumbers of ID
  is_rownumber  <- which(dataset_retrorsine_pred$ID == ID & is.na(dataset_retrorsine_pred$DOSE))
  is_rownumber2 <- which(dataset_ppd$ID == ID & is.na(dataset_ppd$DOSE))
  
   # Convert Y
  dataset_retrorsine_pred$Y_PERCENT[is_rownumber] <- dataset_retrorsine_pred$Y[is_rownumber]*100/(dose*bodyweight)  

  # Convert YPRED
  dataset_retrorsine_pred$YPRED_PERCENT[is_rownumber] <- dataset_retrorsine_pred$YPRED[is_rownumber]*100/(dose*bodyweight)

  #Convert LOQ
  dataset_retrorsine_pred$LOQ_PERCENT[is_rownumber] <- dataset_retrorsine_pred$LOQ[is_rownumber]*100/(dose*bodyweight)
  
  #Convert LOD
  dataset_retrorsine_pred$LOD_PERCENT[is_rownumber] <- dataset_retrorsine_pred$LOD[is_rownumber]*100/(dose*bodyweight)
  
  }


###############################################################################
# PART: Plots                                                            
###############################################################################

# Set working directory to subfolder
setwd(paste0(PATH,"/Figures"))

# Subset of datasets
dataset_retrorsine_plot     <- subset(dataset_retrorsine_pred,
                                      is.na(dataset_retrorsine_pred$DOSE))

#contains also data for which no simulation is done

# Dataset for plotting in long format (TYPE YPRED or Y)
dataset_retrorsine_long_to_be   <- subset(dataset_retrorsine_plot, 
                                          select = c(ID,COMPARTMENT,SPECIES,X,Y,YPRED))
dataset_retrorsine_long_to_be_2 <- gather(dataset_retrorsine_long_to_be, TYPE, Y, 5:6, factor_key=TRUE)                                  
dataset_retrorsine_long         <- subset(dataset_retrorsine_long_to_be_2, (!(is.na(Y))))

# Dataset for plotting in long format (TYPE YPRED_PERCENT or Y_PERCENT)
dataset_retrorsine_long_to_be_percent   <- subset(dataset_retrorsine_plot, 
                                                  select = c(ID,COMPARTMENT,SPECIES,X,Y_PERCENT,YPRED_PERCENT))
dataset_retrorsine_long_to_be_percent_2 <- gather(dataset_retrorsine_long_to_be_percent, TYPE, Y, 5:6, factor_key=TRUE)                                  
dataset_retrorsine_long_percent         <- subset(dataset_retrorsine_long_to_be_percent_2, (!(is.na(Y))))



###############################################################################
# DIAGNOSTIC PLOTS
###############################################################################

### PLOT PREDICTED VS RESIDUAL ################################################

# RTS plasma
dataset <- dataset_retrorsine_plot[which((dataset_retrorsine_plot$ID=="Yang_2017_retrorsine_mouse_plasma")&
                                           (!(is.na(dataset_retrorsine_plot$Y)))),]
aest    <- aes(x=YPRED, y=RESIDUALS)
xlabel  <- "Predicted"
ylabel  <- "Residual"
title   <- "RET - Mouse plasma - Yang et al. (2017)"
residuals_RTS_plasma <- GetPlotResiduals(dataset,aest,xlabel,ylabel,title)
ggsave("residuals_RTS_plasma.png", bg="white",width = 6,height=3)


# RTS serum
dataset <- dataset_retrorsine_plot[which((dataset_retrorsine_plot$ID=="Li_2022_retrorsine_mouse_serum")&
                                           (!(is.na(dataset_retrorsine_plot$Y)))),]
aest    <- aes(x=YPRED, y=RESIDUALS)
xlabel  <- "Predicted"
ylabel  <- "Residual"
title   <- "RET - Mouse plasma - Li et al. (2022)"
residuals_RTS_serum <- GetPlotResiduals(dataset,aest,xlabel,ylabel,title)
ggsave("residuals_RTS_serum.png", bg="white",width = 6,height=3)


# RTS bile
dataset <- dataset_retrorsine_plot[which((dataset_retrorsine_plot$ID=="White_1977_retrorsine_rat_bile")&
                                           (!(is.na(dataset_retrorsine_plot$Y)))),]
aest    <- aes(x=YPRED, y=RESIDUALS)
xlabel  <- "Predicted"
ylabel  <- "Residual"
title   <- "RET - Rat bile - White (1977) "
residuals_RTS_bile <- GetPlotResiduals(dataset,aest,xlabel,ylabel,title)
ggsave("residuals_RTS_bile.png", bg="white",width = 6,height=3)


# RTS urine mouse
dataset <- dataset_retrorsine_plot[which((dataset_retrorsine_plot$ID=="Chu_1991_retrorsine_mouse_urine")&
                                           (!(is.na(dataset_retrorsine_plot$Y)))),]
aest    <- aes(x=YPRED, y=RESIDUALS)
xlabel  <- "Predicted"
ylabel  <- "Residual"
title   <- "RET - Mouse urine - Chu and Segall (1991)"
residuals_RTS_urine_m <- GetPlotResiduals(dataset,aest,xlabel,ylabel,title)
ggsave("residuals_RTS_urine_m.png", bg="white",width = 6,height=3)


# RTS urine rat
dataset <- dataset_retrorsine_plot[which((dataset_retrorsine_plot$ID=="Chu_1991_retrorsine_rat_urine")&
                                           (!(is.na(dataset_retrorsine_plot$Y)))),]
aest    <- aes(x=YPRED, y=RESIDUALS)
xlabel  <- "Predicted"
ylabel  <- "Residual"
title   <- "RET - Rat urine - Chu and Segall (1991)"
residuals_RTS_urine_r <- GetPlotResiduals(dataset,aest,xlabel,ylabel,title)
ggsave("residuals_RTS_urine_r.png", bg="white",width = 6,height=3)

# GSHconjugates

Wliv     <- 1.3725*10^3 #mg            (Brown et al. 1997)
Cprotliv <- 23.84       #mg/100 mg liv (Verma and Chakraborty 2008)
WIS      <- 1.25*10^-8  #g             (personal correspondance)
MIS      <- 391.5       #g/mol         (NCBI 2021)  
Rprot    <- 0.787       #ratio of protein amount of supernatant 
#to liver homogenate (personal correspondance)


dataset <- dataset_retrorsine_plot[which((dataset_retrorsine_plot$ID=="Yang_2017_GSHconjugates_mouse_liver")&
                                           (!(is.na(dataset_retrorsine_plot$Y)))),]
aest    <- aes(x=YPRED/Cprotliv/Wliv*100/WIS*MIS/Rprot/10^9, y=RESIDUALS)
xlabel  <- "Predicted"
ylabel  <- "Residual"
title   <- "DHR:GSH - Mouse liver - Yang et al. (2017)"
residuals_GSHconjugates_liver <- GetPlotResiduals(dataset,aest,xlabel,ylabel,title)
ggsave("residuals_GSHconjugates_liver.png", bg="white",width = 6,height=3)


# Protein adducts

Wliv     <- 1.3725*10^3 #mg            (Brown et al. 1997)
Cprotliv <- 23.84       #mg/100 mg liv (Verma and Chakraborty 2008)
WIS      <- 1.25*10^-8  #g             (personal correspondance)
MIS      <- 284.74      #g/mol         (NCBI 2021)  
Rprot    <- 0.787       #ratio of protein amount of supernatant 
#to liver homogenate (personal correspondance)


dataset <- dataset_retrorsine_plot[which((dataset_retrorsine_plot$ID=="Yang_2017_proteinadducts_mouse_liver")&
                                           (!(is.na(dataset_retrorsine_plot$Y)))),]
aest    <- aes(x=YPRED/Cprotliv/Wliv*100/WIS*MIS/Rprot/10^9, y=RESIDUALS)
xlabel  <- "Predicted"
ylabel  <- "Residual"
title   <- "DHR:PROT - Mouse liver - Yang et al. (2017)"
residuals_proteinadducts_liver <- GetPlotResiduals(dataset,aest,xlabel,ylabel,title)
ggsave("residuals_proteinadducts_liver.png", bg="white",width = 6,height=3)


# DNAadducts

Na       <- 6.02214076*10^23 #1/mol                         (IUPAC 1997)
SFgenome <- 2728222451       #base pairs/mouse haploid cell (NCBI 2021) 
SFliv    <- 128*10^6         #cells/g liver                 (Ring et al. 2011) 
Wliv     <- 1.3725           #g                             (Brown et al. 1997)


dataset <- dataset_retrorsine_plot[which((dataset_retrorsine_plot$ID=="Zhu_2017_DNAadducts_mouse_liver")&
                                           (!(is.na(dataset_retrorsine_plot$Y)))),]
aest    <- aes(x=YPRED/SFgenome/2/2*10^8/SFliv/Wliv*Na/10^9, y=RESIDUALS)
xlabel  <- "Predicted"
ylabel  <- "Residual"
title   <- "DHR:DNA - Mouse liver - Yang et al. (2017)"
residuals_DNAadducts_liver <- GetPlotResiduals(dataset,aest,xlabel,ylabel,title)
ggsave("residuals_DNAadducts_liver.png", bg="white",width = 6,height=3)

diagnosticplot <- plot_grid(residuals_RTS_plasma,residuals_RTS_serum,
                            residuals_RTS_bile,residuals_RTS_urine_m,
                            residuals_RTS_urine_r,
                            residuals_GSHconjugates_liver,
                            residuals_proteinadducts_liver, 
                            residuals_DNAadducts_liver,ncol=2)

diagnosticplot



#--- PLOT FRACTION OF DOSE

#-- COLOR
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2",
  #"#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", #"khaki2",
  "maroon", "orchid1",
  #"deeppink1",
  "blue1", "steelblue4",
  "darkturquoise", #"darkorange4",
  "yellow4", "yellow3",
  "brown", "green1"
)



for (i in 1:length(runtable_retrorsine$RUNID)) {

  results_long[[i]] <- subset(results_long[[i]],
                              compartment!="fsum"
                              & compartment!="bag"
                              & compartment!="Aven"
                              & compartment!="Aart"
                              & compartment!="Aper"
                              & compartment!="GSH"
                              & compartment!="Aliv"
                              & compartment!="DHRGSH"
                              & compartment!="DHRPROT"
                              & compartment!="DHRDNA")


  results_long[[i]]$compartment <- factor(results_long[[i]]$compartment ,
                                          levels = c("Aadi",
                                                     "Abil",
                                                     "Abon",
                                                     "Abra",
                                                     "Agut",
                                                     "Ahea",
                                                     "Akid",
                                                     "Alivc",
                                                     "Alivvi",
                                                     "Alum",
                                                     "Alun",
                                                     "Amus",
                                                     "Apla",
                                                     "Aski",
                                                     "Aspl",
                                                     "Auri",
                                                     "AMetgut",
                                                     "AMetliv"))

  dataset          <- results_long[[i]]
  title            <- paste("RUNID",i,sep = "")
  xlimits          <- c(0,14)
  xbreaks          <- seq(0,14,by = 2)
  ylimits_log      <- c(1e-2,100)
  pdfname_log      <- paste("RUNID",i,".pdf", sep = "")
  pngname_log      <- paste("RUNID",i,".png", sep = "")


  results_plot <-
    ggplot(data=dataset,mapping=aes(color=factor(compartment)))+
    geom_line(aes(x=time, y=value_PERCENT), size=1.5)+
    labs(x="Time (h)",y="Fraction of dose (%)")+
    #ggtitle(title)+
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=18, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 26),
      axis.title        = element_text(size = 28),
      # axis.title.y      = element_text(margin=margin(r=15)),
      # axis.title.x      = element_text(margin=margin(t=12)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=18),
      legend.title      = element_text(size =18),
      legend.position   = "none",
      legend.background = element_rect(fill = NA),
      aspect.ratio=1,
      legend.key.width = unit(2, "cm"),
      legend.key=element_blank(),
      legend.direction = "vertical"

    )+
    guides(col = guide_legend(nrow = 6,title.position = "top")) +
    scale_y_log10(
      # breaks = scales::trans_breaks("log10", function(x) 10^x),
      # labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = ylimits_log,
      labels = function(x) sprintf("%g", x)
    )+
    scale_x_continuous(
      limits = xlimits,
      breaks = xbreaks)+
    annotation_logticks(sides="l")+scale_color_manual(values = c25)+
    guides(color=guide_legend(ncol=3,override.aes = list(size = 3))) 

  assign(paste("results_plot_runid",i,sep=""), results_plot)

}

results_plot_runid11 #po mouse
ggsave("mouse_po.png", bg="white", width = 6, height =6)
results_plot_runid12 # po rat
ggsave("rat_po.png", bg="white", width = 6, height =6)




###############################################################################
# MODEL TRAINING PLOTS
###############################################################################

#--- Retrorsine in plasma, bile, urine

dataset_retrorsine_long_RTS <-  subset(dataset_retrorsine_long, 
                                       ID=="Yang_2017_retrorsine_mouse_plasma"|
                                         ID=="White_1977_retrorsine_rat_bile"|
                                         ID=="Chu_1991_retrorsine_mouse_urine"|
                                         ID=="Chu_1991_retrorsine_rat_urine"|
                                         ID=="Li_2022_retrorsine_mouse_serum")


dataset1    <- subset(dataset_retrorsine_long_RTS, 
                      dataset_retrorsine_long_RTS$TYPE=="Y")
dataset2    <- subset(dataset_retrorsine_long_RTS, 
                      dataset_retrorsine_long_RTS$TYPE=="YPRED")
aest        <- aes(x=X, y=Y)
xlabel      <- "Time (h)"
ylabel      <- "Amount (nmol)"
title       <- ""
ylimits_log <- c(10^-2,10^5)
ylimits     <- c(0,100)
ybreaks     <- seq(0,100,by = 10)
xlimits     <- c(0,36)
xbreaks     <- seq(0,36,by = 4)

retrorsine <- 
  ggplot(data=dataset1,mapping=aes(color=ID))+
  geom_point(aes(x=X, y=Y, shape=factor(COMPARTMENT),fill=factor(ID)),size=6, 
             stroke=1.5, alpha=0.4)+
  geom_line(data=dataset2,aes(x=X, y=Y,linetype=factor(ID)), size=2)+
  labs(x=xlabel,y=ylabel)+
  ggtitle(title)+ 
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=26, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=22),
    legend.title      = element_text(size =22),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=1
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = ylimits_log
  )+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  annotation_logticks(sides="l")+
  scale_color_manual(values=c("deepskyblue3","darkorange3","black","darkorange",
                              "deepskyblue4"))+ 
  scale_fill_manual(values=c("deepskyblue3","darkorange3","black","darkorange",
                             "deepskyblue4","darkorange", "deepskyblue4"))+ 
  scale_shape_manual(values=c(22,21,1,24))+
  scale_linetype_manual(values=c("solid","longdash","solid","longdash","solid"))+
  geom_hline(yintercept=0.001749734, linetype="dashed", color = "black", 
             size=1.5)+
  geom_hline(yintercept=3.499469e-05, linetype="dashed", color = "black", 
             size=1.5)+
  annotate("text", x = 36, y = 0.006, label = "LLOQ",size=7, hjust=1, vjust=1)+
  annotate("text", x = 36, y = 0.00012, label = "LLOD",size=7, hjust=1, vjust=1)+
  annotate("text", x = Inf, y = Inf, label = "RET (plasma, urine, bile) ",
           size=10, hjust=1, vjust=1.3)

retrorsine
ggsave("retrorsine.png", bg="white", width = 7, height =7)



#--- GSH conjugates

Wliv     <- 1.3725*10^3 #mg            (Brown et al. 1997)
Cprotliv <- 23.84       #mg/100 mg liv (Verma and Chakraborty 2008)
WIS      <- 1.25*10^-8  #g             (personal correspondance)
MIS      <- 391.5       #g/mol         (NCBI 2021)  
Rprot    <- 0.787       #ratio of protein amount of supernatant 
#to liver homogenate (personal correspondance)


dataset_retrorsine_long_metabolites <-  subset(dataset_retrorsine_long, 
                                               ID=="Yang_2017_GSHconjugates_mouse_liver")

dataset1    <- subset(dataset_retrorsine_long_metabolites, 
                      dataset_retrorsine_long_metabolites$TYPE=="Y")
dataset2    <- subset(dataset_retrorsine_long_metabolites, 
                      dataset_retrorsine_long_metabolites$TYPE=="YPRED")
aest        <- aes(x=X, y=Y)
xlabel      <- "Time (h)"
ylabel      <- bquote("Area ratio "(A[RET]/A[IS])/"mg protein") 
#"Aanalyte/AIS/mg protein" 
title       <- ""
ylimits_log <- c(0.001,10)
xlimits     <- c(0,8)
xbreaks     <- seq(0,8,by = 2)

DHRGSH_orig <-   ggplot(data=dataset1,mapping=aes(color=ID))+
  geom_point(aes(x=X, y=Y/Cprotliv/Wliv*100/WIS*MIS/Rprot/10^9,fill=factor(ID)),
             shape=21,size=6, stroke=1.5,alpha=0.4)+
  geom_line(data=dataset2,aes(x=X, y=Y/Cprotliv/Wliv*100/WIS*MIS/Rprot/10^9), 
            size=2)+
  labs(x=xlabel,y=ylabel)+
  ggtitle(title)+ 
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=26, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 28),
    axis.title        = element_text(size = 26),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=22),
    legend.title      = element_text(size =22),
    legend.position   = "none", #bottom",
    legend.background = element_rect(fill = NA),
    aspect.ratio=1
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    #breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = ylimits_log
  )+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  annotation_logticks(sides="l")+
  scale_color_manual(values=c("#D55E00", "#CC79A7", "#999999"))+ 
  scale_fill_manual(values=c("#D55E00", "#CC79A7", "#999999"))+ 
  annotate("text", x = Inf, y = Inf, label = "DHR:GSH (liver) ",size=10,hjust=1, 
           vjust=1.3)

DHRGSH_orig
ggsave("DHRGSH_orig.png", bg="white", width = 7, height =7)


#--- protein adducts

Wliv     <- 1.3725*10^3 #mg            (Brown et al. 1997)
Cprotliv <- 23.84       #mg/100 mg liv (Verma and Chakraborty 2008)
WIS      <- 1.25*10^-8  #g             (personal correspondance)
MIS      <- 284.74      #g/mol         (NCBI 2021)  
Rprot    <- 0.787       #ratio of protein amount of supernatant 
#to liver homogenate (personal correspondance)

dataset_retrorsine_long_metabolites <-  subset(dataset_retrorsine_long, 
                                               ID=="Yang_2017_proteinadducts_mouse_liver")

dataset1    <- subset(dataset_retrorsine_long_metabolites, 
                      dataset_retrorsine_long_metabolites$TYPE=="Y")
dataset2    <- subset(dataset_retrorsine_long_metabolites, 
                      dataset_retrorsine_long_metabolites$TYPE=="YPRED")
aest        <- aes(x=X, y=Y)
xlabel      <- "Time (h)"
ylabel      <- bquote("Area ratio "(A[RET]/A[IS])/"mg protein") 
#"Aanalyte/AIS/mg protein" 
title       <- ""
ylimits_log <- c(0.01,100)
xlimits     <- c(0,36)
xbreaks     <- seq(0,36,by = 4)

DHRPROT_orig <-   ggplot(data=dataset1,mapping=aes(color=ID))+
  geom_point(aes(x=X, y=Y/Cprotliv/Wliv*100/WIS*MIS/Rprot/10^9,fill=factor(ID)),
             shape=21,size=6, stroke=1.5,alpha=0.4)+
  geom_line(data=dataset2,aes(x=X, y=Y/Cprotliv/Wliv*100/WIS*MIS/Rprot/10^9), 
            size=2)+
  labs(x=xlabel,y=ylabel)+
  ggtitle(title)+ 
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=26, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 28),
    axis.title        = element_text(size = 26),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=22),
    legend.title      = element_text(size =22),
    legend.position   = "none", #bottom",
    legend.background = element_rect(fill = NA),
    aspect.ratio=1
    
  )+
  # guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = ylimits_log
  )+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  annotation_logticks(sides="l")+
  scale_color_manual(values=c("#CC79A7", "#999999"))+ 
  scale_fill_manual(values=c("#CC79A7", "#999999"))+ 
  annotate("text", x = Inf, y = Inf, label = "DHR:PROT (liver) ",size=10,hjust=1, 
           vjust=1.3)

DHRPROT_orig
ggsave("DHRPROT_orig.png", bg="white", width = 7, height =7)



#--- DNA adducts

Na       <- 6.02214076*10^23 #1/mol                         (IUPAC 1997)
SFgenome <- 2728222451       #base pairs/mouse haploid cell (NCBI 2021) 
SFliv    <- 128*10^6         #cells/g liver                 (Ring et al. 2011) 
Wliv     <- 1.3725           #g                             (Brown et al. 1997)


dataset_retrorsine_long_metabolites <-  subset(dataset_retrorsine_long, 
                                               ID=="Zhu_2017_DNAadducts_mouse_liver")

dataset1    <- subset(dataset_retrorsine_long_metabolites, 
                      dataset_retrorsine_long_metabolites$TYPE=="Y")
dataset2    <- subset(dataset_retrorsine_long_metabolites, 
                      dataset_retrorsine_long_metabolites$TYPE=="YPRED")
xlabel      <- "Time (days)"
ylabel      <- bquote(paste("DNA adducts/",10^8," nucleotides"))   
title       <- ""
ylimits_log <- c(10^0,10^2)
# xlimits     <- c(0,672)
# xbreaks     <- seq(0,672,by = 200)
xlimits     <- c(0,28)
xbreaks     <- seq(0,28,by = 7)

DNAadducts_orig <-   ggplot(data=dataset1,mapping=aes(color=ID))+
  geom_point(aes(x=X/24, y=Y/SFgenome/2/2*10^8/SFliv/Wliv*Na/10^9,
                 fill=factor(ID)),shape=21,size=6, stroke=1.5,alpha=0.4)+
  geom_line(data=dataset2,aes(x=X/24, y=Y/SFgenome/2/2*10^8/SFliv/Wliv*Na/10^9), 
            size=2)+
  labs(x=xlabel,y=ylabel)+
  ggtitle(title)+ 
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=26, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 28),
    axis.title        = element_text(size = 26),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=22),
    legend.title      = element_text(size =22),
    legend.position   = "none", #bottom",
    legend.background = element_rect(fill = NA),
    aspect.ratio=1
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = ylimits_log
  )+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  annotation_logticks(sides="l")+
  scale_color_manual(values=c("#999999"))+ 
  scale_fill_manual(values=c("#999999"))+ 
  annotate("text", x = Inf, y = Inf, label = "DHR:DNA (liver) ",size=10, hjust=1, vjust=1.3)+
  annotation_logticks(sides="l",short = unit(2.5,"mm"),mid = unit(5,"mm"))

DNAadducts_orig
ggsave("DNAadducts_orig.png", bg="white", width = 7, height =7)



###############################################################################
# MODEL EVALUATION PLOTS
###############################################################################

#--- RTS in plasma and liver --------------------------------------------------

dataset_retrorsine_long_RTS <-  subset(dataset_retrorsine_long, 
                                       ID=="Yang_2018_retrorsine_mouse_liver"|
                                       ID=="Yang_2018_retrorsine_mouse_plasma"|
                                       ID=="Li_2022b_retrorsine_mouse_serum")


dataset1    <- subset(dataset_retrorsine_long_RTS, 
                      dataset_retrorsine_long_RTS$TYPE=="Y")
dataset2    <- subset(dataset_retrorsine_long_RTS, 
                      dataset_retrorsine_long_RTS$TYPE=="YPRED")
aest        <- aes(x=X, y=Y)
xlabel      <- "Time (h)"
ylabel      <- "Amount (nmol)"
title       <- ""
ylimits_log <- c(10^-2,10^5)
ylimits     <- c(0,100)
ybreaks     <- seq(0,100,by = 10)
xlimits     <- c(0,8)
xbreaks     <- seq(0,8,by = 2)


retrorsine_eval <- 
  ggplot(data=dataset1,mapping=aes(color=ID))+
  geom_point(aes(x=X, y=Y, shape=factor(COMPARTMENT),fill=factor(ID)),size=6, 
             stroke=1.5, alpha=0.4)+
  geom_line(data=dataset2,aes(x=X, y=Y), size=2)+
  labs(x=xlabel,y=ylabel)+
  ggtitle(title)+ 
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=26, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=22),
    legend.title      = element_text(size =22),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=1
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = ylimits_log
  )+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  annotation_logticks(sides="l")+
  scale_color_manual(values=c("#9370db","#009E73","deepskyblue4"))+ 
  scale_fill_manual(values=c("#9370db","#009E73","deepskyblue4"))+ 
  scale_shape_manual(values=c(23,21,1))+
  geom_hline(yintercept=0.001749734, linetype="dashed", color = "black", 
             size=1.5)+
  geom_hline(yintercept=3.499469e-05, linetype="dashed", color = "black", 
             size=1.5)+
  annotate("text", x = 8, y = 0.006, label = "LLOQ",size=7, hjust=1, vjust=1)+
  annotate("text", x = 8, y = 0.00012, label = "LLOD",size=7, hjust=1, vjust=1)+
  annotate("text", x = Inf, y = Inf, label = "RET (plasma, liver) ",size=10, 
           hjust=1, vjust=1.3)

retrorsine_eval
ggsave("retrorsine_eval.png", bg="white", width = 7, height =7)



#--- GSH conjugates 
Wliv     <- 1.3725*10^3 #mg            (Brown et al. 1997)
Cprotliv <- 23.84       #mg/100 mg liv (Verma and Chakraborty 2008)
WIS      <- 1.6*10^-8  #g             
MIS      <- 391.5       #g/mol         (NCBI 2021)  
Rprot    <- 0.787       #ratio of protein amount of supernatant 
#to liver homogenate (personal correspondance)


#Scaling factor of predictions
SF <- 12

dataset_retrorsine_long_metabolites <-  subset(dataset_retrorsine_long, 
                                               ID=="Yang_2018_GSHconjugates_mouse_liver")

dataset1    <- subset(dataset_retrorsine_long_metabolites, 
                      dataset_retrorsine_long_metabolites$TYPE=="Y")
dataset2    <- subset(dataset_retrorsine_long_metabolites, 
                      dataset_retrorsine_long_metabolites$TYPE=="YPRED")

aest        <- aes(x=X, y=Y)
xlabel      <- "Time (h)"
ylabel      <- bquote("Area ratio "(A[RET]/A[IS])/"mg protein") 
#"Aanalyte/AIS/mg protein" 
y2label     <- "PBTK model prediction"
title       <- ""
ylimits_log <- c(10^-3,10^1)
xlimits     <- c(0,8)
xbreaks     <- seq(0,8,by = 2)

DHRGSH_orig_eval <-   ggplot(data=dataset1,mapping=aes(color=ID))+
  geom_point(aes(x=X, y=Y/Cprotliv/Wliv*100/WIS*MIS/Rprot/10^9,fill=factor(ID)),
             shape=21,size=6, stroke=1.5,alpha=0.4)+
  geom_line(data=dataset2,aes(x=X, y=Y*SF/Cprotliv/Wliv*100/WIS*MIS/Rprot/10^9), 
            size=2)+
  labs(x=xlabel,y=ylabel)+
  ggtitle(title)+ 
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=23, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 24),
    axis.title        = element_text(size = 23),
    axis.title.y.right = element_text(size=23,angle = 90),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=23),
    legend.title      = element_text(size =23),
    legend.position   = "none", #bottom",
    legend.background = element_rect(fill = NA),
    aspect.ratio=1
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = ylimits_log,
    sec.axis = sec_axis( ~./SF, name=y2label,
                         #labels = function(x) sprintf("%g", x)))+
                         breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x))))+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  annotation_logticks(sides="l")+
  scale_color_manual(values=c("#D55E00","#CC79A7"))+ 
  scale_fill_manual(values=c("#D55E00","#CC79A7"))+ 
  annotate("text", x = Inf, y = Inf, label = "DHR:GSH (liver) ",size=8, hjust=1, 
           vjust=1.3)

DHRGSH_orig_eval
ggsave("DHRGSH_orig_eval.png", bg="white", width = 7, height =7)






#--- protein adducts
Wliv     <- 1.3725*10^3 #mg            (Brown et al. 1997)
Cprotliv <- 23.84       #mg/100 mg liv (Verma and Chakraborty 2008)
WIS      <- 1.25*10^-8  #g             (personal correspondance)
MIS      <- 284.74      #g/mol         (NCBI 2021)  
Rprot    <- 0.787       #ratio of protein amount of supernatant 
#to liver homogenate (personal correspondance)

#Scaling factor of predictions
SF <- 150

dataset_retrorsine_long_metabolites <-  subset(dataset_retrorsine_long, 
                                               ID=="Yang_2018_proteinadducts_mouse_liver")

dataset1    <- subset(dataset_retrorsine_long_metabolites, 
                      dataset_retrorsine_long_metabolites$TYPE=="Y")
dataset2    <- subset(dataset_retrorsine_long_metabolites, 
                      dataset_retrorsine_long_metabolites$TYPE=="YPRED")

xlabel      <- "Time (h)"
ylabel      <- bquote("Area ratio "(A[RET]/A[IS])/"mg protein") 
#"Aanalyte/AIS/mg protein" 
y2label     <- "PBTK model prediction"
title       <- ""
ylimits_log <- c(0.01,100)
xlimits     <- c(0,8)
xbreaks     <- seq(0,8,by = 2)

DHRPROT_orig_eval <-   ggplot(data=dataset1,mapping=aes(color=ID))+
  geom_point(aes(x=X, y=Y/Cprotliv/Wliv*100/WIS*MIS/Rprot/10^9,fill=factor(ID)),
             shape=21,size=6, stroke=1.5,alpha=0.4)+
  geom_line(data=dataset2,aes(x=X, y=Y*SF/Cprotliv/Wliv*100/WIS*MIS/Rprot/10^9), 
            size=2)+
  labs(x=xlabel,y=ylabel)+
  ggtitle(title)+ 
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=23, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 24),
    axis.title        = element_text(size = 23),
    axis.title.y.right = element_text(size=23,angle = 90),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=23),
    legend.title      = element_text(size =23),
    legend.position   = "none", #bottom",
    legend.background = element_rect(fill = NA),
    aspect.ratio=1
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = ylimits_log,
    sec.axis = sec_axis( ~./SF, name=y2label,
                         breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x))))+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  annotation_logticks(sides="l")+
  scale_color_manual(values=c("#CC79A7"))+
  scale_fill_manual(values=c("#CC79A7"))+ 
  annotate("text", x = Inf, y = Inf, label = "DHR:PROT (liver) ",size=8, hjust=1, vjust=1.3)

DHRPROT_orig_eval
ggsave("DHRPROT_orig_eval.png", bg="white", width = 7, height =7)

#--- DNA adducts

Na       <- 6.02214076*10^23 #1/mol                         (IUPAC 1997)
SFgenome <- 2728222451       #base pairs/mouse haploid cell (NCBI 2021) 
SFliv    <- 128*10^6         #cells/g liver                 (Ring et al. 2011) 
Wliv     <- 1.3725           #g                             (Brown et al. 1997)


dataset_retrorsine_long_DNAadducts <-  subset(dataset_retrorsine_long, 
                                              ID=="Zhu_2017_DNAadducts_mouse_liver_10"|
                                                ID=="Zhu_2017_DNAadducts_mouse_liver_20"|
                                                ID=="Zhu_2017_DNAadducts_mouse_liver_40"|
                                                ID=="Zhu_2017_DNAadducts_mouse_liver_60")

dataset1    <- subset(dataset_retrorsine_long_DNAadducts, 
                      dataset_retrorsine_long_DNAadducts$TYPE=="Y")
dataset2    <- subset(dataset_retrorsine_long_DNAadducts, 
                      dataset_retrorsine_long_DNAadducts$TYPE=="YPRED")
xlabel      <- "Time (days)"
ylabel      <- bquote(paste("DNA adducts/",10^8," nucleotides"))   
title       <- ""
ylimits_log <- c(1,100)
# xlimits     <- c(0,672)
# xbreaks     <- seq(0,672,by = 200)
xlimits     <- c(0,28)
xbreaks     <- seq(0,28,by = 7)

DNAadducts_orig_eval <-   ggplot(data=dataset1,mapping=aes(color=ID))+
  geom_point(aes(x=X/24, y=Y/SFgenome/2/2*10^8/SFliv/Wliv*Na/10^9,
                 shape=factor(ID),fill=factor(ID)),size=6, stroke=1.5,alpha=0.4)+
  geom_line(data=dataset2,aes(x=X/24, y=Y/SFgenome/2/2*10^8/SFliv/Wliv*Na/10^9),
            size=2)+
  labs(x=xlabel,y=ylabel)+
  ggtitle(title)+ 
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 2),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=26, hjust = 0.5),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 28),
    axis.title        = element_text(size = 26),
    # axis.title.y      = element_text(margin=margin(r=15)),
    # axis.title.x      = element_text(margin=margin(t=12)),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=22),
    legend.title      = element_text(size =22),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=1
    
  )+
  guides(col = guide_legend(nrow = 6,title.position = "top")) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = ylimits_log
  )+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks)+
  annotation_logticks(sides="l")+
  scale_color_manual(values=c("#e0e0e0","#c1c1c1","#999999","#5b5b5b"))+ 
  scale_fill_manual(values=c("#e0e0e0","#c1c1c1","#999999","#5b5b5b"))+ 
  scale_shape_manual(values=c(24,22,21,25))+
  annotate("text", x = Inf, y = Inf, label = "DHR:DNA (liver) ",size=10,hjust=1,
           vjust=1.3)+
  annotation_logticks(sides="l",short = unit(2.5,"mm"),mid = unit(5,"mm"))

DNAadducts_orig_eval
ggsave("DNAadducts_orig_eval.png", bg="white", width = 7, height =7)



###############################################################################
# Calculate Parameters
###############################################################################

#g
Vliv_r <- theta_rat10weeks["Vliv"]*1000
#mL/min/g liver
CLbile_r <- est_mode["CLcanefr"]*1000/60/Vliv_r

#L/h
CLliv_m <- theta_mouse9weeks["Qliv"]*theta_mouse9weeks["fuP"]*theta_mouse9weeks["CLactin"]*((theta_mouse9weeks["Vmaxliv"]/theta_mouse9weeks["Km"])+0)/
  (theta_mouse9weeks["Qliv"]*((theta_mouse9weeks["Vmaxliv"]/theta_mouse9weeks["Km"])+0) + theta_mouse9weeks["fuP"]*theta_mouse9weeks["CLactin"]*((theta_mouse9weeks["Vmaxliv"]/theta_mouse9weeks["Km"])+0 ))

CLliv_r <- theta_rat10weeks["Qliv"]*theta_rat10weeks["fuP"]*theta_rat10weeks["CLactin"]*((theta_rat10weeks["Vmaxliv"]/theta_rat10weeks["Km"])+est_mode["CLcanefr"])/
  (theta_rat10weeks["Qliv"]*((theta_rat10weeks["Vmaxliv"]/theta_rat10weeks["Km"])+est_mode["CLcanefr"]) + theta_rat10weeks["fuP"]*theta_rat10weeks["CLactin"]*((theta_rat10weeks["Vmaxliv"]/theta_rat10weeks["Km"])+est_mode["CLcanefr"] ))

CLliv_ws_m <- (theta_mouse9weeks["Qliv"]*theta_mouse9weeks["fuP"]*theta_mouse9weeks["Vmaxliv"]/theta_mouse9weeks["Km"])/
  (theta_mouse9weeks["Qliv"]+theta_mouse9weeks["fuP"]*theta_mouse9weeks["Vmaxliv"]/theta_mouse9weeks["Km"])

CLliv_ws_r <- (theta_rat10weeks["Qliv"]*theta_rat10weeks["fuP"]*theta_rat10weeks["Vmaxliv"]/theta_rat10weeks["Km"])/
  (theta_rat10weeks["Qliv"]+theta_rat10weeks["fuP"]*theta_rat10weeks["Vmaxliv"]/theta_rat10weeks["Km"])


Fh_m <- 1-(CLliv_m/theta_mouse9weeks["Qliv"])
Fh_r <- 1-(CLliv_r/theta_rat10weeks["Qliv"])

CLgut_m <- (theta_mouse9weeks["Qgut"]*theta_mouse9weeks["fuP"]*0.1*theta_mouse9weeks["Vmaxliv"]/theta_mouse9weeks["Km"])/
  (theta_mouse9weeks["Qgut"]+theta_mouse9weeks["fuP"]*0.1*theta_mouse9weeks["Vmaxliv"]/theta_mouse9weeks["Km"])
CLgut_r <- (theta_rat10weeks["Qgut"]*theta_rat10weeks["fuP"]*0.1*theta_rat10weeks["Vmaxliv"]/theta_rat10weeks["Km"])/
  (theta_rat10weeks["Qgut"]+theta_rat10weeks["fuP"]*0.1*theta_rat10weeks["Vmaxliv"]/theta_rat10weeks["Km"])

Fg_m <- 1-(CLgut_m/theta_mouse9weeks["Qgut"])
Fg_r <- 1-(CLgut_r/theta_rat10weeks["Qgut"])

F_m <- theta_mouse9weeks["Fa"]*Fh_m*Fg_m
F_r <- theta_rat10weeks["Fa"]*Fh_r*Fg_r

#L/h
CLuri_m <- theta_mouse9weeks["fuP"]*theta_mouse9weeks["GFR"]
CLuri_r <- theta_rat10weeks["fuP"]*theta_rat10weeks["GFR"]
# #mL/min/kg bw
# CLuri_m <- theta_mouse9weeks["fuP"]*theta_mouse9weeks["GFR"]*1000/60/theta_mouse9weeks["bw"]
# CLuri_r <- theta_rat10weeks["fuP"]*theta_rat10weeks["GFR"]*1000/60/theta_rat10weeks["bw"]


#L/h
CLtot_m <- CLliv_m + CLgut_m + theta_mouse9weeks["fuP"]*theta_mouse9weeks["GFR"]
CLtot_r <- CLliv_r + CLgut_r + theta_rat10weeks["fuP"]*theta_rat10weeks["GFR"]

#Fraction of CLuri that is total retrorsine CL
fCLuri_m <- 100*CLuri_m/CLtot_m
fCLuri_r <- 100*CLuri_r/CLtot_r

#Fraction of CLgut that is total retrorsine CL
fCLgut_m <- 100*CLgut_m/CLtot_m
fCLgut_r <- 100*CLgut_r/CLtot_r

#Fraction of CLliv that is total retrorsine CL
fCLliv_m <- 100*CLliv_m/CLtot_m
fCLliv_r <- 100*CLliv_r/CLtot_r



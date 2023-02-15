###############################################################################
# This script 
# - performs reverse dosimetry
# Author:  Anja Lehmann                                                       
# Date: 2023-02-14
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
source("DrugDatabase.R")
source("SpeciesDatabase.R")
source("ParameterVector.R")
source("PartitionCoefficients.R")
source("Pred.R")


# Load packages
library(readr)
library(readxl)
library(crayon)
library(RxODE)
library(ggplot2)
library(tidyr)

# Define color-blind friendly palette (Okate and Ito 2002)
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", 
                      "#CC79A7", "#999999", "#F0E442") 

###############################################################################
# PLOT IN VITRO CONCENTRATION RESPOSNE                                                      
###############################################################################

# read file
dataset <- read_excel("Data_Cytotoxicity Assay_Mouse_Rat.xlsx", 
                      sheet = "DataAnalysis_mouse",
                      col_names = c("ID",
                                    "Viability",
                                    "Date"),skip=1)
View(dataset)

# separate column ID
dataset <- dataset %>% separate('ID', 
                                c("Compound", "Species","Concentration", "Unit",
                                  "BiologicalReplicate","TechnicalReplicate"), 
                                "_")
dataset$Concentration = as.numeric(as.character(dataset$Concentration))

#subset dataframe
dataset <- subset(dataset, BiologicalReplicate != 1)

#-- PLOT IN VITRO DATA

ylimits      <- c(0,110)
ybreaks      <- seq(0,100,by = 20)
xlimits      <- c(0,1000)

plot_invitro <-
  ggplot(data=dataset,aes(x=Concentration, y=Viability))+
  geom_point(mapping=aes(color=BiologicalReplicate,fill=BiologicalReplicate),
             shape=21,size=6,stroke=1.5, alpha=0.4)+
  labs(x="RET (µM)",y="Cell viability \n (% of solvent control)")+
  scale_color_manual(values = c("#009E73","#E69F00","#56B4E9"))+
  scale_fill_manual(values =  c("#009E73","#E69F00","#56B4E9"))+
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 1.5),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=22),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 20),
    axis.title        = element_text(size = 26),
    axis.title.y      = element_text(margin=margin(r=15), angle=90, vjust=0.5),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(22),
    legend.title      = element_text(size =22),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=1
  )+
  guides(col = guide_legend(nrow = 2,title.position = "top"))+
  scale_x_log10(limits = c(0.9,1000),
                labels = function(x) sprintf("%g", x)) +
  scale_y_continuous(
    limits = ylimits,
    breaks = ybreaks)+
  annotate("text", x = Inf-1, y = Inf-1, label = "Mouse",size=10, hjust=1.1, vjust=1.2)+
  annotation_logticks(sides="b")

plot_invitro
ggsave("plot_mouse_invitro.png", bg="white", width = 5, height = 5)

###############################################################################
# INITIALIZATION OF PREDICTION                                              
###############################################################################

#--- DEFINE SPECIESTYPES

# speciestype is a string defining species, sex (only human 15 and 35 years), and age, e.g. "rat10weeks"
# Speciestypes, for which data is available in SpeciesData.R: 
# "rat10weeks",   "rat70weeks", 
# "mouse9weeks",  "mouse70weeks",
# "humannewborn", "human1year", "human5years",  "human10years", 
# "humanm15years","humanf15years","humanm35years","humanf35years"


#Extract speciestypes
speciestypes <- rbind("mouse9weeks","rat10weeks")

#--- LOAD THETA FOR EACH SPECIESTYPE and PRINT PLAUSIBILITY CHECKS 

for(i in seq_along(speciestypes)) {
  speciestype       <- speciestypes[i,]
  drug <- "Retrorsine"
  
  theta             <- GetParameterVector(drug,speciestype,
                                          method_partcoeff="RodgersAndRowland")
  # "Retrorsine" has to be included in DrugDatabase.R
  # method_partcoeff choose between "RodgersAndRowland" / "Schmitt", 
  # NOTE: also change in Pred.R
  
  SpeciesData       <- GetSpeciesData(speciestype)
  DrugData          <- GetDrugData(drug="Retrorsine",speciestype)  
  RodgersAndRowland <- UseMethodRodgersAndRowland(DrugData,SpeciesData)
  Schmitt           <- UseMethodSchmitt(DrugData,SpeciesData)
  
  #Assign variable names xxx_speciestype
  assign(paste("theta",speciestype,sep="_"), theta)
 
  #Print PlausibilityCheck
  print("Plausibilitycheck")
  cat(cyan(speciestype),
      "\n RodgersAndRowland$PlausibilityCheckKA_PR:",  
      yellow(RodgersAndRowland$PlausibilityCheckKA_PR),
      "\n RodgersAndRowland$PlausibilityCheckKA_AP:",  
      yellow(RodgersAndRowland$PlausibilityCheckKA_AP),
      "\n Schmitt$PlausibilityCheckfuC:",            
      yellow(Schmitt$PlausibilityCheckfuC),"\n")
}

#--- CREATE RUNTABLE AND RUNID 
DOSEORIGINAL <- c(0,exp(log(10)*seq(log10(0.01),log10(1000),by=0.01)))
DOSEORIGINALUNIT <- rep("mgperkgbw",length(DOSEORIGINAL))
DOSE <- DOSEORIGINAL/351.399*10^6
DOSEUNIT <- rep("nmolperkgbw",length(DOSE))
ADM <- rep("po",length(DOSE))
CMAXULIV <- rep(NA,length(DOSE))
CMAXULIVUNIT <- rep("µM",length(DOSE))
runtable_reversedosimetry <- data.frame(DOSEORIGINAL,DOSEORIGINALUNIT,DOSE,DOSEUNIT,ADM,CMAXULIV,CMAXULIVUNIT)
View(runtable_reversedosimetry)
  

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

#################################################################################
# PBTK MODEL PREDICTION                                     
#################################################################################

#--- SOLVE ODE SYSTEM USING OPTIMIZED PARAMETERS 
#mouse
runtable_reversedosimetry <- PredRxODE(speciestype="mouse9weeks",
                                     t=seq(0,2,by=0.01),
                                     runtable=runtable_reversedosimetry)

View(runtable_reversedosimetry)

################################################################################  
# REVERSE DOSIMETRY - COMPARISON OF Cmaxuliv and Cuinvitro
################################################################################

#initialize column
dataset$DOSEORIGINAL <- as.numeric(NA)
dataset$DOSEORIGINALUNIT <- as.character(NA)
dataset$CMAXULIV <- as.numeric(NA)

# identify Cmax,u,liv that is closest to Cu,invitro
for(row1 in 1:nrow(dataset)){
  
  row2 <- which.min(abs(runtable_reversedosimetry$CMAXULIV - dataset[row1,]$Concentration))
  dataset[row1,]$DOSEORIGINAL <- runtable_reversedosimetry[row2,]$DOSEORIGINAL
  dataset[row1,]$DOSEORIGINALUNIT <- runtable_reversedosimetry[row2,]$DOSEORIGINALUNIT
  dataset[row1,]$CMAXULIV <- runtable_reversedosimetry[row2,]$CMAXULIV
  
}

#add colum response
dataset$Response <- 100 -dataset$Viability

###############################################################################
# PLOT IN VIVO DOSE RESPONSE                                                      
###############################################################################

ylimits      <- c(-2,110)
ybreaks      <- seq(0,100,by = 25)

plot_invivo <-
  ggplot(data=dataset,aes(x=DOSEORIGINAL, y=Viability))+
  geom_point(mapping=aes(color=BiologicalReplicate,fill=BiologicalReplicate),shape=21,size=6,stroke=1.5, alpha=0.4)+
  labs(x="RET dose (mg/kg bw)",y="Liver viability (%)")+
  scale_color_manual(values = c("#009E73","#E69F00","#56B4E9"))+
  scale_fill_manual(values =  c("#009E73","#E69F00","#56B4E9"))+
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 1.5),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=22),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 20),
    axis.title        = element_text(size = 26),
    axis.title.y      = element_text(margin=margin(r=15), angle=90, vjust=0.5),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(22),
    legend.title      = element_text(size =22),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=1
  )+
  guides(col = guide_legend(nrow = 2,title.position = "top"))+
  scale_x_log10(limits = c(0.1,1000),
                labels = function(x) sprintf("%g", x)) +
  scale_y_continuous(
    limits = ylimits,
    breaks = ybreaks)+
  annotate("text", x = Inf-1, y = Inf-1, label = "Mouse",size=10, hjust=1.1, vjust=1.2)+
  annotation_logticks(sides="b")+
  geom_segment(aes(x = 21.6, y = 95, xend = 21.6, yend = 0),size=1,arrow = arrow(),color="darkgray")+
  geom_hline(yintercept=95, linetype="dashed", color = "darkgray", size=1)+
  geom_segment(aes(x = 21.6, y = 95, xend = 93.9, yend = 95),size=6)
  
plot_invivo
ggsave("plot_mouse_invivo.png", bg="white", width = 5, height = 5)


###############################################################################
# EXPORT IN VIVO DOSE RESPONSE DATA                                                      
###############################################################################

# rename columns
colnames(dataset)[which(names(dataset) == "Viability")] <- "Liver integrity (%)"
colnames(dataset)[which(names(dataset) == "DOSEORIGINAL")] <- "RET dose (mg/kg bw)"

write.table(x=dataset,
            file= "dataset_mouse.txt",
            row.names = FALSE,
            sep="\t")

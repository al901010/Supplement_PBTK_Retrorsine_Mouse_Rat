
###############################################################################
# Author: Anja Lehmann
# Date: 2023-12-02
###############################################################################

#--- INITIALIZE ---------------------------------------------------------------

# Remove all objects from current workspace
rm(list=ls())

# Reset graphics
graphics.off()
cat("\014")

# Set working directory to current path
PATH <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(PATH)

#load dependencies
source("MyFunctions.R")

# install packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load("rstudioapi","ggplot2","readxl","readr","tidyr","stringr",
               "latex2exp","FME","nlmrt","coda")

# Load packages
library(ggplot2)
library(readxl)
library(readr)
library(tidyr)
library(stringr)
library(latex2exp)
library(FME)
library(dfoptim)
library(nlmrt)
library(coda)
library(RxODE)
library(cowplot)
library(bayestestR)


#functions
# for pairs plot
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor)
}

#mode
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
# define colour palette
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00",
                      "#CC79A7", "#999999") #"#F0E442"

#set seed
set.seed(10)

###############################################################################
# LOAD AND PROCESS DATA 
###############################################################################

# set WD to parent folder
setwd('..')

# read file
dataset <- read_excel("Data_Liver Microsomal Assay_Mouse_Rat.xlsx", 
                      sheet = "R data analysis",
                      col_names = c("ID","PeakArea",
                                    "CalculatedConcentration",
                                    "NormalizedConcentration",
                                    "ActualConcentration"), skip=1)
View(dataset)

# reset WD
setwd(PATH)


# separate column ID
dataset <- dataset %>% separate('ID', 
                                c("Compound", "Concentration","ConcentrationUnit", 
                                  "Species", "Assay", "Time", "TimeUnit",
                                  "Replicate"))

dataset$Concentration = as.numeric(as.character(dataset$Concentration))
dataset$Concentration = as.character(dataset$Concentration)
dataset$Time = as.numeric(as.character(dataset$Time))

# Log transform data
dataset$NormalizedConcentration <- TransData(dataset$NormalizedConcentration)



###############################################################################
# FIT AND MCMC
###############################################################################

# filter dataset for time points
dataset <- subset(dataset, dataset$Time == 0 | dataset$Time == 8 | 
                    dataset$Time == 15| dataset$Time == 20  | 
                    dataset$Time == 30| dataset$Time == 40| 
                    dataset$Time == 50| dataset$Time == 60)

# filter dataset for species
dataset <- subset(dataset, dataset$Species == "mouse")


# filter dataset for livermicrosomes and time points
dataset <- subset(dataset, dataset$Assay == "livermicrosomes")

#Create new colums
dataset$ConcentrationReplicate   <- NA
dataset$ConcentrationReplicate   <- paste(dataset$Concentration,
                                          dataset$Replicate,sep="")
dataset$Fit                      <- NA
dataset$InhibConc                <- NA
dataset$MaxDepletionRate         <- NA
dataset$MichaelisMentenConstant  <- NA
dataset$ModelSD                  <- NA

#--- COMPILE ODE SYSTEM RxODE -------------------------------------------------

mModel <-RxODE({
k_1   = Vmax/Km*(1-(1/(1+Km)))*(1-(P_1/(EC50+P_1)));
k_15  = Vmax/Km*(1-(15/(15+Km)))*(1-(P_15/(EC50+P_15)));
k_50  = Vmax/Km*(1-(50/(50+Km)))*(1-(P_50/(EC50+P_50)));
k_200 = Vmax/Km*(1-(200/(200+Km)))*(1-(P_200/(EC50+P_200)));
d/dt (S_1)   = -k_1*S_1;
d/dt (S_15)  = -k_15*S_15;
d/dt (S_50)  = -k_50*S_50;
d/dt (S_200) = -k_200*S_200;
d/dt (P_1)   = k_1*S_1;
d/dt (P_15)  = k_15*S_15;
d/dt (P_50)  = k_50*S_50;
d/dt (P_200) = k_200*S_200;
})

#--- DEFINE INITIAL ESTIMATES -------------------------------------------------

est0 <- c(sd=0.1,EC50=1,Vmax=1,Km=40)

#-- DEFINE BOUNDARIES FOR MLE -------------------------------------------------

lboundary <- c(1e-10,1e-10,0,1e-10)
uboundary <- c(Inf,Inf,Inf,Inf)

#--- MLE ----------------------------------------------------------------------

# out_hjkb   <- hjkb(par=est0,fn=GetMinus2LL,
#                    lower=lboundary,
#                    upper=uboundary,
#                    control=list(),
#                    dataset=dataset,
#                    t=c(0,8,15,20,30,40,50,60))
# 
# out_hjkb
# save(out_hjkb, file="out_hjkb.rda")
load(file = "out_hjkb.rda")

est_hat <- out_hjkb$par

#--- MCMC ---------------------------------------------------------------------

# numberiter <- 5000  
# burnin     <- 0.5*numberiter
# update     <- 0.1*numberiter
# dr         <- 3
# dispersion <- sample(seq(0.7,1.3,by=0.1),length(est_hat),replace=TRUE)
# dispersion2 <- sample(seq(0.7,1.3,by=0.1),length(est_hat),replace=TRUE)
# 
# 
# out_MCMC1 <- modMCMC(f=GetMinus2LL,
#                      p=est_hat,   # initial values for parameters
#                      jump=NULL,          # SD of proposal normal distribution; if NULL 10% of p
#                                          # it can be efficient to use covar from model fit 
#                      lower=lboundary,
#                      upper=uboundary,
#                      prior=NULL,         # -2log(parameter prior probability); NULL: non-informative prior, all parameters are equally likely 
#                      var0=NULL,          # NULL: it is assumned that model variance is 1 and the return element from f is -2logL
#                      wvar0=NULL,         # "weight" for initial model variance, NULL: error variances are assumed to be fixed
#                      n0=NULL,            # parameter used for weihghing initial model variance, NULL: n0=wvar0*n
#                      niter=numberiter,   # no of iterations for the MCMC
#                      updatecov=update,   # setting updatecov smaller than niter will trigger adaptive MH, 
#                                          # proposal distribution is only updated during burnin when burninlenth is positive 
#                      burninlength=burnin,# no of initial iterations to be removed, about 10% of niter
#                      ntrydr=dr,          # max no of tries for delayed rejection procedure  
#                      dataset=dataset,
#                      t=c(0,8,15,20,30,40,50,60))
#    
# # generate second chain for Gelman-Rubin convergence diagnostic
# out_MCMC2 <- modMCMC(f=GetMinus2LL,
#                      p=dispersion*est_hat, 
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
#                      dataset=dataset,
#                      t=c(0,8,15,20,30,40,50,60))
#    
# # generate third chain for Gelman-Rubin convergence diagnostic
# out_MCMC3 <- modMCMC(f=GetMinus2LL,
#                       p=dispersion2*est_hat, 
#                       jump=NULL,         
#                       lower=lboundary,
#                       upper=uboundary,
#                       prior=NULL,        
#                       var0=NULL,         
#                       wvar0=NULL,        
#                       n0=NULL,           
#                       niter=numberiter,     
#                       updatecov=update,
#                       burninlength=burnin,  
#                       ntrydr=dr,          
#                       dataset=dataset,
#                       t=c(0,8,15,20,30,40,50,60))
# 
# save(out_MCMC1, file="out_MCMC1.rda")
# save(out_MCMC2, file="out_MCMC2.rda")
# save(out_MCMC3, file="out_MCMC3.rda")

load(file = "out_MCMC1.rda")
load(file = "out_MCMC2.rda")
load(file = "out_MCMC3.rda")

# join 3 chains
out_MCMC <- rbind(out_MCMC1$pars,out_MCMC2$pars,out_MCMC3$pars)

#assign(paste("out_MCMC",species,sep="_"), out_MCMC)

# MCMC output can be used as functions from the coda package
MC1 <- as.mcmc(out_MCMC1$pars)
MC2 <- as.mcmc(out_MCMC2$pars)
MC3 <- as.mcmc(out_MCMC3$pars)
   
#assign(paste("MC",species,sep="_"), MC)

# Gelman-Rubin convergence dignostc
combinedchains <- mcmc.list(MC1,MC2,MC3)
gelman.plot(combinedchains)
gelmandiag     <- gelman.diag(combinedchains, confidence = 0.95, 
                              autoburnin = TRUE, multivariate = TRUE)
gelmandiag

MCMC_plot <- subset(out_MCMC, select = c(EC50,Vmax,Km))

pairs(MCMC_plot, diag.panel = panel.hist,upper.panel = panel.cor,
      panel=panel.smooth,cex.labels = 1.5, font.labels = 2, cex = 1.5,pch = 21)


# Summary of posterior distribution

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



#--- INCLUDE PREDICTIONS INTO DATASET -----------------------------------------
dataset_pred <- GetPred(est=est_hat,dataset=dataset,t=c(0,8,15,20,30,40,50,60))
dataset$Fit  <- dataset_pred$Fit


#--- INCLUDE PARAMETERS INTO DATASET ------------------------------------------
#mouse
dataset$ModelSD[which(dataset$Species == "mouse")]   <- est_hat["sd"]
dataset$InhibConc[which(dataset$Species == "mouse")] <- est_hat["EC50"]
dataset$MaxDepletionRate[which(dataset$Species == "mouse")] <- est_hat["Vmax"]
dataset$MichaelisMentenConstant[which(dataset$Species=="mouse")]<-est_hat["Km"]

#--- INCLUDE ADDITIONAL TIMEPOINTS INTO DATASET -------------------------------
new_rows <- GetPred2(est=est_hat,dataset=dataset,t=seq(0,60,by=0.1))
dataset  <- rbind(dataset, new_rows)
     

# 
# 
# #--- DRAW SAMPLES FROM PREDICTIVE POSTERIOR DISTRIBUTION --------------------
# # Initialize
# dataset_ppd <- data.frame(matrix(NA,nrow = 1, ncol = ncol(dataset)))
# colnames(dataset_ppd) <- colnames(dataset)
# 
# for(z in 1:nrow(sample_params)) {
#   
#   print(z)
#    
#   dataset_ppd_temp <- GetPred3(est=sample_params[z,],dataset=dataset,t=seq(0,60,by=0.1))
#   # join datasets
#   dataset_ppd <- rbind(dataset_ppd,dataset_ppd_temp)
# 
# }
# 
# 
# # add 0 in Concentration column to order them
# dataset_ppd$Concentration <- str_pad(dataset_ppd$Concentration, 3, pad = "0")
# 
# 
# #--- PREPARE DATA FOR PLOTTING ----------------------------------------------
# dataset$ConcentrationReplicate <- str_pad(dataset$ConcentrationReplicate, 4, pad = "0")
# dataset$Concentration          <- as.character(dataset$Concentration)
# dataset$Concentration          <- str_pad(dataset$Concentration, 3, pad = "0")
# 
# #Detrans observed and predicited data and simulated data
# dataset$NormalizedConcentration <- DeTransData(dataset$NormalizedConcentration) 
# dataset$Fit <- DeTransData(dataset$Fit)
# dataset_ppd$Fit <- DeTransData(dataset_ppd$Fit)
# 
# 
# 
# save(dataset, file="dataset.rda")
# save(dataset_ppd, file="dataset_ppd.rda")

load(file = "dataset.rda")
load(file = "dataset_ppd.rda")


#--- PLOT DATA AND CI OF SIMULATED DATA  --------------------------------------

dataset$NormalizedConcentration <- dataset$NormalizedConcentration*100
dataset$Fit <- dataset$Fit*100
dataset_ppd$Fit <- dataset_ppd$Fit*100

dataset_long <- gather(dataset, Type, Y, "NormalizedConcentration","Fit", 
                       factor_key=TRUE, na.rm = TRUE)
datasetone   <- subset(dataset_long, Type =="NormalizedConcentration")
datasettwo   <- subset(dataset_long, Type =="Fit")
aest         <- aes(x=Time, y=Y)
title        <- ""

xlimits      <- c(0,60)
xbreaks      <- seq(0,60,by = 10)
ylimits      <- c(30,110)


plot <-
  ggplot(data=datasetone,mapping=aes(color=Concentration, fill= Concentration,
                                     shape=Replicate))+
  stat_summary(data=dataset_ppd,
               aes(x=as.numeric(as.character(Time)), 
                   y=as.numeric(as.character(Fit))), 
               geom= "ribbon",
               fun.min = function(x) quantile(x, 0.05),
               fun.max = function(x) quantile(x, 0.95),
               alpha=0.2,
               colour=NA)+
  stat_summary(data=dataset_ppd,
               aes(x=as.numeric(as.character(Time)), 
                   y=as.numeric(as.character(Fit))), 
               geom= "line",
               fun = median, size=1.5)+
  geom_point(aest,size=5, stroke=1.5,alpha=0.6)+
  labs(x="Time (min)",y="RET (%)")+
  scale_color_manual(values = palette_OkabeIto)+
  scale_fill_manual(values =  palette_OkabeIto)+
  scale_shape_manual(values=c(19,17))+
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                     size = 1.5),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=26),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    axis.title.y      = element_text(angle=90),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=26),
    legend.title      = element_text(size =26),
    legend.position   = "none", #"bottom",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66
  )+
  guides(col = guide_legend(nrow = 1,title.position = "top"))+
  scale_y_log10(
     limits = ylimits)+
  annotate("text", x = Inf-1, y = Inf-1, label = "Mouse  ",size=12, hjust=1.2, 
           vjust=1.4)+
  annotation_logticks(sides="l",short = unit(2.5,"mm"),mid = unit(5,"mm"))

plot
ggsave("plotlivermicrosomes_mouse.png", bg="white", width = 9, height = 9)



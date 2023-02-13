
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
                      "#CC79A7", "#999999") 

#set seed
set.seed(10)

###############################################################################
# LOAD AND PROCESS DATA 
###############################################################################

# set WD to parent folder
setwd('..')

# read file
dataset <- read_excel("Data_Medium Loss Assay_Mouse_Rat.xlsx", 
                      sheet = "rat123",
                      col_names = c("ID",
                                    "Concentration",
                                    "NormalizedConcentration"),skip=1)
View(dataset)

# reset WD
setwd(PATH)


# separate column ID
dataset <- dataset %>% separate('ID', 
                                c("Compound", "Concentration","ConcentrationUnit", 
                                  "Condition", "Cells","Time", "TimeUnit", 
                                  "BiologicalReplicate","TechnicalReplicate"), "_")

dataset$Time = as.numeric(as.character(dataset$Time))

#subset dataset
dataset <- subset(dataset, Time==0.0000|
                           Time==0.0417|
                           Time==0.0833|
                           Time==0.1670|
                           Time==0.5000|
                           Time==1.0000)

# Log transform data
dataset$NormalizedConcentration <- TransData(dataset$NormalizedConcentration)


###############################################################################
# FIT AND MCMC
###############################################################################

#Create new colums
dataset$ConcentrationReplicate   <- NA
dataset$ConcentrationReplicate   <- paste(dataset$Concentration,
                                          dataset$Replicate,sep="")
dataset$Fit                      <- NA
dataset$DepletionRate1           <- NA
dataset$ModelSD                  <- NA

#--- COMPILE ODE SYSTEM RxODE -------------------------------------------------

ode <- "
d/dt (S_w) = -k_w*S_w;
d/dt (S_c) = -k_c*S_c;

"

mMonoexponential <- RxODE(model=ode,modName="mMonoexponential")


#--- DEFINE INITIAL ESTIMATES -----------------------------------------------

est0 <- c(sd_w=0.1,k_w=1,
          sd_c=0.1,k_c=1)


#-- DEFINE BOUNDARIES FOR MLE -----------------------------------------------
lboundary <- c(1e-10,0,
               1e-10,0)
uboundary <- c(Inf,Inf,
               Inf,Inf)

#--- MLE --------------------------------------------------------------------
# out_hjkb   <- hjkb(par=est0,fn=GetMinus2LL,
#                    lower=lboundary,
#                    upper=uboundary,
#                    control=list(),
#                    dataset=dataset,
#                    t=c(0,0.0417,0.0833,0.1670,0.5000,1))
# 
# out_hjkb
# save(out_hjkb, file="out_hjkb.rda")

load(file = "out_hjkb.rda")


est_hat <- out_hjkb$par

minus2LL <- out_hjkb$value
minus2LL
AIC <- 2*length(est0) + minus2LL
AIC


#--- MCMC -------------------------------------------------------------------
numberiter <- 5000 
burnin     <- 0.5*numberiter
update     <- 0.1*numberiter
dr         <- 3
dispersion <- sample(seq(0.7,1.3,by=0.1),length(est_hat),replace=TRUE)
dispersion2 <- sample(seq(0.7,1.3,by=0.1),length(est_hat),replace=TRUE)


# out_MCMC1 <- modMCMC(f=GetMinus2LL,
#                      p=est_hat,   # initial values for parameters
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
#                      dataset=dataset,
#                      t=c(0,0.0417,0.0833,0.1670,0.5000,1))
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
#                      t=c(0,0.0417,0.0833,0.1670,0.5000,1))
# 
# # generate third chain for Gelman-Rubin convergence diagnostic
# out_MCMC3 <- modMCMC(f=GetMinus2LL,
#                      p=dispersion2*est_hat, 
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
#                      t=c(0,0.0417,0.0833,0.1670,0.5000,1))
# 
# save(out_MCMC1, file="out_MCMC1.rda")
# save(out_MCMC2, file="out_MCMC2.rda")
# save(out_MCMC3, file="out_MCMC3.rda")

load(file = "out_MCMC1.rda")
load(file = "out_MCMC2.rda")
load(file = "out_MCMC3.rda")

# join 3 chains
out_MCMC <- rbind(out_MCMC1$pars,out_MCMC2$pars,out_MCMC3$pars)

# MCMC output can be used as functions from the coda package
MC1 <- as.mcmc(out_MCMC1$pars)
MC2 <- as.mcmc(out_MCMC2$pars)
MC3 <- as.mcmc(out_MCMC3$pars)

# Gelman-Rubin convergence dignostc
combinedchains <- mcmc.list(MC1,MC2,MC3)
gelman.plot(combinedchains)
gelmandiag     <- gelman.diag(combinedchains, confidence = 0.95, 
                              autoburnin = TRUE, multivariate = TRUE)
gelmandiag

# pairs plot
MCMC_plot <- subset(out_MCMC, select = c(k_w,k_c))

MCMC_new   <- MCMC_plot[seq(1, NROW(MCMC_plot), 10), ]
pairs(MCMC_new, diag.panel = panel.hist,upper.panel = panel.cor,
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


#--- INCLUDE PREDICTIONS INTO DATASET ---------------------------------------
dataset_pred <- GetPred(est=est_hat,dataset=dataset,t=c(0,0.0417,0.0833,0.1670,
                                                        0.5000,1))
dataset$Fit  <- dataset_pred$Fit


#--- INCLUDE PARAMETERS INTO DATASET ------------------------------------------
#warm
dataset$ModelSD[which(dataset$Condition == "warm")]        <- est_hat["sd_w"]
dataset$DepletionRate1[which(dataset$Condition == "warm")]  <- est_hat["k_w"]

#cold
dataset$ModelSD[which(dataset$Condition == "cold")]        <- est_hat["sd_c"]
dataset$DepletionRate1[which(dataset$Condition == "cold")]  <- est_hat["k_c"]


#--- INCLUDE ADDITIONAL TIMEPOINTS INTO DATASET -------------------------------
new_rows <- GetPred2(est=est_hat,dataset=dataset,t=seq(0,1,by=0.01))
dataset  <- rbind(dataset, new_rows)


# #--- DRAW SAMPLES FROM COMBINED MCMC CHAINS -----------------------------------
# sample_params <- out_MCMC[sample(nrow(out_MCMC),size=1000,replace = TRUE),]
# 
# # pairs plot of sample_params
# sample_params_plot <- subset(sample_params, select = c(k_w,k_c))
# 
# sample_params_plot   <- sample_params_plot[seq(1, NROW(sample_params_plot), 1), ]
# pairs(sample_params_plot, diag.panel = panel.hist,upper.panel = panel.cor,panel=panel.smooth,
#       cex.labels = 1.5, font.labels = 2, cex = 1.5,pch = 21)
# 
# 
# #--- DRAW SAMPLES FROM PREDICTIVE POSTERIOR DISTRIBUTION ------------------
# # Initialize
# dataset_ppd <- data.frame(matrix(NA,nrow = 1, ncol = ncol(dataset)))
# colnames(dataset_ppd) <- colnames(dataset)
# 
# for(z in 1:nrow(sample_params)) {
# 
#   print(z)
# 
#   dataset_ppd_temp <- GetPred3(est=sample_params[z,],dataset=dataset,t=seq(0,1,by=0.001))
#   # join datasets
#   dataset_ppd <- rbind(dataset_ppd,dataset_ppd_temp)
# 
# }
# 
# #--- PREPARE DATA FOR PLOTTING ------------------------------------------------
# 
# #Detrans observed and predicited data and simulated data
# dataset$NormalizedConcentration <- DeTransData(dataset$NormalizedConcentration) 
# dataset$Fit <- DeTransData(dataset$Fit)
# dataset_ppd$Fit <- DeTransData(dataset_ppd$Fit)
# 
# save(dataset, file="dataset.rda")
# save(dataset_ppd, file="dataset_ppd.rda")

load(file = "dataset.rda")
load(file = "dataset_ppd.rda")



#--- PLOT DATA AND CI OF SIMULATED DATA  --------------------------------------
library(ggnewscale)

dataset$NormalizedConcentration <- dataset$NormalizedConcentration*100
dataset$Fit <- dataset$Fit*100
dataset$Time <- dataset$Time*60
dataset_ppd$Fit <- dataset_ppd$Fit*100
dataset_ppd$Time <- dataset_ppd$Time*60



#new column 
dataset$ConditionBiologicalTechnicalReplicate <- paste(dataset$Condition,dataset$BiologicalReplicate,dataset$TechnicalReplicate)


dataset_long <- gather(dataset, Type, Y, "NormalizedConcentration","Fit", 
                       factor_key=TRUE, na.rm = TRUE)
datasetone   <- subset(dataset_long, Type =="NormalizedConcentration")
datasettwo   <- subset(dataset_long, Type =="Fit")
aest         <- aes(x=Time, y=Y)
xlabel       <- "Time (min)"
ylabel       <- "RET (%)"
title        <- ""
ylimits      <- c(68,110)

xlimits      <- c(0,60)
xbreaks      <- seq(0,60,by = 20)


plot2 <-
  ggplot()+
  stat_summary(data=dataset_ppd,
               aes(x=as.numeric(as.character(Time)), y=as.numeric(as.character(Fit)),fill=Condition),
               geom= "ribbon",
               fun.min = function(x) quantile(x, 0.05),
               fun.max = function(x) quantile(x, 0.95),
               alpha=0.5,
               colour=NA)+
  stat_summary(data=dataset_ppd,
               aes(x=as.numeric(as.character(Time)), y=as.numeric(as.character(Fit)),color=Condition),
               geom= "line",
               fun = median, 
               size=1.2)+
  scale_color_manual(values = c("#68abb8","#c1766f"))+
  new_scale_colour() +
  geom_point(data=datasetone,
             aes(x=Time, y=Y,color=ConditionBiologicalTechnicalReplicate,shape=TechnicalReplicate),
             alpha=0.8, 
             stroke=1.5,
             size=5)+
  labs(x=xlabel,y=ylabel)+
  scale_color_manual(values = c("#a8dbd9","#a8dbd9","#a8dbd9","#68abb8","#68abb8","#68abb8","#3b738f","#3b738f","#3b738f",
                                "#e0c2a2","#e0c2a2","#e0c2a2","#c1766f","#c1766f","#c1766f","#813753","#813753","#813753"))+
  scale_fill_manual(values = c("#a8dbd9","#e0c2a2"))+
  scale_shape_manual(values =  c(19,17,15))+
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 1.5),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=26),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 26),
    axis.title        = element_text(size = 28),
    axis.title.y      = element_text(margin=margin(r=15), angle=90, vjust=0.5),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=26),
    legend.title      = element_text(size =26),
    legend.position   = "none",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66
  )+
  guides(col = guide_legend(nrow = 1,title.position = "top"))+
  scale_y_log10(
    limits = ylimits
  )+
  scale_x_continuous(
    limits = xlimits,
    breaks = xbreaks,
    labels = function(x) sprintf("%g", x))+
  annotate("text", x = Inf-1, y = Inf-1, label = "Rat  ",size=12, hjust=1.1, vjust=1.3)+
  annotation_logticks(sides="l",short = unit(2.5,"mm"),mid = unit(5,"mm"))

plot2
ggsave("plot2_rat.png", bg="white", width = 9, height = 9)


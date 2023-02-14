
###############################################################################
# Author: Anja Lehmann
# Date: 2023-02-11
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
dataset <- read_excel("Data_Cytotoxicity Assay_Mouse_Rat.xlsx", 
                      sheet = "DataAnalysis_rat",
                      col_names = c("ID",
                                    "Viability",
                                    "Date"),skip=1)
View(dataset)

# reset WD
setwd(PATH)

# separate column ID
dataset <- dataset %>% separate('ID', 
                                c("Compound", "Species","Concentration", "Unit",
                                  "BiologicalReplicate","TechnicalReplicate"), 
                                "_")

dataset$Concentration = as.numeric(as.character(dataset$Concentration))

# Log transform data
dataset$Concentration <- TransData(dataset$Concentration)

###############################################################################
# FIT AND MCMC
###############################################################################

#Create new colums
dataset$Fit                      <- NA

# Emax model
sigmoidal <- function(S) 100/(1+10^(S-log(IC50)))

# initial estimates

est0 <- c(sd=1,
          IC50=1)

# boundaries
lboundary <- c(1e-10,
               0)
uboundary <- c(Inf,
               Inf)

# extract concentrations from dataset
# replicate each= no of experimental replicates
C <- rep(unique(dataset$Concentration),each=9)
#C <- unique(dataset$Concentration)

# out_hjkb   <- hjkb(par=est0,
#                    fn=GetMinus2LL,
#                    lower=lboundary,
#                    upper=uboundary,
#                    control=list(),
#                    dataset=dataset,
#                    C=C)
# 
# out_hjkb
# save(out_hjkb, file="out_hjkb.rda")
load(file = "out_hjkb.rda")

est_hat <- out_hjkb$par

#minus 2LL and AIC
minus2LL <- out_hjkb$value
minus2LL
AIC <- 2*length(est0) + minus2LL
AIC

# prediction for best-fit parameters
IC50        <- est_hat["IC50"]

best_fit         <- sigmoidal(S=C)
dataset_pred     <- dataset
dataset_pred$Fit <- best_fit

View(dataset_pred)


# prediction for additional timepoints
new_rows <- data.frame(matrix(NA,nrow = NROW(seq(0,1000,0.01)), 
                              ncol = ncol(dataset)))
colnames(new_rows)              <- colnames(dataset)
new_rows$Concentration          <- TransData(seq(0,1000,0.01))
new_rows$Species                <- "mouse"
new_rows$Fit                    <- sigmoidal(S=TransData(seq(0,1000,by=0.01)))

#add new rows to dataset
dataset_pred <- rbind(dataset_pred, new_rows)

#Detrans observed and predicited data 
dataset_pred$Concentration <- DeTransData(dataset_pred$Concentration)



#--- MCMC -------------------------------------------------------------------
numberiter <- 50000 
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
#                      C=C)
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
#                      C=C)
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
#                      C=C)
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

# pairs plot
MCMC_plot <- subset(out_MCMC, select = c(IC50))

MCMC_new   <- MCMC_plot[seq(1, NROW(MCMC_plot), 10), ]
pairs(MCMC_new, diag.panel = panel.hist,upper.panel = panel.cor,panel=panel.smooth,
      cex.labels = 1.5, font.labels = 2, cex = 1.5,pch = 21)


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


# #--- DRAW SAMPLES FROM COMBINED MCMC CHAINS -----------------------------------
# sample_params <- out_MCMC[sample(nrow(out_MCMC),size=1000,replace = TRUE),]
# 
# #--- DRAW SAMPLES FROM PREDICTIVE POSTERIOR DISTRIBUTION ------------------
# # Initialize
# c1 <- exp(log(10)*seq(log10(0.1),log10(1000),by=0.001))
# c2 <- seq(0,1000,by=0.1)
# 
# Concentration         <- TransData(c1)
# dataset_ppd           <- data.frame(matrix(NA,nrow = 1, ncol = 4))
# colnames(dataset_ppd) <- c("Compound","Species","Concentration","Fit")
# 
# Compound            <- rep("Retrorsine",length(Concentration))
# Species             <- rep("mouse",length(Concentration))
# 
# for(z in 1:nrow(sample_params)) {
# 
#   print(z)
#   sd         <- sample_params[z,"sd"]
#   IC50       <- sample_params[z,"IC50"]
#   
#   
#   Fit        <- sigmoidal(S=Concentration) + rnorm(n=length(Concentration),mean=0,sd=sd)
#   dataset_ppd_temp <- cbind(Compound,Species,Concentration,Fit)
#   # join datasets
#   dataset_ppd <- rbind(dataset_ppd,dataset_ppd_temp)
# }
# 
# 
# #Detrans  simulated data
# dataset_ppd$Concentration <- as.numeric(dataset_ppd$Concentration)
#  dataset_ppd$Concentration       <- DeTransData(dataset_ppd$Concentration)
# 
# 
# 
# save(dataset_pred, file="dataset_pred.rda")
# save(dataset_ppd, file="dataset_ppd.rda")

load(file = "dataset_pred.rda")
load(file = "dataset_ppd.rda")

#--- PLOT DATA AND CI OF SIMULATED DATA  --------------------------------------

dataset_long <- gather(dataset_pred, Type, Y, "Viability","Fit", 
                       factor_key=TRUE, na.rm = TRUE)
datasetone   <- subset(dataset_long, Type =="Viability")
datasettwo   <- subset(dataset_long, Type =="Fit")
aest         <- aes(x=Concentration, y=Y)
xlabel       <- "RET (µM)"
ylabel       <- "Cell viability \n (% of solvent control)" 
title        <- ""
ylimits      <- c(0,110)
ybreaks      <- seq(0,100,by = 25)
#xlimits      <- c(0,1000)
# xbreaks      <- seq(0,24,by = 6)

  plot2 <-
    ggplot(data=datasetone,aes(x=Concentration, y=Y))+
    stat_summary(data=dataset_ppd,
                 aes(x=as.numeric(as.character(Concentration)), 
                     y=as.numeric(as.character(Fit))),
                 geom= "ribbon",
                 fun.min = function(x) quantile(x, 0.05),
                 fun.max = function(x) quantile(x, 0.95),
                 alpha=0.2,
                 colour=NA)+
    geom_point(mapping=aes(color=BiologicalReplicate,fill=BiologicalReplicate),
               shape=21,size=6,stroke=1.5, alpha=0.4)+
    stat_summary(data=dataset_ppd,
                 aes(x=as.numeric(as.character(Concentration)), 
                     y=as.numeric(as.character(Fit))),
                 geom= "line",
                 fun = median,
                 size=0.5)+
    labs(x=xlabel,y=ylabel)+
    scale_color_manual(values = c("#009E73","#E69F00","#56B4E9"))+
    scale_fill_manual(values =  c("#009E73","#E69F00","#56B4E9"))+
    #scale_shape_manual(values =  c(22,21,25))+
    theme(
      plot.background   = element_rect(fill = "white"),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",
                                       size = 1.5),
      panel.background  = element_rect(fill="white"),
      plot.title        = element_text(size=26),
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
    geom_hline(yintercept=0.666, linetype="dashed", color = "black", size=1.5)+
    annotate("text", x = 1.5, y = 8, label = "10% DMSO",size=5, hjust=1, vjust=1, 
             color="black")+
    annotate("text", x = Inf-1, y = Inf-1, label = "Rat",size=10, hjust=1.1, 
             vjust=1.2)+
    annotation_logticks(sides="b")
      
plot2
ggsave("plot_rat_sim.png", bg="white", width = 5, height = 5)





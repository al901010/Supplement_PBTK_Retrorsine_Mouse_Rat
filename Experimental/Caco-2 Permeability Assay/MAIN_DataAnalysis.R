
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

# install packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load("rstudioapi","ggplot2","readxl","readr","tidyr","stringr",
               "latex2exp","FME","nlmrt","coda")

# Load packages
library(ggplot2)
library(readxl)
library(readr)

# define colour palette
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00",
                      "#CC79A7", "#999999") #"#F0E442"

#--- LOAD AND PROCESS DATA ---------------------------------------------------- 

# set WD to parent folder
setwd('..')

# read file
dataset <- read_excel("Data_Caco-2 Permeability Assay.xlsx", 
                      sheet = "DataAnalysis",
                      col_names = c("ID",
                                    "Creceiver",
                                    "Unit",
                                    "Percentreceiver"),skip=1)
View(dataset)

# reset WD
setwd(PATH)

# separate column ID
dataset <- dataset %>% separate('ID', 
                                c("Compound", "IncubationSide","Time", "TimeUnit",
                                  "TechnicalReplicate"), "_")

dataset$Time <- as.numeric(as.character(dataset$Time))



#--- PLOT --------------------------------------------------------------------- 
aest         <- aes(x=Time, y=Percentreceiver)
xlabel       <- "Time (h)"
ylabel       <- c("% RET in acceptor \ncompartment")
title        <- ""
ylimits      <- c(0,50)
ybreaks      <- seq(0,50,by = 10)
xlimits      <- c(0,24)
xbreaks      <- seq(0,24,by = 4)

#summary for error bars
library(dplyr)
dataset_summary <- dataset %>%
  group_by(IncubationSide, Time) %>%
  summarise(
    sd = sd(Percentreceiver),
    mean = mean(Percentreceiver)
  )
dataset_summary

plot <-
  ggplot(data=dataset,mapping=aes(color=IncubationSide,fill=IncubationSide,
                                  shape=IncubationSide))+
  geom_point(aest,size=6,stroke=1.5, alpha=0.4)+
  stat_summary(data=dataset,aes(x=Time,y=Percentreceiver,group=IncubationSide,
                                linetype=IncubationSide),
               fun=mean,geom="line", size=1.5)+
  geom_errorbar( data = dataset_summary,aes(x=Time,y=sd,ymin = mean-sd, 
                                            ymax = mean+sd),width=1.2,size=1.5)+  
  labs(x=xlabel,y=ylabel)+
  scale_shape_manual(values = c(21,23))+
  scale_color_manual(values = palette_OkabeIto)+
  scale_fill_manual(values =  palette_OkabeIto)+
  theme(
    plot.background   = element_rect(fill = "white"),
    panel.border      = element_rect(fill = "NA",colour="black",linetype="solid",
                                     size = 1),
    panel.background  = element_rect(fill="white"),
    plot.title        = element_text(size=26),
    axis.line         = element_line(size = 0.5),
    axis.text         = element_text(colour = "black",size = 20),
    axis.title        = element_text(size = 24),
    axis.title.y      = element_text(angle=90),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(2.5,"mm"),
    legend.text       = element_text(size=26),
    legend.title      = element_text(size =26),
    legend.position   = "none", #"bottom",
    legend.background = element_rect(fill = NA),
    aspect.ratio=0.66
  )+
  guides(col = guide_legend(nrow = 2,title.position = "top"))+
  scale_x_continuous(
    limits = xlimits, 
    breaks = xbreaks)+
  scale_y_continuous(
    limits = ylimits,
    breaks = ybreaks)+
  annotate("text", x = Inf-1, y = Inf-1, label = "Caco-2",size=10, hjust=1.1,
           vjust=1.2)

plot
ggsave("plot_Caco2.png", bg="white", width = 6, height = 6)






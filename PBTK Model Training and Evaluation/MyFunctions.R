
###############################################################################
# Additional basic functions            
# Authors: Anja Lehmann, Christoph Hethey                                                       
# Date: 2023-02-13                                            
###############################################################################

###############################################################################
# Function to log transform data
###############################################################################
TransData <- function(x) {
  
  x_hat <- log(x+ 1)
  
  return(x_hat)
  
}

###############################################################################
# Function to reverse log transformation of data
###############################################################################

DeTransData <- function(x_hat) {
  
  x <- exp(x_hat) - 1
  
  return(x)
  
}


###############################################################################
# Functions to create histogram and pairs plot for MCMC results
###############################################################################

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

###############################################################################
# Function to calculate mode 
###############################################################################

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}


###############################################################################
# Function to measure runtime
###############################################################################

tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

toc <- function()
{
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  print(toc - tic)
  invisible(toc)
}

###############################################################################
# function returns head of a list 
###############################################################################

head.list <- function(obj, n = 6L, ...)
{
  stopifnot(length(n) == 1L)
  origN <- n
  n <- if (n < 0L)
    max(length(obj) + n, 0L)
  else min(n, length(obj))
  lapply(obj[seq_len(n)], function(x)
  {
    tryCatch({
      head(x, origN, ...)
    }, error = function(e) {
      x
    })
  })
}



################################################################################
# Residuals
################################################################################


#plot data normal scale 
GetPlotResiduals <- function(dataset,aest,xlabel,ylabel,title) {
  library(ggplot2)
  plot <- 
    ggplot(dataset,aest)+
    geom_point(shape=16,size=2)+
    labs(x=xlabel,y=ylabel)+
    ggtitle(title)+ 
    theme(
      plot.background   = element_rect(fill = NA),
      panel.border      = element_rect(fill = NA,colour="black",linetype="solid",size = 0.5),
      panel.background  = element_rect(fill=NA),
      plot.title        = element_text(size=10, hjust = 0.5),
      axis.line         = element_line(size = 0.5),
      axis.text         = element_text(colour = "black",size = 10),
      axis.title        = element_text(size = 10),
      axis.title.y      = element_text(margin=margin(r=15)),
      axis.title.x      = element_text(margin=margin(t=15)),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(2.5,"mm"),
      legend.text       = element_text(size=10),
      legend.title      = element_text(size =10),
      legend.position   = "bottom",
      legend.background = element_rect(fill = NA),
      aspect.ratio = 0.5
    )+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    scale_y_continuous(
      limits = c(-3,3),
      breaks = seq(-3,3,by = 1))
  
  return(plot)
}


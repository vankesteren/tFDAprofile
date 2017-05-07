# Profiling dataset tFDA helper functions

# Erik-Jan van Kesteren
# Thesis - UMC Utrecht Julius Centre
# License: https://opensource.org/licenses/MIT

# Version information:
# platform       x86_64-w64-mingw32
# arch           x86_64
# os             mingw32
# system         x86_64, mingw32
# status
# major          3
# minor          3.2
# year           2016
# month          10
# day            31
# svn rev        71607
# language       R
# version.string R version 3.3.2 (2016-10-31)
# nickname       Sincere Pumpkin Patch

n_out <- function(x) x*(x+1)/2-x

densityderiv <- function(y, precision = 2^16){
  # Function to calculate kernel density and first and second derivatives in
  # order to estimate first and second inflection points of the tFDA density.
  # With help from:
  # http://stackoverflow.com/questions/12568715/derivative-of-kernel-density
  
  x <- seq(min(y), max(y), length.out = precision)
  f <- density(y, n = precision)
  
  # first derivative
  dy <- diff(f$y)
  dx <- diff(f$x)[1]
  f1 <- ((c(dy, tail(dy, 1)) + c(head(dy, 1), dy))/2)/dx
  
  # second derivative
  dy <- diff(f1)
  dx <- diff(f$x)[1]
  f2 <- ((c(dy, tail(dy, 1)) + c(head(dy, 1), dy))/2)/dx
  
  # calculate points of interest
  inflect <- x[which.min(f1)]
  inflect2 <- x[which.min(f2):length(x)][which.max(f2[which.min(f2):length(x)])]
  
  # plot the three functions
  opt <- par(mfrow=c(3,1), mar = c(1,4,1,1), oma = c(4,0,0,0))
  plot(x, f$y, type = "l", xlim = c(1,1+(3*(inflect2-1))), ylab = "density", 
       axes = F, xlab = "")
  abline(v = c(inflect, inflect2), lty = c(1,3))
  axis(2)
  plot(x, f1, type = "l", xlim = c(1,1+(3*(inflect2-1))), ylab = "first deriv",
       axes = F, xlab = "")
  abline(v = c(inflect, inflect2), lty = c(1,3))
  axis(2)
  plot(x, f2, type = "l", xlim = c(1,1+(3*(inflect2-1))), ylab = "second deriv",
       axes = F, xlab = "")
  abline(v = c(inflect, inflect2), lty = c(1,3))
  axis(2)
  axis(1)
  mtext("tFDA",1,2,TRUE)
  par(opt)
  
  return(c(inflect = inflect, inflect2 = inflect2))
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
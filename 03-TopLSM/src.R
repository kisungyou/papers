# Generate a network of SBM structure
# n : number of nodes
# k : number of clusters
# p_high : within-community connectivity
# p_low  : between-community linkage
require(igraph)

generate_sbm <- function(n, k, p_high, p_low){
  nk    = floor(n/k)
  nsize = c(rep(nk, k-1), round(n-(nk*(k-1))))
  probs = array(p_low, c(k,k))
  diag(probs) = p_high
  
  return(igraph::sample_sbm(n, probs, nsize))
  # obj_matrix = as.matrix(igraph::as_adjacency_matrix(out_igraph))
  # return(obj_matrix)
}

# extract network features into (n x 6) : 6-features
#
# Wikipedia : Network Science 
# https://en.wikipedia.org/wiki/Network_science#Network_properties
# 
# 1. density (igraph:graph density)
# 2. average degree
# 3. average shortest path (igraph::mean_distance)
# 4. diameter
# 5. average clustering coefficient (transitivity)

network_features <- function(list_igraph){
  n = length(list_igraph)
  for (i in 1:n){
    if (!inherits(list_igraph[[i]], "igraph")){
      stop(paste0(i,"-th network is not igraph object"))
    }
  }
  
  output = array(0,c(n, 5))
  for (i in 1:n){
    graph_i = list_igraph[[i]] # select
    
    output[i,1] = igraph::edge_density(graph_i, loops = FALSE)       # 1. density
    output[i,2] = base::mean(igraph::degree(graph_i, loops = FALSE)) # 2. average degree
    output[i,3] = igraph::mean_distance(graph_i, directed = FALSE)   # 3. average shortest path
    output[i,4] = igraph::diameter(graph_i, directed = FALSE)        # 4. diameter
    output[i,5] = igraph::transitivity(graph_i, type="global")   # 5. average clustering coefficient
  }
  colnames(output) = c("density","AD","ASD","diameter","CC")
  return(output)
}

# PERSONAL FUNCTIONS 
# (1) draw_complex
# (2) draw_distmatrix
#
# Use 'vignette("ggplot2-specs")' for ggplot2 specs


library(ggplot2)
library(ggforce)
library(reshape2)

# (1) draw_complex --------------------------------------------------------
#     x,y          : coordinate vectors of points
#     r            : fixed size of radius
#     dotsize      : size of the center dots
#     linewidth    : width of the line segments
#     xlims, ylims : range of displays
#    
# default parameters
# alpha    = 0.1     : facet transparency
# colour   = "black" : color of dots and lines 
# linetype = 0       : no line around the circles
draw_complex <- function(x, y, r=0, dotsize=1, linewidth=1, xlims=NULL, ylims=NULL){
  # data wrangling
  xx = as.vector(x)
  yy = as.vector(y)
  nn = length(xx)
  
  rval = max(0, as.double(r)) # radius
  rr   = rep(rval, nn)
  if (length(yy)!=length(xx)){
    stop("* draw_complex : x, y should be of same length.")
  }
  cc   = rep(1, nn) # fill color
  lval = max(sqrt(.Machine$double.eps), as.double(linewidth))
  
  # plot points and circles circles
  circles <- data.frame(x0=xx, y0=yy, rr=rr, cc=cc)
  g <- ggplot2::ggplot() + 
    ggplot2::geom_point(aes(x=x0, y=y0), data=circles, 
                        size=dotsize, colour="black") # parameter : dotsize for center dots
  if (rval > 0){
    g <- g + ggforce::geom_circle(aes(x0=x0, y0=y0, r=rr, fill=cc), data=circles, 
                                  show.legend = FALSE, linetype=0, alpha=0.1)
  }
  
  # add lines
  x1s = c()
  x2s = c()
  y1s = c()
  y2s = c()
  for (i in 1:(nn-1)){
    for (j in (i+1):nn){
      dist_ij = sqrt(((xx[i]-xx[j])^2) + ((yy[i]-yy[j])^2))
      if (dist_ij/2 <= rval){
        x1s = c(x1s, xx[i])
        x2s = c(x2s, xx[j])
        y1s = c(y1s, yy[i])
        y2s = c(y2s, yy[j])
        
        # df.tmp = data.frame(x1=xx[i], x2=xx[j], y1=yy[i], y2=yy[j])
      }
    }
  }
  if (length(x1s) > 0){
    df.tmp = data.frame(x1=x1s, x2=x2s, y1=y1s, y2=y2s)
    g <- g + geom_segment(data=df.tmp, colour="black", size=lval,
                          aes(x = x1, y = y1, xend = x2, yend = y2))  
  }
  # for (i in 1:(nn-1)){
  #   for (j in (i+1):nn){
  #     dist_ij = sqrt(((xx[i]-xx[j])^2) + ((yy[i]-yy[j])^2))
  #     if (dist_ij/2 <= rval){
  #       df.tmp = data.frame(x1=xx[i], x2=xx[j], y1=yy[i], y2=yy[j])
  #       g <- g + geom_segment(data=df.tmp, colour="black", size=lval,
  #                             aes(x = x1, y = y1, xend = x2, yend = y2))
  #     }
  #   }
  # }
  
  # range
  if (is.null(xlims)){
    xrange = c(min(xx)-rval, max(xx)+rval)
    xwidth = xrange[2]-xrange[1]
    xlims  = c(xrange[1]-0.05*xwidth, xrange[2]+0.05*xwidth)
  }
  if (is.null(ylims)){
    yrange = c(min(yy)-rval, max(yy)+rval)
    ywidth = yrange[2]-yrange[1]
    ylims  = c(yrange[1]-0.05*ywidth, yrange[2]+0.05*ywidth)
  }
  
  # final adjust
  g <- g + 
    coord_fixed(xlim = xlims, ylim=ylims) + 
    theme_void()
}
# #   example
# x = rep(1:3, 3)
# y = rep(1:3, each = 3)
# hey = draw_complex(x,y,r=0.5, dotsize = 2, linewidth = 0.5)
# plot(hey)



# (2) draw_distmatrix -----------------------------------------------------
draw_distmatrix <- function(distmat, maxcol="red", no.legend=TRUE){
  if (inherits(distmat,"dist")){
    distmat = as.matrix(distmat)
  }
  longData  = reshape2::melt(distmat)
  if (no.legend){
    gplot <- ggplot(longData, aes(x=Var2, y=Var1)) + 
      geom_raster(aes(fill=value), show.legend = FALSE) + 
      scale_fill_gradient(low="grey90", high=maxcol) + 
      coord_fixed() + 
      theme_void() + 
      scale_y_reverse()  
  } else {
    gplot <- ggplot(longData, aes(x=Var2, y=Var1)) + 
      geom_raster(aes(fill=value), show.legend = TRUE) + 
      scale_fill_gradient(low="grey90", high=maxcol) + 
      coord_fixed() + 
      theme_void() + 
      scale_y_reverse()
  }
  return(gplot)
}
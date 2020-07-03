
#'cmrSocPlot
#'
#'igraph plotting for CMRnet object
#'
#'@param nets A list of igraph networks output by cmr_igraph
#'@param fixed_locs (TRUE/FALSE). Whether to keep the the locations of the nodes fixed between network windows or allow them to change
#'@param locs A matrix or dataframe giving x and y coordinates of all nodes in the full network. Can be used to provide user-defined layouts to the graph
#'@param dynamic (TRUE/FALSE). If true this function will plot networks for the different network windows sequentially in the same plot window with a delay. If false then this function will plot the networks in a single, multipanel figure
#'@param rows Sets the number of rows in the multipanel figure is dynamic=FALSE
#'@param \dots Extra arguments from \code{\link[igraph]{plot.igraph}} if desired
#'@details This function plots co-capture (social) networks produced by CMRnet. As well as the arguments detailed above, additional arguments from plot.igraph can be used to produce plots desired by the user
#'@return A series of network plots showing the co-capture (social) network from each network window
#'@examples
#'\dontrun{
#'
#'# load in data
#'data(cmrData)
#'
#'# set parameters
#'mindate<-"2010-01-01"
#'maxdate<-"2015-01-01"
#'intwindow<-60
#'netwindow<-12
#'overlap<-0
#'spacewindow <- 0
#'
#'# Dynamic networks
#'netdat <- DynamicNetCreate(
#'  data = cmrData,
#'  intwindow = intwindow,
#'  mindate = mindate,
#'  maxdate = maxdate,
#'  netwindow = netwindow,
#'  overlap = overlap,
#'  spacewindow = spacewindow
#')
#'
#'# run permutations
#'networks<-cmr_igraph(netdat,type="social")
#'cmrSocPlot(nets=networks,fixed_locs=TRUE,dynamic=FALSE,rows=2,vertex.label=NA)
#'}
#'@export

cmrSocPlot<-function(nets,fixed_locs=c(TRUE,FALSE),locs=NULL,dynamic=c(TRUE,FALSE),rows=3, ...){
  L<-length(nets[[1]])

  if(fixed_locs==TRUE&length(locs)==0){
    lo<-igraph::layout_with_fr(nets[[2]])
  }

  if(length(locs)>0){
    lo<-locs
  }

  if(dynamic==FALSE){
    nws<-length(nets[[1]])
    graphics::par(mar=c(1,1,1,1),mfrow=c(rows,ceiling(nws/rows)))
    for(i in 1:nws){
      if(fixed_locs==FALSE&length(locs)==0){
        lo2<-igraph::layout_with_fr(nets[[1]][[i]])
      }
      lo2<-lo[igraph::vertex_attr(nets[[2]])$name%in%igraph::vertex_attr(nets[[1]][[i]])$name,]
      if(nrow(lo2)>0){
        igraph::plot.igraph(nets[[1]][[i]],layout=lo2,main=paste("Network Window",i), ...)
      }
      if(nrow(lo2)==0){
        graphics::plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",bty="n",main=paste("Network Window",i))
      }
    }
  }
  if(dynamic==TRUE){
    nws<-length(nets[[1]])
    graphics::par(mar=c(1,1,1,1),mfrow=c(1,1))
    for(i in 1:nws){
      if(fixed_locs==FALSE&length(locs)==0){
        lo2<-igraph::layout_with_fr(nets[[1]][[i]])
      }
      lo2<-lo[igraph::vertex_attr(nets[[2]])$name%in%igraph::vertex_attr(nets[[1]][[i]])$name,]
      if(nrow(lo2)>0){
        igraph::plot.igraph(nets[[1]][[i]],layout=lo2,main=paste("Network Window",i), ...)
      }
      if(nrow(lo2)==0){
        graphics::plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",bty="n",main=paste("Network Window",i))
      }
      Sys.sleep(1)
    }
  }
}

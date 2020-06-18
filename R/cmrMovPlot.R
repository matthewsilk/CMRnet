
#'cmrMovPlot
#'#'
#'igraph plotting for CMRnet object
#'
#'@param nets A list of igraph networks output by cmr_igraph
#'@param fixed_locs (TRUE/FALSE). Whether to keep the the locations of the nodes fixed between network windows or allow them to change
#'@param locs A matrix or dataframe giving x and y coordinates of all nodes in the full network. Can be used to provide user-defined layouts to the graph - helpful for plotting locations in geographic space
#'@param dynamic (TRUE/FALSE). If true this function will plot networks for the different network windows sequentially in the same plot window with a delay. If false then this function will plot the networks in a single, multipanel figure
#'@param rows Sets the number of rows in the multipanel figure is dynamic=FALSE
#'@details This function plots movement networks produced by CMRnet. As well as the arguments detailed above, additional arguments from plot.igraph can be used to produce plots desired by the user
#'#'@return A series of network plots showing the movement network from each network window
#' #'@examples
#'\dontrun{
#'# example without multiple layers ####
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
#'
#'# create network
#'movenetdat<-MoveNetCreate(data=cmrData,
#'intwindow=intwindow,
#'mindate=mindate,
#'maxdate=maxdate,
#'netwindow=netwindow,
#'overlap=overlap,
#'nextonly=TRUE)
#'
#'# convert to networks
#'networks<-cmr_igraph(movenetdat,type="movement")
#'cmrSocPlot(nets=networks,fixed_locs=TRUE,dynamic=FALSE,rows=2)

cmrMovPlot<-function(nets,fixed_locs=c(TRUE,FALSE),locs=NULL,dynamic=c(TRUE,FALSE),rows=3, ...){
  L<-length(nets[[1]])

  if(fixed_locs==TRUE&length(locs)==0){
    lo<-igraph::layout_with_fr(nets[[2]])
  }

  if(length(locs)>0){
    lo<-locs
  }

  if(dynamic==FALSE){
    nws<-length(nets[[1]])
    par(mar=c(1,1,1,1),mfrow=c(rows,ceiling(nws/rows)))
    for(i in 1:nws){
      if(fixed_locs==FALSE&length(locs)==0){
        lo2<-igraph::layout_with_fr(nets[[1]][[i]])
      }
      lo2<-lo[igraph::vertex_attr(nets[[2]])$name%in%igraph::vertex_attr(nets[[1]][[i]])$name,]
      igraph::plot.igraph(nets[[1]][[i]],layout=lo2,main=paste("Network Window",i), ...)
    }
  }
  if(dynamic==TRUE){
    nws<-length(nets[[1]])
    par(mar=c(1,1,1,1),mfrow=c(1,1))
    for(i in 1:nws){
      if(fixed_locs==FALSE&length(locs)==0){
        lo2<-igraph::layout_with_fr(nets[[1]][[i]])
      }
      lo2<-lo[igraph::vertex_attr(nets[[2]])$name%in%igraph::vertex_attr(nets[[1]][[i]])$name,]
      igraph::plot.igraph(nets[[1]][[i]],layout=lo2,main=paste("Network Window",i), ...)
      Sys.sleep(1)
    }
  }
}
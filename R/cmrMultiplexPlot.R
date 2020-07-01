
#'cmrMultiplexPlot
#'#'
#'igraph plotting for CMRnet object
#'
#'@param nets A list of multiplex igraph networks output by \code{cmr_igraph}
#'@param fixed_locs (TRUE/FALSE). Whether to keep the the locations of the nodes fixed between network windows or allow them to change
#'@param locs A matrix or dataframe giving x and y coordinates of all nodes in the full network. Can be used to provide user-defined layouts to the graph - helpful for plotting locations in geographic space
#'@param dynamic (TRUE/FALSE). If true this function will plot networks for the different network windows sequentially in the same plot window with a delay. If false then this function will plot the networks in a single, multipanel figure
#'@param rows Sets the number of rows in the multipanel figure is dynamic=FALSE
#'@param layer_colours Set the colour of the plane and nodes for each layer of the multiplex network. Needs to be a vector the same length as the number of layers in the network
#'@details This function plots multiplex movement networks produced by CMRnet. It is developmental and it has fairly basic functionality currently. Networks can be plotted sequentially (dynamic=TRUE) or in a single multi-panelled figure (dynamic=FALSE)
#'@return A series of network plots showing the multiplex movement network from each network window
#'@examples
#'\dontrun{
#'
#'# load in data
#'data(cmrData2)
#'
#'# set parameters
#'mindate<-"2010-01-01"
#'maxdate<-"2015-01-01"
#'intwindow<-60
#'netwindow<-12
#'overlap<-0
#'
# Create movement networks
#'multimovenetdat <-
#'  MultiMoveNetCreate(
#'    data = cmrData2,
#'    intwindow = intwindow,
#'    mindate = mindate,
#'    maxdate = maxdate,
#'    netwindow = netwindow,
#'    overlap = overlap,
#'    nextonly = TRUE
#'  )
#'
#'# convert to networks
#'networks<-cmr_igraph(movenetdat,type="movement")
#'cmrMultiplexPlot(nets=networks,
#'                 fixed_locs=TRUE,
#'                 dynamic=FALSE,
#'                 rows=2,
#'                 layer_colours=c("firebrick","dodgerblue"))
#'}
#'@export

cmrMultiplexPlot<-function(nets,fixed_locs=c(TRUE,FALSE),locs=NULL,dynamic=c(TRUE,FALSE),rows=3,layer_colours){

  L<-length(nets[[1]])
  l<-length(nets[[1]][[1]])


  if(fixed_locs==TRUE&length(locs)==0){
    los<-list()
    tn<-nets[[2]][[1]]
    for(i in 1:l){
      tn<-igraph::union(tn,nets[[2]][[i]])
    }
    lo<-igraph::layout_with_fr(tn)
    for(i in 1:l){
      los[[i]]<-lo
    }
  }

  if(length(locs)>0){
    los<-list()
    for(i in 1:l){
      los[[i]]<-locs
    }
  }

  heights<-8*seq(0,l-1,1)+1

    if(dynamic==FALSE){
      graphics::par(mar=c(1,1,1,1),mfrow=c(rows,ceiling(L/rows)))
    }
    if(dynamic==TRUE){
      graphics::par(mar=c(1,1,1,1),mfrow=c(1,1))
    }
    for(i in 1:L){
      for(j in 1:l){
        if(fixed_locs==FALSE&length(locs)==0){
          lo2<-igraph::layout_with_fr(nets[[1]][[i]])
        }
        lo2<-los[[j]][igraph::vertex_attr(nets[[2]][[j]])$name%in%igraph::vertex_attr(nets[[1]][[i]][[j]])$name,]
        if(j==1){
          cD<-cbind(lo2,rep(heights[j],nrow(lo2)))
          C<-scatterplot3d::scatterplot3d(cD,xlim=c(min(lo2[,1])-1,max(lo2[,1])+1),ylim=c(min(lo2[,2])-1,max(lo2[,2])+1),zlim=c(0,max(heights)+3),color=layer_colours[j],pch=16,box=F,grid=F,cex.symbols=1.5,angle=70,axis=F,scale.y=0.5)
          theta<-seq(0,2*pi,length=1000)
          alpha<-pi/10
          ell.top.x<- (max(lo2[,1])+1)*cos(theta)*cos(alpha)-(min(lo2[,1])-1)*sin(theta)*sin(alpha)
          ell.top.y<- (max(lo2[,2])+1)*cos(theta)*sin(alpha)+(min(lo2[,2])-1)*sin(theta)*cos(alpha)
          C$points3d(ell.top.x,ell.top.y,rep(heights[j],length(ell.top.x)),type="l",col=grDevices::adjustcolor(layer_colours[j],alpha=1))
          graphics::polygon(C$xyz.convert(ell.top.x,ell.top.y,rep(heights[j],length(ell.top.x)))$x,C$xyz.convert(ell.top.x,ell.top.y,rep(heights[j],length(ell.top.x)))$y,col=grDevices::adjustcolor(layer_colours[j],0.08),border=NA)
          tmat<-igraph::as_adjacency_matrix(nets[[1]][[i]][[j]],sparse=FALSE)
          for(a in 1:(ncol(tmat)-1)){
            for(b in (i+1):ncol(tmat)){
              if(tmat[a,b]>0){
                C$points3d(c(cD[a,1],cD[b,1]),c(cD[a,2],cD[b,2]),c(cD[a,3],cD[b,3]),type="l",col="dark grey",lwd=3*tmat/max(tmat),lty=1)
              }
            }
          }
          C$points3d(cD[,1],cD[,2],cD[,3],col=layer_colours[j],pch=16,cex=1.5)
        }
        if(j>1){
          cD<-cbind(lo2,rep(heights[j],nrow(lo2)))
          theta<-seq(0,2*pi,length=1000)
          alpha<-pi/10
          ell.top.x<- (max(lo2[,1])+1)*cos(theta)*cos(alpha)-(min(lo2[,1])-1)*sin(theta)*sin(alpha)
          ell.top.y<- (max(lo2[,2])+1)*cos(theta)*sin(alpha)+(min(lo2[,2])-1)*sin(theta)*cos(alpha)
          C$points3d(ell.top.x,ell.top.y,rep(heights[j],length(ell.top.x)),type="l",col=grDevices::adjustcolor(layer_colours[j],alpha=1))
          graphics::polygon(C$xyz.convert(ell.top.x,ell.top.y,rep(heights[j],length(ell.top.x)))$x,C$xyz.convert(ell.top.x,ell.top.y,rep(heights[j],length(ell.top.x)))$y,col=grDevices::adjustcolor(layer_colours[j],0.08),border=NA)
          tmat<-igraph::as_adjacency_matrix(nets[[1]][[i]][[j]],sparse=FALSE)
          for(a in 1:(ncol(tmat)-1)){
            for(b in (i+1):ncol(tmat)){
              if(tmat[a,b]>0){
                C$points3d(c(cD[a,1],cD[b,1]),c(cD[a,2],cD[b,2]),c(cD[a,3],cD[b,3]),type="l",col="dark grey",lwd=3*tmat/max(tmat),lty=1)
              }
            }
          }
          C$points3d(cD[,1],cD[,2],cD[,3],col=layer_colours[j],pch=16,cex=1.5)
          if(dynamic==TRUE){Sys.sleep(1)}
        }
      }
    }

}

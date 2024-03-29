#'cmr_igraph
#'
#'Converts CMRnet networks to lists of igraph networks
#'
#'@param cmrnet A CMRnet object from \code{DynamicNetCreate}, \code{MoveNetCreate} or \code{MultiMoveNetCreate} functions
#'@param type The type of CMRnet object. Three possible arguments: "social" (object from \code{DynamicNetCreate}), "movement" (object from \code{MoveNetCreate}) and "multiplex" (object from \code{MultiMoveNetCreate})
#'@details This function converts CMRnet objects into lists of igraph networks for onward analysis and plotting
#'@return igraph networks for onward analysis. For all types of network the function returns a list of two outputs. For social and movement networks the first element is itself a list containing the igraph network corresponding to each network window. The second element of the list is the overall network calculated for all network windows combined. For multiplex movement networks the first list is more complex with an additional layer. The list of networks for different network windows is nested so that for each network window the movement network for each layer is stored separately.

#'@examples
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
#'
#'# example with multiplex networks ####
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
#'# create network
#'multimovenetdat<-MultiMoveNetCreate(data=cmrData2,
#'intwindow=intwindow,
#'mindate=mindate,
#'maxdate=maxdate,
#'netwindow=netwindow,
#'overlap=overlap,
#'nextonly=TRUE)
#'
#'# convert to networks
#'networks<-cmr_igraph(multimovenetdat,type="multiplex")
#'}
#'@export

cmr_igraph<-function(cmrnet,type=c("social","movement","multiplex")){

  if(type=="social"|type=="movement"){

    Fmatrix<-apply(cmrnet[[2]],1:2,sum)
    rownames(Fmatrix)<-colnames(Fmatrix)<-cmrnet[[3]][,1]

    matrices<-list()
    for(i in 1:dim(cmrnet[[2]])[3]){
      matrices[[i]]<-cmrnet[[2]][cmrnet[[3]][,i+1]==1,cmrnet[[3]][,i+1]==1,i]
      if(length(dim(matrices[[i]]))>0&dim(matrices[[i]])[1]>0){
        rownames(matrices[[i]])<-colnames(matrices[[i]])<-cmrnet[[3]][,1][cmrnet[[3]][,i+1]==1]
      }
      if(length(dim(matrices[[i]]))==0|dim(matrices[[i]])[1]==0){
        warning(paste("No nodes in network window",i))
      }
    }

    if(type=="social"){
      networks<-lapply(matrices,igraph::graph.adjacency,mode="undirected",weighted=TRUE)
      networkF<-igraph::graph.adjacency(Fmatrix,mode="undirected",weighted=TRUE)
    }
    if(type=="movement"){
      networks<-lapply(matrices,igraph::graph.adjacency,mode="directed",weighted=TRUE)
      networkF<-igraph::graph.adjacency(Fmatrix,mode="directed",weighted=TRUE)
    }

    output<-list(networks,networkF)
    return(output)

  }

  if(type=="multiplex"){

    nt<-length(cmrnet[[2]])
    ly<-dim(cmrnet[[3]])[3]
    Fmatrix<-array(0,dim=c(dim(cmrnet[[3]])[1],dim(cmrnet[[3]])[1],ly))

    for(i in 1:nt){
      for(j in 1:ly){
        Fmatrix[,,j]<-Fmatrix[,,j]+cmrnet[[2]][[i]][,,j]
        Fmatrix[,,j]<-Fmatrix[,,j]+cmrnet[[2]][[i]][,,j]
      }
    }
    rownames(Fmatrix)<-colnames(Fmatrix)<-rownames(cmrnet[[3]])

    networksF<-apply(Fmatrix,3,igraph::graph.adjacency,mode="directed",weighted=TRUE)

    networks<-list()
    for(i in 1:nt){
      networks[[i]]<-list()
      for(j in 1:ly){
        tmat1<-cmrnet[[2]][[i]][cmrnet[[3]][,i,j]==1,cmrnet[[3]][,i,j]==1,j]
        if(is.matrix(tmat1)){
          networks[[i]][[j]]<-igraph::graph.adjacency(tmat1,mode="directed",weighted=TRUE)
        }
        if(is.matrix(tmat1)==FALSE){
          warning("Less than two locations present in network layer")
        }
      }
    }

    output<-list(networks,networksF)
    return(output)

  }

}

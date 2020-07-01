
#'cmrNodeswap
#'
#'Conducts node-swap randomisations for CMRnet outputs
#'
#'@param cmrnet A CMRnet object from \code{DynamicNetCreate}, \code{MoveNetCreate} or \code{MultiMoveNetCreate}
#'@param n.rand The number of permuted matrices to be conducted by the rmperm() algorithm
#'@param multi (TRUE/FALSE). Indicates whether the input object is a monolayer network (from \code{DynamicNetCreate} or \code{MoveNetCreate}) or multiplex network (\code{MultiMoveNetCreate})
#'@details This function conducts network permutations of CMRnet objects, acting as a wrapper for \code{sna::rmperm()}
#'@return The randomised networks. If multi=FALSE then this consists of a list in which each element of the list corresponds to a particular network window and contains an array consisting of all of the randomised networks. If multi=TRUE then this consists of a nested list arrangement in which each element of the first list corresponds to a particular network window, and each element of the second list a layer of the multiplex for that network window. Each element of the second level contains an array consisting of all of the randommised versions of that layer in that network window

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
#'# run permutations
#'A<-cmrNodeswap(movenetdat,n.rand=1000)
#'
#'# example with multiple layered networks ####
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
#'# run permutations
#'B<-cmrNodeswap(multimovenetdat,n.rand=1000,multi=TRUE)
#'}
#'@export

cmrNodeswap<-function(cmrnet,n.rand,multi=FALSE){



  if(multi==FALSE){
    nets<-cmrnet[[2]]
    ne<-cmrnet[[3]]
    rnets<-list()

    # set up progress bar
    pb <- progress::progress_bar$new(total = dim(nets)[3]*n.rand, clear = FALSE)
    pb$tick(0)

    for(r in 1:dim(nets)[3]){
      rnets[[r]]<-array(NA,dim=c(dim(nets)[1],dim(nets)[2],n.rand))
      for( i in 1:n.rand){

        pb$tick()

        tmp_rnet<-nets[which(ne[,r+1]==1),which(ne[,r+1]==1),r]
        tmp_rnet<-sna::rmperm(tmp_rnet)
        tmp_rnet2<-array(0,dim(rnets[[r]])[1:2])
        tmp_rnet2[which(ne[,r+1]==1),which(ne[,r+1]==1)]<-tmp_rnet
        rnets[[r]][,,i]<-tmp_rnet2
      }
      rownames(rnets[[r]])<-colnames(rnets[[r]])<-rownames(nets[,,r])
    }
    return(rnets)
  }
  if(multi==TRUE){
    nets<-cmrnet[[2]]
    ne<-cmrnet[[3]]
    rnets<-list()

    # calculate total number of s
    s2 = 0
    for(r in 1:length(nets)){
      s2 <-  s2 + dim(nets[[r]])[3]
    }

    # set up progress bar
    pb <- progress::progress_bar$new(total = s2*n.rand, clear = FALSE)
    pb$tick(0)

    for(r in 1:length(nets)){
      rnets[[r]]<-list()
      for(s in 1:dim(nets[[r]])[3]){
        rnets[[r]][[s]]<-array(NA,dim=c(dim(nets[[r]])[1],dim(nets[[r]])[2],n.rand))
        for( i in 1:n.rand){

          pb$tick()

          tmp_rnet<-nets[[r]][which(ne[,r,s]==1),which(ne[,r,s]==1),s]
          tmp_rnet<-sna::rmperm(tmp_rnet)
          tmp_rnet2<-array(0,dim(rnets[[r]][[s]])[1:2])
          tmp_rnet2[which(ne[,r,s]==1),which(ne[,r,s]==1)]<-tmp_rnet
          rnets[[r]][[s]][,,i]<-tmp_rnet2
        }
        rownames(rnets[[r]][[s]])<-colnames(rnets[[r]][[s]])<-rownames(nets[[r]][,,s])
      }
    }
    return(rnets)
  }
}

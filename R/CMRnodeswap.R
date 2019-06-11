
#'cmrNodeswap
#'This function conducts network permutations of CMRnet objects, acting as a wrapper for the rmperm() function in the R package sna
#'Conducts node-swap randomisations for CMRnet outputs
#'
#'@param cmrnet A cmrnet object from the DynamicNetCreate, MoveNetCreate or MultiMoveNetCreate functions
#'@param n.rand The number of permuted matrices to be conducted by the rmperm() algorithm
#'@param multi (TRUE/FALSE). Indicates whether the input object is a monolayer network (from DynamicNetCreate or MoveNetCreate) or multiplex network (MultiMoveNetCreate)

#'@output The randomised networks. If multi=FALSE then this consists of a list in which each element of the list corresponds to a particular network window and contains an array consisting of all of the randomised networks. If multi=TRUE then this consists of a nested list arrangement in which each element of the first list corresponds to a particular network window, and each element of the second list a layer of the multiplex for that network window . Each element of the second level contains an array consisting of all of the randommised versions of that layer in that network window

#'@examples
#'data(cmr_dat)
#'mindate<-"2010-01-01"
#'maxdate<-"2015-01-01"
#'intwindow<-60
#'netwindow<-12
#'overlap<-0
#'movenetdat<-MoveNetCreate(data=cmr_dat,intwindow=intwindow,mindate=mindate,maxdate=maxdate,netwindow=netwindow,overlap=overlap,nextonly=TRUE)
#'A<-cmrNodeswap(movenetdat,n.rand=1000)
#'
#'data(cmr_dat2)
#'mindate<-"2010-01-01"
#'maxdate<-"2015-01-01"
#'intwindow<-60
#'netwindow<-12
#'overlap<-0
#'multimovenetdat<-MultiMoveNetCreate(data=cmr_dat,intwindow=intwindow,mindate=mindate,maxdate=maxdate,netwindow=netwindow,overlap=overlap,nextonly=TRUE)
#'B<-cmrNodeswap(multimovenetdat,n.rand=1000,multi=TRUE)
#'
#'@export

cmrNodeswap<-function(cmrnet,n.rand,multi=FALSE){

  require(sna)

  if(multi==FALSE){
    nets<-cmrnet[[2]]
    rnets<-list()
    for( r in 1:dim(nets)[3]){
      rnets[[r]]<-array(NA,dim=c(dim(nets)[1],dim(nets)[2],n.rand))
      for( i in 1:n.rand){
        rnets[[r]][,,i]<-sna::rmperm(nets[,,r])
      }
      rownames(rnets[[r]])<-colnames(rnets[[r]])<-rownames(nets[,,r])
    }
    return(rnets)
  }
  if(multi==TRUE){
    nets<-cmrnet[[2]]
    rnets<-list()
    for( r in 1:length(nets)){
      rnets[[r]]<-list()
      for(s in 1:dim(nets[[r]])[3]){
        rnets[[r]][[s]]<-array(NA,dim=c(dim(nets[[r]])[1],dim(nets[[r]])[2],n.rand))
        for( i in 1:n.rand){
          rnets[[r]][[s]][,,i]<-sna::rmperm(nets[[r]][,,s])
        }
        rownames(rnets[[r]][[s]])<-colnames(rnets[[r]][[s]])<-rownames(nets[[r]][,,s])
      }
    }
    return(rnets)
  }
}

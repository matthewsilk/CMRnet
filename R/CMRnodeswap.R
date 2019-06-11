###Node swap permutations for CMRnet objects
cmr.nodeswap<-function(cmrnet,n.rand,multi=FALSE){
  require(sna)
  if(multi==FALSE){
    nets<-cmrnet[[2]]
    rnets<-list()
    for( r in 1:dim(nets)[3]){
      rnets[[r]]<-array(NA,dim=c(dim(nets)[1],dim(nets)[2],n.rand))
      for( i in 1:n.rand){
        rnets[[r]][,,i]<-rmperm(nets[,,r])
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
          rnets[[r]][[s]][,,i]<-rmperm(nets[[r]][,,s])
        }
        rownames(rnets[[r]][[s]])<-colnames(rnets[[r]][[s]])<-rownames(nets[[r]][,,s])
      }
    }
    return(rnets)
  }
}
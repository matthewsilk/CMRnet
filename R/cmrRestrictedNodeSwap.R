
#'cmrRestrictedNodeswap
#'
#'Conducts node-swap randomisations for CMRnet outputs with restrictions on which nodes can be swapped
#'
#'@param cmrnet A cmrnet object from the DynamicNetCreate, MoveNetCreate or MultiMoveNetCreate functions
#'@param restrict An array showing which swaps are allowed and which are not within each network window (needs to be the same dimensions as the array being permuted). Value of 1 allows a swap. Value of 0 prevents a swap
#'@param n.rand The number of permutations
#'@param n.burnin The number of permutations to be discarded prior to recording results
#'@param n.swaps The number of swaps conducted before each output matrix is saved
#'@param multi (TRUE/FALSE). Indicates whether the input object is a monolayer network (from DynamicNetCreate or MoveNetCreate) or multiplex network (MultiMoveNetCreate)
#'@details This function conducts network permutations of CMRnet objects, acting as a wrapper for \code{sna::rmperm()}

#'@return The randomised networks. If multi=FALSE then this consists of a list in which each element of the list corresponds to a particular network window and contains an array consisting of all of the randomised networks. If multi=TRUE then this consists of a nested list arrangement in which each element of the first list corresponds to a particular network window, and each element of the second list a layer of the multiplex for that network window. Each element of the second level contains an array consisting of all of the randommised versions of that layer in that network window

#'@examples
#'
# example without multiple layers ####
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
#'res_arr<-array(sample(c(0,1),162000,replace=TRUE,prob=c(0.9,0.1)),dim=c(180,180,5))
#'diag(res_arr[,,1])<-diag(res_arr[,,2])<-0
#'diag(res_arr[,,3])<-diag(res_arr[,,4])<-diag(res_arr[,,5])<-0
#'
#'# run permutations
#'A<-cmrRestrictedNodeswap(cmrnet=movenetdat,
#'                                restrict=res_arr,
#'                                n.rand=1000,
#'                                n.burnin=500,
#'                                n.swaps=10,
#'                                multi=FALSE)
#'
#'@export

cmrRestrictedNodeswap<-function(cmrnet,restrict,n.rand,n.burnin,n.swaps,multi=FALSE){

  if(multi==FALSE){
    nets<-cmrnet[[2]]
    ne<-cmrnet[[3]]
    rnets<-list()

    tot.swaps<-n.burnin+n.rand*n.swaps
    tot.chain<-n.rand*n.swaps
    saves<-seq(n.burnin+n.swaps,tot.swaps,n.swaps)

    # set up progress bar
    pb <- progress::progress_bar$new(total = dim(nets)[3]*tot.swaps, clear = FALSE)
    pb$tick(0)

    for(r in 1:dim(nets)[3]){
      rnets[[r]]<-array(NA,dim=c(dim(nets)[1],dim(nets)[2],n.rand))

      tmp_rnet<-nets[which(ne[,r+1]==1),which(ne[,r+1]==1),r]
      tmp_restrict<-restrict[which(ne[,r+1]==1),which(ne[,r+1]==1),r]

      c<-1
      for( i in 1:tot.swaps){

        pb$tick()

        smp1<-sample(1:nrow(tmp_rnet),1)
        smp2<-sample(seq(1,nrow(tmp_rnet),1)[-smp1],1)

        if(tmp_restrict[smp1,smp2]==0){
          tmp_rnet2<-array(0,dim(rnets[[r]])[1:2])
          tmp_rnet2[which(ne[,r+1]==1),which(ne[,r+1]==1)]<-tmp_rnet
          if(i%in%saves){
            rnets[[r]][,,c]<-tmp_rnet2
            c<-c+1
          }
        }
        if(tmp_restrict[smp1,smp2]==1){
          tmp_rnetT<-tmp_rnet
          tmp_rnet[smp1,]<-tmp_rnetT[smp2,]
          tmp_rnet[smp2,]<-tmp_rnetT[smp1,]
          tmp_rnet[,smp1]<-tmp_rnetT[,smp2]
          tmp_rnet[,smp2]<-tmp_rnetT[,smp1]
          tmp_rnet2<-array(0,dim(rnets[[r]])[1:2])
          tmp_rnet2[which(ne[,r+1]==1),which(ne[,r+1]==1)]<-tmp_rnet
          if(i%in%saves){
            rnets[[r]][,,c]<-tmp_rnet2
            c<-c+1
          }
        }

      }
      rownames(rnets[[r]])<-colnames(rnets[[r]])<-rownames(nets[,,r])
    }
    return(rnets)
  }

  if(multi==TRUE){
    nets<-cmrnet[[2]]
    ne<-cmrnet[[3]]
    rnets<-list()

    tot.swaps<-n.burnin+n.rand*n.swaps
    tot.chain<-n.rand*n.swaps
    saves<-seq(n.burnin+n.swaps,tot.swaps,n.swaps)

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

        tmp_rnet<-nets[[r]][which(ne[,r,s]==1),which(ne[,r,s]==1),s]
        tmp_restrict<-restrict[[r]][which(ne[,r,s]==1),which(ne[,r,s]==1),s]

          c<-1
          for( i in 1:tot.swaps){

            pb$tick()

            smp1<-sample(1:nrow(tmp_rnet),1)
            smp2<-sample(seq(1,nrow(tmp_rnet),1)[-smp1],1)

            if(tmp_restrict[smp1,smp2]==0){
              tmp_rnet2<-array(0,dim(rnets[[r]])[1:2])
              tmp_rnet2[which(ne[,r,s]==1),which(ne[,r,s]==1)]<-tmp_rnet
              if(i%in%saves){
                rnets[[r]][[s]][,,c]<-tmp_rnet2
                c<-c+1
              }
            }
            if(tmp_restrict[smp1,smp2]==1){
              tmp_rnetT<-tmp_rnet
              tmp_rnet[smp1,]<-tmp_rnetT[smp2,]
              tmp_rnet[smp2,]<-tmp_rnetT[smp1,]
              tmp_rnet[,smp1]<-tmp_rnetT[,smp2]
              tmp_rnet[,smp2]<-tmp_rnetT[,smp1]
              tmp_rnet2<-array(0,dim(rnets[[r]])[1:2])
              tmp_rnet2[which(ne[,r,s]==1),which(ne[,r,s]==1)]<-tmp_rnet
              if(i%in%saves){
                rnets[[r]][[s]][,,c]<-tmp_rnet2
                c<-c+1
              }
            }


        }
        rownames(rnets[[r]][[s]])<-colnames(rnets[[r]][[s]])<-rownames(nets[[r]][,,s])
      }
    }
    return(rnets)
  }
}

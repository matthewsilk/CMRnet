#'cmrperm.soc
#'An internal function that operate within the DatastreamPermSpat function
#'
#'@param D The input dataset to be randomised
#'@param same.time (TRUE/FALSE) Whether swaps should be restricted to only occur betwen individuals trapped on the same date or not
#'@param time.restrict Provided as a number of months. Imposes time restrictions on when swaps can take place so that individuals can only be swapped with those a fixed time before or after being captured
#'@param same.spat (TRUE/FALSE) Whether swaps should be restricted to only occur between individuals trapped at the same location
#'@param spat.restrict Provided on the same scale as the coordinates in the input dataset. Imposes space restrictions on when swaps can take place so that individuals can only be swapped with those captued within a fixed distance
#'@param n.swaps The number of swaps between each random network being extracted (e.g. n.swaps = 10 would equate to 10 swaps taking place between each random network being saved)
#'@param n.rand The number of randomised networks to be generated
#'@param n.burnin The number of swaps to discard as burn-in before the first random network is created. The total number of swaps conducted is thus n.burnin+n.swaps*n.rand
#'@param warn.thresh The number of times no matches are found (i.e. constraints on randomisations are too restrictive) before the function is stopped and an error message returned
#'
#'@output A randomised dataset with the same dimensions as the original input dataset
#'
#'@examples
#'data(cmr_dat)
#'mindate<-"2010-01-01"
#'maxdate<-"2015-01-01"
#'intwindow<-60
#'netwindow<-12
#'overlap<-0
#'movenetdat<-MoveNetCreate(data=cmr_dat,intwindow=intwindow,mindate=mindate,maxdate=maxdate,netwindow=netwindow,overlap=overlap,nextonly=TRUE)
#'
#'same.time=FALSE
#'time.restrict=6
#'same.id=FALSE
#'n.swaps=10
#'n.rand=100
#'n.burnin=100
#'warn.thresh=100
#'
#'Rs<-DatastreamPermSpat(data=cmr_dat,intwindow,mindate,maxdate,netwindow,overlap,spacewindow,same.time,time.restrict,same.id,n.swaps,n.rand,n.burnin,warn.thresh)

#'@export

cmrperm.spat<-function(D,same.time,time.restrict,same.id,n.swaps,n.rand,n.burnin,warn.thresh){

  D.rand<-list()

  warns<-0

  ctr<-1

  for(sw in 1:(n.rand*n.swaps+n.burnin)){
    tmp1<-sample(1:nrow(D),1)
    if(same.time==TRUE){
      tmpdays<-D$Jdays[tmp1]
    } else if(time.restrict!="n"){
      tmpmaxT<-AddMonths(D$date[tmp1],time.restrict)
      tmpminT<-AddMonths(D$date[tmp1],-time.restrict)
      tmpdays<-julian(seq(as.Date(tmpminT),as.Date(tmpmaxT),by="day"),origin=as.Date("1970-01-01"))
    } else{
      tmpdays<-sort(unique(D$Jdays))
    }
    if(same.id==TRUE){
      tmpids<-D$id[tmp1]
    } else{
      tmpids<-sort(unique(D$id))
    }

    poss<-which(D$id%in%tmpids&D$Jdays%in%tmpdays)
    poss<-poss[which(poss%in%tmp1==FALSE)]

    if(length(poss)>0){

      tmp2<-sample(poss,1)

      tmp.loc1<-D$loc[tmp1]
      tmp.loc2<-D$loc[tmp2]

      D$loc[tmp1]<-tmp.loc2
      D$loc[tmp2]<-tmp.loc1

      if(ctr>n.burnin){
        if((ctr-n.burnin)%%n.swaps==0){
          D.rand[[(ctr-n.burnin)/n.swaps]]<-D
        }
      }

      ctr<-ctr+1

    } #end if loop



    if(length(poss)==0){

      warning('No matches found. May need to reconsider permutation restrictions')
      warns<-warns+1

    }

    if (warns>warn.thresh){
      stop('Permutation restrictions need changing')
    }


  }

  return(D.rand)

}

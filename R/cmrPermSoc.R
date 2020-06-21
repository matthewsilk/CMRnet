#'cmrPermSoc
#'
#'An internal function that operate within CMRnet::DatastreamPermSoc
#'
#'@param D The input dataset to be randomised
#'@param locmat A matrix indicating whether distance between the capture locations is less than spat.restrict if not "n"
#'@param same.time (TRUE/FALSE) Whether swaps should be restricted to only occur betwen individuals trapped on the same date or not
#'@param time.restrict Provided as a number of months. Imposes time restrictions on when swaps can take place so that individuals can only be swapped with those a fixed time before or after being captured
#'@param same.spat (TRUE/FALSE) Whether swaps should be restricted to only occur between individuals trapped at the same location
#'@param spat.restrict Provided on the same scale as the coordinates in the input dataset. Imposes space restrictions on when swaps can take place so that individuals can only be swapped with those captued within a fixed distance
#'@param n.swaps The number of swaps between each random network being extracted (e.g. n.swaps = 10 would equate to 10 swaps taking place between each random network being saved)
#'@param n.rand The number of randomised networks to be generated
#'@param burnin (TRUE/FALSE) Whether burnin is required
#'@param n.burnin The number of swaps to discard as burn-in before the first random network is created. The total number of swaps conducted is thus n.burnin+n.swaps*n.rand
#'@param buffer The number of days around a capture event that an individual can't be swapped into. Defaults to zero. This exists for when trapping schedules would prevent an individual trapped on one day being trapped on adjacent days
#'

#'@return A randomised dataset with the same dimensions as the original input dataset
#'@export

cmrPermSoc<-function(D,locmat,same.time,time.restrict,same.spat,spat.restrict,n.swaps,n.rand,burnin,n.burnin,buffer=0){

  D.rand<-list()

  ctr<-1

  if(burnin==FALSE){n.burnin<-0}

  for(sw in 1:(n.rand*n.swaps+n.burnin)){
    tmp1<-sample(1:nrow(D),1)
    if(same.time==TRUE){
      tmpdays<-D$Jdays[tmp1]
    } else if(time.restrict!="n"){
      tmpmaxT<-DescTools::AddMonths(D$date[tmp1],time.restrict)
      tmpminT<-DescTools::AddMonths(D$date[tmp1],-time.restrict)
      tmpdays<-julian(seq(as.Date(tmpminT),as.Date(tmpmaxT),by="day"),origin=as.Date("1970-01-01"))
    } else{
      tmpdays<-sort(unique(D$Jdays))
    }
    if(same.spat==TRUE){
      tmpsites<-D$loc[tmp1]
    } else if(spat.restrict!="n"){
      tmpsites<-rownames(locmat)[which(locmat[,which(colnames(locmat)==D$loc[tmp1])]==TRUE)]
    } else{
      tmpsites<-sort(unique(D$loc))
    }

    poss<-which(D$loc%in%tmpsites&D$Jdays%in%tmpdays)
    poss<-poss[which(poss%in%tmp1==FALSE)]

    tmp2<-sample(poss,1)

    tmp.id1<-D$id[tmp1]
    tmp.id2<-D$id[tmp2]

    #check that tmp.id1 is not trapped on the same day that it will be swapped to already (or within a buffer number of days)
    #check that tmp.id2 is not trapped on the same day that it will be swapped to already (or within a buffer number of days)

    check<-0

    tmp.date1<-D$Jdays[tmp1]
    tmp.date2<-D$Jdays[tmp2]
    tmp.days1<-D$Jdays[D$id==tmp.id1]
    tmp.days2<-D$Jdays[D$id==tmp.id2]

    if(buffer>0){
      extra.days<-seq(-buffer,buffer,1)
      tmp.days1b<-tmp.days1-buffer
      tmp.days2b<-tmp.days2-buffer
      for(ii in 2:length(extra.days)){
        tmp.days1b<-c(tmp.days1b,tmp.days1+extra.days[ii])
        tmp.days2b<-c(tmp.days2b,tmp.days2+extra.days[ii])
      }
      tmp.days1<-sort(unique(tmp.days1b))
      tmp.days2<-sort(unique(tmp.days2b))
    }

    if((sum(tmp.days1%in%tmp.date2)>0|sum(tmp.days2%in%tmp.date1)>0)==FALSE){check<-1}

    if(check==1){

      D$id[tmp1]<-tmp.id2
      D$id[tmp2]<-tmp.id1

    }

    if(ctr>n.burnin){
      if((ctr-n.burnin)%%n.swaps==0){
        D.rand[[(ctr-n.burnin)/n.swaps]]<-D
      }
    }

    ctr<-ctr+1

  }

  return(D.rand)

}

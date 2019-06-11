
##Conducts datastream permutations for movement network analysis of CMR data

##Takes original dataset and swaps sites within particular time restrictions (set in months) and potentially constrained to only swap within an individual

##n.rand is the number of networks to return
##n.swaps is the number of swaps between each randomisation


###swap function internal

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

###main function

datastream.perm.spat<-function(data,intwindow,mindate,maxdate,netwindow,overlap,next.only=FALSE,same.time,time.restrict,same.id,n.swaps,n.rand,n.burnin,warn.thresh){
  
  require(chron)
  require(DescTools)
  
  D<-data
  names(D)<-c("id","loc","x","y","date")
  Jdays<-julian(as.Date(D$date),origin=as.Date("1970-01-01"))
  D<-data.frame(D,Jdays)
  D<-D[order(D$Jdays, D$loc,D$id),]
  D$id<-as.factor(D$id)
  D$loc<-as.factor(D$loc)
  D$x<-as.numeric(D$x)
  D$y<-as.numeric(D$y)
  
  L<-netwindow
  O<-overlap
  start<-julian(as.Date(mindate),origin=as.Date("1970-01-01"))
  end<-julian(as.Date(maxdate),origin=as.Date("1970-01-01"))
  days<-seq(start,end,1)
  length<-length(days)
  
  month_seq<-seq(as.Date(mindate),as.Date(maxdate), by = "month")
  
  starts<-month_seq[seq(1,length(month_seq)-L,L-O)]
  ends<-month_seq[which(month_seq%in%starts)+(L)]
  
  starts<-julian(as.Date(starts),origin=as.Date("1970-01-01"))
  ends<-julian(as.Date(ends),origin=as.Date("1970-01-01"))
  
  #Provide warning message if the netwindows stop early
  print(paste0("stopped",end-ends[length(ends)],"days early"))
  
  #Counts the number of windows over which networks are built
  Ws<-length(starts)
  
  #size of sliding window
  X<-intwindow
  
  #get only the data of interest
  #Having less than end 
  D2<-D[which(D$Jdays>=start&D$Jdays<end),]
  
  #extract unique individuals and record how many there are
  ids<-sort(unique(D2$id))
  n.ids<-length(ids)
  
  #extract unique locations and record how many there are
  locs<-sort(unique(D2$loc))
  n.locs<-length(locs)
  
  n.caps<-length(D2[,1])
  
  EDGES<-list()
  NET<-list()
  
  E1<-rep(locs,each=n.locs-1)
  
  E1<-factor(E1,levels=levels(locs))
  
  E2<-locs[2:length(locs)]
  for(i in 2:n.locs){
    E2<-c(as.character(E2),as.character(locs[-i]))
  }
  
  E2<-factor(E2,levels=levels(locs))
  
  for(ts in 1:Ws){
    EDGES[[ts]]<-array(0,dim=c(((n.locs-1)*(n.locs)),3,n.rand))
    EDGES[[ts]][,1,]<-E1
    EDGES[[ts]][,2,]<-E2
    
    NET[[ts]]<-array(0,dim=c(n.ids,n.ids,n.rand))
    colnames(NET[[ts]])<-ids
    rownames(NET[[ts]])<-ids
  }
  
  NODE.EXIST<-matrix(0,nr=n.ids,nc=Ws)
  
  #Less than ends
  for (ts in 1:Ws){
    
    D3<-D2[which(D2$Jdays>=starts[ts]&D2$Jdays<ends[ts]),]
    D3$id<-factor(D3$id,levels=levels(D2$id))
    D3$loc<-factor(D3$loc,levels=levels(D2$loc))
    n.Caps2<-length(D3$id)
    
    #loop through IDs and record whether each one was recorded
    #in this time period with a binary response
    for (i in 1:n.locs){
      
      ifelse(locs[i]%in%D3$loc>0,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts]+1,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts])
      #print(paste(ts,"-",i,"-tickB"))
    }
    
    rands<-cmrperm.spat(D=D3,same.time=same.time,time.restrict=time.restrict,same.id=same.id,n.swaps=n.swaps,n.rand=n.rand,n.burnin=n.burnin)
    
    for(r in 1:length(rands)){
      
      for (i in 1:n.Caps2){
        D3<-rands[[r]]
        range<-seq(D3$Jdays[i],D3$Jdays[i]+X)
        timematch<-which(D3$Jdays%in%range==TRUE)
        idmatch<-which(D3$id%in%D3$id[i]==TRUE)
        twomatch<-timematch[which(timematch%in%idmatch==TRUE)]
        MATCH<-twomatch[-which(twomatch==i)]
        
        if(next.only==TRUE){MATCH<-MATCH[which.min(D3$Jdays[MATCH]-D3$Jdays[i])]}
        
        
        for (j in 1:length(MATCH)){
          EDGES[[ts]][which(EDGES[[ts]][,1,r]%in%as.numeric(D3$loc[i])==TRUE&EDGES[[ts]][,2,r]%in%as.numeric(D3$loc[MATCH[j]])==TRUE),3,r]<-EDGES[[ts]][which(EDGES[[ts]][,1,r]%in%as.numeric(D3$loc[i])==TRUE&EDGES[[ts]][,2,r]%in%as.numeric(D3$loc[MATCH[j]])==TRUE),3,r]+1
        }
        
      } #end loop over n.caps
      
      ##and now turn the edge list into an association matrix as well to put network in double format
      NET.rows<-as.numeric(factor(rownames(NET[[ts]]),levels=levels(D$loc)))
      
      for (i in 1:length(EDGES[[ts]][,3,r])){
        NET[[ts]][which(NET.rows%in%EDGES[[ts]][i,1,r]==TRUE),which(NET.rows%in%EDGES[[ts]][i,2,r]==TRUE),r]<-NET[[ts]][which(NET.rows%in%EDGES[[ts]][i,1,r]==TRUE),which(NET.rows%in%EDGES[[ts]][i,2,r]==TRUE),r]+EDGES[[ts]][i,3,r]
      }
      
      print(paste(ts,"-",r))
    } #end r loop over randomisations
    
    print(paste(ts,"done"))
    #end loop over ts/Ws
  }
  
  NODE.EXIST<-data.frame(ids,NODE.EXIST)
  
  results<-list(EDGES,NET,NODE.EXIST)
  
  return(results)
  
}
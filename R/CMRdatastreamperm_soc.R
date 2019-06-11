
##Conducts datastream permutations for social network analysis of CMR data

##Takes original dataset and swaps individuals within particular time restrictions (set in months) and space restrictions (only within the same site if same.pat=TRUE or based on the coordinate system provided if spat.restrict is given)

##n.rand is the number of networks to return
##n.swaps is the number of swaps between each randomisation


###swap function internal

cmrperm.soc<-function(D,locmat,same.time,time.restrict,same.spat,spat.restrict,n.swaps,n.rand,n.burnin){

  D.rand<-list()
  
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
    if(same.spat==TRUE){
      tmpsites<-D$loc[tmp1]
    } else if(spat.restrict!="n"){
      tmpsites<-rownames(locmat2)[which(locmat[,which(colnames(locmat)==D$loc[tmp1])]==TRUE)]
    } else{
      tmpsites<-sort(unique(D$loc))
    }
    
    poss<-which(D$loc%in%tmpsites&D$Jdays%in%tmpdays)
    poss<-poss[which(poss%in%tmp1==FALSE)]
    
    tmp2<-sample(poss,1)
    
    tmp.id1<-D$id[tmp1]
    tmp.id2<-D$id[tmp2]
    
    D$id[tmp1]<-tmp.id2
    D$id[tmp2]<-tmp.id1
    
    if(ctr>n.burnin){
      if((ctr-n.burnin)%%n.swaps==0){
        D.rand[[(ctr-n.burnin)/n.swaps]]<-D
      }
    }
    
    ctr<-ctr+1
    
  }
  
  return(D.rand)
  
}

###Edge list and association matrix creation function internal

#net.create<-function(fullD,subD,locmat,t_EDGES,t_NET){

#  for (i in 1:n.Caps2){
#    range<-seq(D3$Jdays[i]-X,D3$Jdays[i]+X)
#    timematch<-which(D3$Jdays%in%range==TRUE)
#    spacematch<-which(D3$loc%in%D3$loc[i]==TRUE)
#    spacematch<-which(D3$loc%in%rownames(locmat2)[which(locmat2[,which(colnames(locmat2)==D3$loc[i])]==TRUE)]==TRUE)
#    twomatch<-timematch[which(timematch%in%spacematch==TRUE)]
#    MATCH<-twomatch[-which(twomatch==i)]
  
#    for (j in 1:length(MATCH)){
#      EDGES[which(EDGES[,1,ts]%in%as.numeric(D3$id[i])==TRUE&EDGES[,2,ts]%in%as.numeric(D3$id[MATCH[j]])==TRUE),3,ts]<-EDGES[which(EDGES[,1,ts]%in%as.numeric(D3$id[i])==TRUE&EDGES[,2,ts]%in%as.numeric(D3$id[MATCH[j]])==TRUE),3,ts]+1
#    }
#  
#  }

  ##and now turn the edge list into an association matrix as well to put network in double format
#  NET.rows<-as.numeric(factor(rownames(NET),levels=levels(D$id)))
#
#  for (i in 1:length(EDGES[,3,ts])){
#   NET[which(NET.rows%in%EDGES[i,1,ts]==TRUE),which(NET.rows%in%EDGES[i,2,ts]==TRUE),ts]<-NET[which(NET.rows%in%EDGES[i,1,ts]==TRUE),which(NET.rows%in%EDGES[i,2,ts]==TRUE),ts]+EDGES[i,3,ts]
#    NET[which(NET.rows%in%EDGES[i,2,ts]==TRUE),which(NET.rows%in%EDGES[i,1,ts]==TRUE),ts]<-NET[which(NET.rows%in%EDGES[i,2,ts]==TRUE),which(NET.rows%in%EDGES[i,1,ts]==TRUE),ts]+EDGES[i,3,ts]
#  }

#}
  
###main function

datastream.perm.soc<-function(data,intwindow,mindate,maxdate,netwindow,overlap,spacewindow,same.time,time.restrict,same.spat,spat.restrict,n.swaps,n.rand,n.burnin){
  
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
  
  locdat<-aggregate(D[,2:4],by=list(D$loc),unique)[,2:4]
  locmat<-as.matrix(dist(locdat[,2:3]))
  locmat2<-locmat<spat.restrict
  rownames(locmat2)<-colnames(locmat2)<-locs
  
  n.caps<-length(D2[,1])
  
  EDGES<-list()
  NET<-list()
  
  e1<-seq((n.ids-1),0,-1)
  E1<-rep(ids[1],e1[1])
  for(i in 2:(n.ids)){
    E1<-c(as.character(E1),rep(as.character(ids[i]),e1[i]))
  }
  
  E1<-factor(E1,levels=levels(ids))
  
  E2<-ids[2:length(ids)]
  for(i in 3:n.ids){
    E2<-c(as.character(E2),as.character(ids[i:length(ids)]))
  }
  
  E2<-factor(E2,levels=levels(ids))
  
  for(ts in 1:Ws){
    EDGES[[ts]]<-array(0,dim=c(((n.ids-1)*(n.ids))/2,3,n.rand))
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
    n.Caps2<-length(D3$id)
    
    #loop through IDs and record whether each one was recorded
    #in this time period with a binary response
    for (i in 1:n.ids){
      
      ifelse(ids[i]%in%D3$id>0,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts]+1,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts])
      #print(paste(ts,"-",i,"-tickB"))
    }
    
    rands<-cmrperm.soc(D=D3,locmat=locmat2,same.time=same.time,time.restrict=time.restrict,same.spat=same.spat,spat.restrict=spat.restrict,n.swaps=n.swaps,n.rand=n.rand,n.burnin=n.burnin)
      
    for(r in 1:length(rands)){
    
      for (i in 1:n.Caps2){
        D3<-rands[[r]]
        range<-seq(D3$Jdays[i]-X,D3$Jdays[i]+X)
        timematch<-which(D3$Jdays%in%range==TRUE)
        spacematch<-which(D3$loc%in%D3$loc[i]==TRUE)
        spacematch<-which(D3$loc%in%rownames(locmat2)[which(locmat2[,which(colnames(locmat2)==D3$loc[i])]==TRUE)]==TRUE)
        twomatch<-timematch[which(timematch%in%spacematch==TRUE)]
        MATCH<-twomatch[-which(twomatch==i)]
      
        for (j in 1:length(MATCH)){
          EDGES[[ts]][which(EDGES[[ts]][,1,r]%in%as.numeric(D3$id[i])==TRUE&EDGES[[ts]][,2,r]%in%as.numeric(D3$id[MATCH[j]])==TRUE),3,r]<-EDGES[[ts]][which(EDGES[[ts]][,1,r]%in%as.numeric(D3$id[i])==TRUE&EDGES[[ts]][,2,r]%in%as.numeric(D3$id[MATCH[j]])==TRUE),3,r]+1
        }
      
      }
    
      ##and now turn the edge list into an association matrix as well to put network in double format
      NET.rows<-as.numeric(factor(rownames(NET[[ts]]),levels=levels(D$id)))
    
      for (i in 1:length(EDGES[[ts]][,3,r])){
        NET[[ts]][which(NET.rows%in%EDGES[[ts]][i,1,r]==TRUE),which(NET.rows%in%EDGES[[ts]][i,2,r]==TRUE),r]<-NET[[ts]][which(NET.rows%in%EDGES[[ts]][i,1,r]==TRUE),which(NET.rows%in%EDGES[[ts]][i,2,r]==TRUE),r]+EDGES[[ts]][i,3,r]
        NET[[ts]][which(NET.rows%in%EDGES[[ts]][i,2,r]==TRUE),which(NET.rows%in%EDGES[[ts]][i,1,r]==TRUE),ts]<-NET[[ts]][which(NET.rows%in%EDGES[[ts]][i,2,r]==TRUE),which(NET.rows%in%EDGES[[ts]][i,1,r]==TRUE),ts]+EDGES[[ts]][i,3,r]
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
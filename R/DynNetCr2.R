
#takes data in the form (ID,Spatialfactor,X,Y,Date) and turns into
#edgelist and association matrix

##please give netwindow as a number of months (e.g. 6 = 6 months, 24 = 24 months etc)
##and overlap as number from 0 (no overlap) to the value of netwindow (i.e. in months))
##Maxdate is always the day after the end of the period of interest

dynamic.net.create<-function(data,intwindow,mindate,maxdate,netwindow,overlap,spacewindow){
  
  #Add a column with Julian Dates
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
  
  #calculate distance matrix between locations
  locdat<-aggregate(D[,2:4],by=list(D$loc),unique)[,2:4]
  locmat<-as.matrix(dist(locdat[,2:3]))
  locmat2<-locmat<=spacewindow
  rownames(locmat2)<-colnames(locmat2)<-locs
  
  n.caps<-length(D2[,1])
  
  EDGES<-array(0,dim=c(((n.ids-1)*(n.ids))/2,3,Ws))
  
  NET<-array(0,dim=c(n.ids,n.ids,Ws))
  colnames(NET)<-ids
  rownames(NET)<-ids
  
  e1<-seq((n.ids-1),0,-1)
  E1<-rep(ids[1],e1[1])
  for(i in 2:(n.ids)){
    E1<-c(as.character(E1),rep(as.character(ids[i]),e1[i]))
  }
  
  E1<-factor(E1,levels=levels(ids))
  EDGES[,1,]<-E1
  
  E2<-ids[2:length(ids)]
  for(i in 3:n.ids){
    E2<-c(as.character(E2),as.character(ids[i:length(ids)]))
  }
  
  E2<-factor(E2,levels=levels(ids))
  EDGES[,2,]<-E2

  NODE.EXIST<-matrix(0,nr=n.ids,nc=Ws)
  
  EDGE.EXIST<-matrix(0,nr=length(EDGES[,1,1]),nc=Ws)
  
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
    
    
    for (i in 1:n.Caps2){
      range<-seq(D3$Jdays[i]-X,D3$Jdays[i]+X)
      timematch<-which(D3$Jdays%in%range==TRUE)
      spacematch<-which(D3$loc%in%D3$loc[i]==TRUE)
      spacematch<-which(D3$loc%in%rownames(locmat2)[which(locmat2[,which(colnames(locmat2)==D3$loc[i])]==TRUE)]==TRUE)
      twomatch<-timematch[which(timematch%in%spacematch==TRUE)]
      MATCH<-twomatch[-which(twomatch==i)]
      
      for (j in 1:length(MATCH)){
        EDGES[which(EDGES[,1,ts]%in%as.numeric(D3$id[i])==TRUE&EDGES[,2,ts]%in%as.numeric(D3$id[MATCH[j]])==TRUE),3,ts]<-EDGES[which(EDGES[,1,ts]%in%as.numeric(D3$id[i])==TRUE&EDGES[,2,ts]%in%as.numeric(D3$id[MATCH[j]])==TRUE),3,ts]+1
      }
      
      #print(paste(ts,"-",i,"-tickC"))
      
    }
    
    ##and now turn the edge list into an association matrix as well to put network in double format
    NET.rows<-as.numeric(factor(rownames(NET),levels=levels(D$id)))
    
    for (i in 1:length(EDGES[,3,ts])){
      NET[which(NET.rows%in%EDGES[i,1,ts]==TRUE),which(NET.rows%in%EDGES[i,2,ts]==TRUE),ts]<-NET[which(NET.rows%in%EDGES[i,1,ts]==TRUE),which(NET.rows%in%EDGES[i,2,ts]==TRUE),ts]+EDGES[i,3,ts]
      NET[which(NET.rows%in%EDGES[i,2,ts]==TRUE),which(NET.rows%in%EDGES[i,1,ts]==TRUE),ts]<-NET[which(NET.rows%in%EDGES[i,2,ts]==TRUE),which(NET.rows%in%EDGES[i,1,ts]==TRUE),ts]+EDGES[i,3,ts]
    }
    
    print(paste(ts,"done"))
    
    #end loop over ts/Ws
  }
  
  NODE.EXIST<-data.frame(ids,NODE.EXIST)
  
  results<-list(EDGES,NET,NODE.EXIST,E1,E2)
  
  return(results)
  
  #end function
}
  
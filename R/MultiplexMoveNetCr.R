
#takes data in the form (ID,Spatialfactor,X,Y,Date,LayerFactor) and turns into
#edgelists and association matrices

##please give netwindow as a number of months (e.g. 6 = 6 months, 24 = 24 months etc)
##and overlap as number from 0 (no overlap) to the value of netwindow (i.e. in months))
##Maxdate is always the day after the end of the period of interest

multi.move.net.create<-function(data,intwindow,mindate,maxdate,netwindow,overlap,next.only=FALSE){
  
  #Add a column with Julian Dates
  D<-data
  names(D)<-c("id","loc","x","y","date","layer")
  Jdays<-julian(as.Date(D$date),origin=as.Date("1970-01-01"))
  D<-data.frame(D,Jdays)
  D<-D[order(D$Jdays, D$loc,D$id),]
  D$id<-as.factor(D$id)
  D$loc<-as.factor(D$loc)
  D$layer<-as.factor(D$layer)
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
  
  layers<-sort(unique(D2$layer))
  n.layers<-length(layers)
  
  n.caps<-length(D2[,1])
  
  EDGES<-list()
  for(ts in 1:Ws){
  EDGES[[ts]]<-array(0,dim=c(((n.locs-1)*(n.locs)),3,n.layers))
  }
  
  NET<-list()
  for(ts in 1:Ws){
  NET[[ts]]<-array(0,dim=c(n.locs,n.locs,n.layers))
  colnames(NET[[ts]])<-locs
  rownames(NET[[ts]])<-locs
  }
  
  E1<-rep(locs,each=n.locs-1)
  
  E1<-factor(E1,levels=levels(locs))
  for(ts in 1:Ws){
  EDGES[[ts]][,1,]<-E1
  }
  
  E2<-locs[2:length(locs)]
  for(i in 2:n.locs){
    E2<-c(as.character(E2),as.character(locs[-i]))
  }
  
  E2<-factor(E2,levels=levels(locs))
  for(ts in 1:Ws){
  EDGES[[ts]][,2,]<-E2
  }
  
  NODE.EXIST<-matrix(0,nr=n.locs,nc=Ws)
  
  #Less than ends
  for (ts in 1:Ws){
    
    D3<-D2[which(D2$Jdays>=starts[ts]&D2$Jdays<ends[ts]),]
    D3$id<-factor(D3$id,levels=levels(D2$id))
    D3$loc<-factor(D3$loc,levels=levels(D2$loc))
    D3$layer<-factor(D3$layer,levels=levels(D2$layer))
    n.Caps2<-length(D3$id)
    
    #loop through IDs and record whether each one was recorded
    #in this time period with a binary response
    for (i in 1:n.locs){
      
      ifelse(locs[i]%in%D3$loc>0,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts]+1,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts])
      #print(paste(ts,"-",i,"-tickB"))
    }
    
    for(ls in 1:n.layers){
      
      D4<-D3[which(D3$layer==layers[ls]),]
      D4$id<-factor(D4$id,levels=levels(D3$id))
      D4$loc<-factor(D4$loc,levels=levels(D3$loc))
      D4$layer<-factor(D4$layer,levels=levels(D3$layer))
      n.Caps2<-length(D4$id)
    
      for (i in 1:n.Caps2){
        range<-seq(D4$Jdays[i],D4$Jdays[i]+X)
        timematch<-which(D4$Jdays%in%range==TRUE)
        idmatch<-which(D4$id%in%D4$id[i]==TRUE)
        twomatch<-timematch[which(timematch%in%idmatch==TRUE)]
        MATCH<-twomatch[-which(twomatch==i)]

        if(next.only==TRUE){MATCH<-MATCH[which.min(D4$Jdays[MATCH]-D4$Jdays[i])]}
      
        for (j in 1:length(MATCH)){
          EDGES[[ts]][which(EDGES[[ts]][,1,ls]%in%as.numeric(D4$loc[i])==TRUE&EDGES[[ts]][,2,ls]%in%as.numeric(D4$loc[MATCH[j]])==TRUE),3,ls]<-EDGES[[ts]][which(EDGES[[ts]][,1,ls]%in%as.numeric(D4$loc[i])==TRUE&EDGES[[ts]][,2,ls]%in%as.numeric(D4$loc[MATCH[j]])==TRUE),3,ls]+1
        }
        
      } #end loop over captures
    
      ##and now turn the edge list into an association matrix as well to put network in double format
      NET.rows<-as.numeric(factor(rownames(NET[[ts]]),levels=levels(D$loc)))
    
      for (i in 1:length(EDGES[[ts]][,3,ls])){
        NET[[ts]][which(NET.rows%in%EDGES[[ts]][i,1,ls]==TRUE),which(NET.rows%in%EDGES[[ts]][i,2,ls]==TRUE),ls]<-NET[[ts]][which(NET.rows%in%EDGES[[ts]][i,1,ls]==TRUE),which(NET.rows%in%EDGES[[ts]][i,2,ls]==TRUE),ls]+EDGES[[ts]][i,3,ls]
      }
    
    } #end loop over layers
    
    print(paste(ts,"done"))
    
    #end loop over ts/Ws
  }
  
  NODE.EXIST<-data.frame(ids,NODE.EXIST)
  
  results<-list(EDGES,NET,NODE.EXIST,E1,E2)
  
  return(results)
  
}

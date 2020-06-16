
#'MoveNetCreate
#'
#'This function creates dynamic, directed movement networks from capture-mark-recapture datasets using information on the capture locations and times of individuals. Networks connect locations that individuals have moved between within a particular interaction window. The time period for each network, together with the temporal and spatial restrictions on the capture window used to infer a movement can be defined by the user
#'
#'@param data A 5 column dataframe with columns for the ID of the captured individual, the location of its capture (a name or number), the x coordinate of its capture location, the y coordinate of the capture location, and the date of capture
#'@param intwindow The maximum period of time (in days) between the capture of an individual at two different locations for it to be added as an edge to the movemenet network (i.e. if intwindow = 10 then an individual captured at two locations 9 days apart could be considered a movement in the network but two individuals captured 11 days apart could not)
#'@param mindate The start date ("YYYY-MM-DD") of the study (i.e. when you want to build networks from)
#'@param maxdate The end date ("YYYY-MM-DD") of the study (i.e. when you want to build networks until). Please provide as the day after the last day of the study.
#'@param netwindow The period of time over which each network is built in months (i.e. netwindow=12 would correspond to yearly networks)
#'@param overlap The amount of overlap between netwindows in months (i.e. overlap=2 would result in a second network window starting 2 months before the end of the first). Overlap=0 ensures no overlap between successive network windows
#'@param nextonly (TRUE/FALSE). Determines whether a network edge is only created to the next capture of an individual or all captures within the intwindow. Defaults to FALSE

#'@return A list of length 3 containing:
#'\enumerate{
#'    \item the edgelist for the network in each of the netwindows as an array
#'    \item the adjacency matrix for the network in each of the netwindows as an array
#'    \item a matrix indicating which individuals occurred in each netwindow
#'}
#'@examples
#'\dontrun{
#'# load example data
#'data(cmrData)
#'
#'# set parameters for network creation
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
#'}
#'@export


MoveNetCreate<-function(data,intwindow,mindate,maxdate,netwindow,overlap,nextonly=FALSE){

  D<-data
  names(D)<-c("id","loc","x","y","date")
  Jdays<-timeDate::julian(as.Date(D$date),origin=as.Date("1970-01-01"))
  D<-data.frame(D,Jdays)
  D<-D[order(D$Jdays, D$loc,D$id),]
  D$id<-as.factor(D$id)
  D$loc<-as.factor(D$loc)
  D$x<-as.numeric(D$x)
  D$y<-as.numeric(D$y)

  L<-netwindow
  O<-overlap
  start<-timeDate::julian(as.Date(mindate),origin=as.Date("1970-01-01"))
  end<-timeDate::julian(as.Date(maxdate),origin=as.Date("1970-01-01"))
  days<-seq(start,end,1)
  length<-length(days)

  month_seq<-seq(as.Date(mindate),as.Date(maxdate), by = "month")

  starts<-month_seq[seq(1,length(month_seq)-L,L-O)]
  ends<-month_seq[which(month_seq%in%starts)+(L)]

  starts<-timeDate::julian(as.Date(starts),origin=as.Date("1970-01-01"))
  ends<-timeDate::julian(as.Date(ends),origin=as.Date("1970-01-01"))

  #Provide warning message if the netwindows stop early
  if(end-ends[length(ends)]==0){
    print("End of final network window aligns with end of study")
  }
  if(end-ends[length(ends)]>0){
    print(paste0("Final network window stops ",end-ends[length(ends)]," before the end of the study"))
  }

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

  EDGES<-array(0,dim=c(((n.locs-1)*(n.locs)),3,Ws))

  NET<-array(0,dim=c(n.locs,n.locs,Ws))
  colnames(NET)<-locs
  rownames(NET)<-locs

  E1<-rep(locs,each=n.locs-1)

  E1<-factor(E1,levels=levels(locs))
  EDGES[,1,]<-E1

  E2<-locs[2:length(locs)]
  for(i in 2:n.locs){
    E2<-c(as.character(E2),as.character(locs[-i]))
  }

  E2<-factor(E2,levels=levels(locs))
  EDGES[,2,]<-E2

  NODE.EXIST<-matrix(0,nrow=n.locs,ncol=Ws)

  EDGE.EXIST<-matrix(0,nrow=length(EDGES[,1,1]),ncol=Ws)

  # this is the longest step - set up a progress bar
  pb <- progress::progress_bar$new(total = nrow(D2), clear = FALSE)
  pb$tick(0)

  #Less than ends
  for (ts in 1:Ws){

    D3<-D2[which(D2$Jdays>=starts[ts]&D2$Jdays<ends[ts]),]
    D3$id<-factor(D3$id,levels=levels(D2$id))
    D3$loc<-factor(D3$loc,levels=levels(D2$loc))
    n.Caps2<-length(D3$id)

    for (i in 1:n.locs){

      ifelse(locs[i]%in%D3$loc>0,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts]+1,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts])

    }

    for (i in 1:n.Caps2){

      # add progress update
      pb$tick()

      range<-seq(D3$Jdays[i],D3$Jdays[i]+X)
      timematch<-which(D3$Jdays%in%range==TRUE)
      idmatch<-which(D3$id%in%D3$id[i]==TRUE)
      twomatch<-timematch[which(timematch%in%idmatch==TRUE)]
      MATCH<-twomatch[-which(twomatch==i)]

      if(nextonly==TRUE){MATCH<-MATCH[which.min(D3$Jdays[MATCH]-D3$Jdays[i])]}

      for (j in 1:length(MATCH)){
        EDGES[which(EDGES[,1,ts]%in%as.numeric(D3$loc[i])==TRUE&EDGES[,2,ts]%in%as.numeric(D3$loc[MATCH[j]])==TRUE),3,ts]<-EDGES[which(EDGES[,1,ts]%in%as.numeric(D3$loc[i])==TRUE&EDGES[,2,ts]%in%as.numeric(D3$loc[MATCH[j]])==TRUE),3,ts]+1
      }

    }

    NET.rows<-as.numeric(factor(rownames(NET),levels=levels(D$loc)))

    EDGES.tmp<-EDGES[which(EDGES[,3,ts]>0),,ts]

    if(is.matrix(EDGES.tmp)){

      for (i in 1:length(EDGES.tmp[,3])){
        NET[which(NET.rows%in%EDGES.tmp[i,1]==TRUE),which(NET.rows%in%EDGES.tmp[i,2]==TRUE),ts]<-NET[which(NET.rows%in%EDGES.tmp[i,1]==TRUE),which(NET.rows%in%EDGES.tmp[i,2]==TRUE),ts]+EDGES.tmp[i,3]
      }

    }

    if(is.vector(EDGES.tmp)){
      NET[which(NET.rows%in%EDGES.tmp[1]==TRUE),which(NET.rows%in%EDGES.tmp[2]==TRUE),ts]<-NET[which(NET.rows%in%EDGES.tmp[1]==TRUE),which(NET.rows%in%EDGES.tmp[2]==TRUE),ts]+EDGES.tmp[3]
    }

    #end loop over ts/Ws

  }

  NODE.EXIST<-data.frame(locs,NODE.EXIST)

  results<-list(EDGES,NET,NODE.EXIST)

  return(results)

}

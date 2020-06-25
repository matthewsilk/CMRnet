#'DynamicNetCreateHi
#'
#'This function creates dynamic social networks from capture-mark-recapture datasets using information on co-capture
#'The time period for each network, together with the temporal and spatial restrictions on the co-capture window can be defined by the user. This function includes times and dates and can be used for data with higher temporal resolution than DynamicNetCreate()
#'
#'@param data A 5 column dataframe with columns for the ID of the captured individual, the location of its capture (a name or number), the x coordinate of its capture location, the y coordinate of the capture location, and the date and time of capture
#'@param intwindow The maximum period of time (in minutes) between two co-captures (i.e. if \code{intwindow = 60} then two individuals captured 60 minutes apart could be considered co-captured but two individuals captured 61 days apart couldn't)
#'@param mindate The start date (format = \code{"YYYY-MM-DD hh:mm:ss"}) of the study (i.e. when you want to build networks from)
#'@param maxdate The end date (format = \code{"YYYY-MM-DD hh:mm:ss"}) of the study (i.e. when you want to build networks until). Please provide as the day after the last day of the study.
#'@param netwindow The period of time over which each network is built in days (i.e. \code{netwindow=30} would correspond to monthly networks)
#'@param overlap The amount of overlap between netwindows in days (i.e. \code{overlap = 5} would result in a second network window starting 5 days before the end of the first). When \code{overlap=0}, there is no overlap between successive network windows
#'@param spacewindow The maximum distance between locations that can be classed as a co-capture (calculated using the coordinate system provided in in the input data). Best used when multiple capture locations occur very close together
#'@param index Defaults to FALSE. If FALSE edges are weighted by the number of co-captures. If TRUE then edges are weighted by the number by an association index (no. co-captures / (total captures of A + total captures of B - no. of captures))
#'@return A list of length 3 containing: \cr
#'\enumerate{
#'    \item The edgelist for the social network in each of the netwindows as an array
#'    \item The adjacency matrix for the social network in each of the netwindows as an array
#'    \item A matrix indicating which individuals occurred in each netwindow
#'}
#'@examples
#'\dontrun{
#'# load example data
#'data(cmrData3)
#'
#'# set parameters for network creation
#'mindate<-"2019-12-01 00:00:00"
#'maxdate<-"2020-07-01 00:00:00"
#'intwindow<-24*60
#'netwindow<-20
#'overlap<-2
#'spacewindow<-0
#'
#'# create dynamic network
#'netdat<-DynamicNetCreate(data=cmrData,
#'intwindow=intwindow,
#'mindate=mindate,
#'maxdate=maxdate,
#'netwindow=netwindow,
#'overlap=overlap,
#'spacewindow=spacewindow,
#'index=FALSE)
#'}
#'@export

DynamicNetCreateHi<-function(data,intwindow,mindate,maxdate,netwindow,overlap,spacewindow,index=FALSE){

  #Add a column with Julian Dates
  D<-data
  names(D)<-c("id","loc","x","y","date")
  Jdays<-timeDate::julian(as.Date(D$date),origin=as.Date("1970-01-01"))
  times<-chron::times(strftime(D$date,"%H:%M:%S"))
  hours<-chron::hours(times)
  minutes<-chron::minutes(times)
  seconds<-chron::seconds(times)
  tts<-Jdays*24*60*60+hours*60*60+minutes*60+seconds
  D<-data.frame(D,Jdays,tts)
  D<-D[order(D$tts, D$loc,D$id),]
  D$id<-as.factor(D$id)
  D$loc<-as.factor(D$loc)
  D$x<-as.numeric(D$x)
  D$y<-as.numeric(D$y)

  L<-netwindow*24*60*60
  O<-overlap*24*60*60
  startdate<-timeDate::julian(as.Date(mindate),origin=as.Date("1970-01-01"))
  starttime<-chron::times(strftime(mindate,"%H:%M:%S"))
  enddate<-timeDate::julian(as.Date(maxdate),origin=as.Date("1970-01-01"))
  endtime<-chron::times(strftime(maxdate,"%H:%M:%S"))
  start<-startdate*24*60*60+chron::hours(starttime)*60*60+chron::minutes(starttime)*60+chron::seconds(starttime)
  end<-enddate*24*60*60+chron::hours(endtime)*60*60+chron::minutes(endtime)*60+chron::seconds(endtime)
  seconds<-seq(start,end,1)
  length<-length(seconds)

  starts<-seq(start,end-L,L-O)
  ends<-starts+L

  # Provide warning message if the netwindows stop early
  if(end-ends[length(ends)]==0){
    print("End of final network window aligns with end of study")
  }
  if(end-ends[length(ends)]>0){
    print(paste0("Final network window stops ",end-ends[length(ends)],"days before the end of the study"))
  }

  # Counts the number of windows over which networks are built
  Ws<-length(starts)

  #size of sliding window
  X<-intwindow*60

  #get only the data of interest
  #Having less than end
  D2<-D[which(D$tts>=start&D$tts<end),]

  #extract unique individuals and record how many there are
  ids<-sort(unique(D2$id))
  n.ids<-length(ids)

  #extract unique locations and record how many there are
  locs<-sort(unique(D2$loc))
  n.locs<-length(locs)

  #calculate distance matrix between locations
  locdat<-stats::aggregate(D2[,2:4],by=list(D2$loc),unique)[,2:4]
  locmat<-as.matrix(stats::dist(locdat[,2:3]))
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

  NODE.EXIST<-matrix(0,nrow=n.ids,ncol=Ws)

  EDGE.EXIST<-matrix(0,nrow=length(EDGES[,1,1]),ncol=Ws)

  # this is the longest step - set up a progress bar
  pb <- progress::progress_bar$new(total = nrow(D2), clear = FALSE)
  pb$tick(0)

  #Less than ends
  for (ts in 1:Ws){

    D3<-D2[which(D2$tts>=starts[ts]&D2$tts<ends[ts]),]
    D3$id<-factor(D3$id,levels=levels(D2$id))
    n.Caps2<-length(D3$id)

    #loop through IDs and record whether each one was recorded
    #in this time period with a binary response
    for (i in 1:n.ids){

      ifelse(ids[i]%in%D3$id>0,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts]+1,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts])
      #print(paste(ts,"-",i,"-tickB"))
    }

    if(nrow(D3)>0){

    for (i in 1:n.Caps2){

      pb$tick()

      range<-seq(D3$tts[i]-X,D3$tts[i]+X)
      timematch<-which(D3$tts%in%range==TRUE)
      spacematch<-which(D3$loc%in%D3$loc[i]==TRUE)
      spacematch<-which(D3$loc%in%rownames(locmat2)[which(locmat2[,which(colnames(locmat2)==D3$loc[i])]==TRUE)]==TRUE)
      twomatch<-timematch[which(timematch%in%spacematch==TRUE)]
      MATCH<-twomatch[-which(twomatch==i)]

      for (j in 1:length(MATCH)){
        EDGES[which(EDGES[,1,ts]%in%as.numeric(D3$id[i])==TRUE&EDGES[,2,ts]%in%as.numeric(D3$id[MATCH[j]])==TRUE),3,ts]<-EDGES[which(EDGES[,1,ts]%in%as.numeric(D3$id[i])==TRUE&EDGES[,2,ts]%in%as.numeric(D3$id[MATCH[j]])==TRUE),3,ts]+1
      }
    }

    }

    if(index==TRUE){
      for(i in 1:nrow(EDGES[,,ts])){
        if(EDGES[i,3,ts]>0){
          EDGES[i,3,ts]<-EDGES[i,3,ts]/(sum(D3$id==EDGES[i,1,ts])+sum(D3$id==EDGES[i,2,ts])-EDGES[i,3,ts])
        }
      }
    }

    ##and now turn the edge list into an association matrix as well to put network in double format
    NET.rows<-as.numeric(factor(rownames(NET),levels=levels(D$id)))

    if(sum(EDGES[,3,ts])>0){

    EDGES.tmp<-EDGES[which(EDGES[,3,ts]>0),,ts]

    if(is.matrix(EDGES.tmp)){

      for (i in 1:length(EDGES.tmp[,3])){
        NET[which(NET.rows%in%EDGES.tmp[i,1]==TRUE),which(NET.rows%in%EDGES.tmp[i,2]==TRUE),ts]<-NET[which(NET.rows%in%EDGES.tmp[i,1]==TRUE),which(NET.rows%in%EDGES.tmp[i,2]==TRUE),ts]+EDGES.tmp[i,3]
        NET[which(NET.rows%in%EDGES.tmp[i,2]==TRUE),which(NET.rows%in%EDGES.tmp[i,1]==TRUE),ts]<-NET[which(NET.rows%in%EDGES.tmp[i,2]==TRUE),which(NET.rows%in%EDGES.tmp[i,1]==TRUE),ts]+EDGES.tmp[i,3]
      }

    }

    if(is.vector(EDGES.tmp)){
      NET[which(NET.rows%in%EDGES.tmp[1]==TRUE),which(NET.rows%in%EDGES.tmp[2]==TRUE),ts]<-NET[which(NET.rows%in%EDGES.tmp[1]==TRUE),which(NET.rows%in%EDGES.tmp[2]==TRUE),ts]+EDGES.tmp[3]
      NET[which(NET.rows%in%EDGES.tmp[2]==TRUE),which(NET.rows%in%EDGES.tmp[1]==TRUE),ts]<-NET[which(NET.rows%in%EDGES.tmp[2]==TRUE),which(NET.rows%in%EDGES.tmp[1]==TRUE),ts]+EDGES.tmp[3]
    }

    }

    #end loop over ts/Ws

  }

  NODE.EXIST<-data.frame(ids,NODE.EXIST)

  results<-list(EDGES,NET,NODE.EXIST)

  return(results)

  #end function
}

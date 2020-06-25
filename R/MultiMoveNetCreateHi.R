
#'MultiMoveNetCreate
#'
#'This function creates dynamic, directed, multiplex movement networks from capture-mark-recapture datasets using information on the capture locations and times of individuals. Multiplex networks connect locations that individuals have moved between within a particular interaction window with different layers being defined by the user. The time period for each network, together with the temporal and spatial restrictions on the capture window used to infer a movement can be defined by the user
#'
#'@param data A 6 column dataframe with columns for the ID of the captured individual, the location of its capture (a name or number), the x coordinate of its capture location, the y coordinate of the capture location, the date and time of capture and information to define multiplex network layers. The "Layers" column can consist of any unique identifiers (e.g. if layers representing movements by males and females are used then they could be represented by "M" and "F" or 1 and 2). If the user wants a layer per individual then the "Layers" column can simply be a copy of the individual ID column.
#'@param intwindow The maximum period of time (in minutes) between two co-captures (i.e. if \code{intwindow = 60} then two individuals captured 60 minutes apart could be considered co-captured but two individuals captured 61 days apart couldn't)
#'@param mindate The start date (format = \code{"YYYY-MM-DD hh:mm:ss"}) of the study (i.e. when you want to build networks from)
#'@param maxdate The end date (format = \code{"YYYY-MM-DD hh:mm:ss"}) of the study (i.e. when you want to build networks until). Please provide as the day after the last day of the study. Please provide as the day after the last day of the study.
#'@param netwindow The period of time over which each network is built in days (i.e. \code{netwindow=30} would correspond to monthly networks)
#'@param overlap The amount of overlap between netwindows in days (i.e. \code{overlap = 5} would result in a second network window starting 5 days before the end of the first). When \code{overlap=0}, there is no overlap between successive network windows
#'@param nextonly (TRUE/FALSE). Determines whether a network edge is only created to the next capture of an individual or all captures within the intwindow. Defaults to FALSE
#'@param index Defaults to FALSE. If FALSE edges are weighted by the number of movements. If TRUE then edges are weighted by the number of movements divided by the number of captures in a group
#'@return A list of length 3 containing:
#'\enumerate{
#'    \item A list of edgelists (the same length as the number of network windows) containing the multiplex network for each of the netwindows as an array
#'    \item a list of adjacency matrices (the same length as the number of ntwork windows) containing the multiplex network for each of the netwindows as an array
#'    \item a matrix indicating which individuals occurred in each netwindow
#'}
#'@export

MultiMoveNetCreateHi<-function(data,intwindow,mindate,maxdate,netwindow,overlap,nextonly=FALSE,index=FALSE){

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

  NODE.EXIST<-array(0,dim=c(n.locs,Ws,n.layers))

  # this is the longest step - set up a progress bar
  pb <- progress::progress_bar$new(total = nrow(D2), clear = FALSE)
  pb$tick(0)

  #Less than ends
  for (ts in 1:Ws){

    D3<-D2[which(D2$tts>=starts[ts]&D2$tts<ends[ts]),]
    D3$id<-factor(D3$id,levels=levels(D2$id))
    D3$loc<-factor(D3$loc,levels=levels(D2$loc))
    D3$layer<-factor(D3$layer,levels=levels(D2$layer))
    n.Caps2<-length(D3$id)



    for(ls in 1:n.layers){

      D4<-D3[which(D3$layer==layers[ls]),]
      D4$id<-factor(D4$id,levels=levels(D3$id))
      D4$loc<-factor(D4$loc,levels=levels(D3$loc))
      D4$layer<-factor(D4$layer,levels=levels(D3$layer))
      n.Caps2<-length(D4$id)

      #loop through IDs and record whether each one was recorded
      #in this time period/layer with a binary response
      for (i in 1:n.locs){

        ifelse(locs[i]%in%D4$loc>0,NODE.EXIST[i,ts,ls]<-NODE.EXIST[i,ts,ls]+1,NODE.EXIST[i,ts,ls]<-NODE.EXIST[i,ts,ls])

      }

      for (i in 1:n.Caps2){

        pb$tick()

        range<-seq(D4$tts[i],D4$tts[i]+X)
        timematch<-which(D4$tts%in%range==TRUE)
        idmatch<-which(D4$id%in%D4$id[i]==TRUE)
        twomatch<-timematch[which(timematch%in%idmatch==TRUE)]
        MATCH<-twomatch[-which(twomatch==i)]

        if(nextonly==TRUE){MATCH<-MATCH[which.min(D4$tts[MATCH]-D4$tts[i])]}

        for (j in 1:length(MATCH)){
          EDGES[[ts]][which(EDGES[[ts]][,1,ls]%in%as.numeric(D4$loc[i])==TRUE&EDGES[[ts]][,2,ls]%in%as.numeric(D4$loc[MATCH[j]])==TRUE),3,ls]<-EDGES[[ts]][which(EDGES[[ts]][,1,ls]%in%as.numeric(D4$loc[i])==TRUE&EDGES[[ts]][,2,ls]%in%as.numeric(D4$loc[MATCH[j]])==TRUE),3,ls]+1
        }

      } #end loop over captures

      if(index==TRUE){
        for(i in 1:nrow(EDGES[,,ts])){
          if(EDGES[i,3,ts]>0){
            EDGES[[ts]][i,3,ls]<-EDGES[[ts]][i,3,ls]/(sum(as.numeric(D3$loc)==EDGES[[ts]][i,1,ls]))
          }
        }
      }

      ##and now turn the edge list into an association matrix as well to put network in double format
      NET.rows<-as.numeric(factor(rownames(NET[[ts]]),levels=levels(D$loc)))

      if(sum(EDGES[[ts]][,3,ls])>0){

      EDGES.tmp<-EDGES[[ts]][which(EDGES[[ts]][,3,ls]>0),,ls]

      if(is.matrix(EDGES.tmp)){

        for (i in 1:length(EDGES.tmp[,3])){
          NET[[ts]][which(NET.rows%in%EDGES.tmp[i,1]==TRUE),which(NET.rows%in%EDGES.tmp[i,2]==TRUE),ls]<-NET[[ts]][which(NET.rows%in%EDGES.tmp[i,1]==TRUE),which(NET.rows%in%EDGES.tmp[i,2]==TRUE),ls]+EDGES.tmp[i,3]
        }

      }

      if(is.vector(EDGES.tmp)){
        NET[[ts]][which(NET.rows%in%EDGES.tmp[1]==TRUE),which(NET.rows%in%EDGES.tmp[2]==TRUE),ls]<-NET[[ts]][which(NET.rows%in%EDGES.tmp[1]==TRUE),which(NET.rows%in%EDGES.tmp[2]==TRUE),ls]+EDGES.tmp[3]
      }

      }

    } #end loop over layers

    #end loop over ts/Ws

  }

  rownames(NODE.EXIST)<-locs

  results<-list(EDGES,NET,NODE.EXIST)

  return(results)

}


#'DatastreamPermSpat
#'This function creates randomised movement networks for each network window using datastream permutations with user-defined restrictions (to constrain swaps according to temporal or individual ID-based restrictions)
#'
#'@param data A 5 column dataframe with columns for the ID of the captured individual, the location of its capture (a name or number), the x coordinate of its capture location, the y coordinate of the capture location, and the date of capture
#'@param intwindow The maximum period of time (in days) between two co-captures (i.e. if intwindow = 10 then two individuals captured 10 days apart could be considered co-captured but two indivviduals captured 11 days apart couldn't)
#'@param mindate The start date ("YYYY-MM-DD") of the study (i.e. when you want to build networks from)
#'@param maxdate The end date ("YYYY-MM-DD") of the study (i.e. when you want to build networks until). Please provide as the day after the last day of the study.
#'@param netwindow The period of time over which each network is built in months (i.e. netwindow=12 would correspond to yearly networks)
#'@param overlap The amount of overlap between netwindows in months (i.e. overlap=2 would result in a second network window starting 2 months before the end of the first). Overlap=0 ensures no overlap between successive network windows
#'@param nextonly (TRUE/FALSE).Determines whether a network edge is only created to the next capture of an individual or all captures within the intwindow. Defaults to false
#'@param same.time (TRUE/FALSE) Whether swaps should be restricted to only occur betwen individuals trapped on the same date or not
#'@param time.restrict Provided as a number of months. Imposes time restrictions on when swaps can take place so that individuals can only be swapped with those a fixed time before or after being captured
#'@param same.id (TRUE/FALSE) Whether swaps should be restricted to only be between captures of the same individual
#'@param n.swaps The number of swaps between each random network being extracted (e.g. n.swaps = 10 would equate to 10 swaps taking place between each random network being saved)
#'@param n.rand The number of randomised networks to be generated
#'@param n.burnin The number of swaps to discard as burn-in before the first random network is created. The total number of swaps conducted is thus n.burnin+n.swaps*n.rand
#'

#'@output A list of length 3 with elements corresponding to 1) Randomised edgelist list: a list of with the same number of elements at the number of network windows, with each element containing an array of the randomised edgelists, 2) Randomised adjacency matrix list: a list of with the same number of elements at the number of network windows, with each element containing an array of the randomised adjacency matrices, and 3) a matrix identifying whether a location was present (i.e. had at least one individual captured there) in each network window.

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

DatastreamPermSpat<-function(data,intwindow,mindate,maxdate,netwindow,overlap,nextonly=FALSE,same.time,time.restrict,same.id,n.swaps,n.rand,n.burnin,warn.thresh){

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

        if(nextonly==TRUE){MATCH<-MATCH[which.min(D3$Jdays[MATCH]-D3$Jdays[i])]}

        for (j in 1:length(MATCH)){
          EDGES[[ts]][which(EDGES[[ts]][,1,r]%in%as.numeric(D3$loc[i])==TRUE&EDGES[[ts]][,2,r]%in%as.numeric(D3$loc[MATCH[j]])==TRUE),3,r]<-EDGES[[ts]][which(EDGES[[ts]][,1,r]%in%as.numeric(D3$loc[i])==TRUE&EDGES[[ts]][,2,r]%in%as.numeric(D3$loc[MATCH[j]])==TRUE),3,r]+1
        }

      } #end loop over n.caps

      ##and now turn the edge list into an association matrix as well to put network in double format
      NET.rows<-as.numeric(factor(rownames(NET[[ts]]),levels=levels(D$loc)))

      for (i in 1:length(EDGES[[ts]][,3,r])){
        NET[[ts]][which(NET.rows%in%EDGES[[ts]][i,1,r]==TRUE),which(NET.rows%in%EDGES[[ts]][i,2,r]==TRUE),r]<-NET[[ts]][which(NET.rows%in%EDGES[[ts]][i,1,r]==TRUE),which(NET.rows%in%EDGES[[ts]][i,2,r]==TRUE),r]+EDGES[[ts]][i,3,r]
      }

    } #end r loop over randomisations

    #end loop over ts/Ws
  }

  NODE.EXIST<-data.frame(ids,NODE.EXIST)

  results<-list(EDGES,NET,NODE.EXIST)

  return(results)

}

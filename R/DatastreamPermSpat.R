#'DatastreamPermSpat
#'
#'This function creates randomised movement networks for each network window using datastream permutations with user-defined restrictions (to constrain swaps according to temporal or individual ID-based restrictions)
#'
#'@param data A 5 column dataframe with columns for the ID of the captured individual, the location of its capture (a name or number), the x coordinate of its capture location, the y coordinate of the capture location, and the date of capture
#'@param intwindow The maximum period of time (in days) between two co-captures (i.e. if \code{intwindow = 10} then two individuals captured 10 days apart could be considered co-captured but two individuals captured 11 days apart couldn't)
#'@param mindate The start date (format = \code{"YYYY-MM-DD"}) of the study (i.e. when you want to build networks from)
#'@param maxdate The end date (format = \code{"YYYY-MM-DD"}) of the study (i.e. when you want to build networks until). Please provide as the day after the last day of the study.
#'@param netwindow The period of time over which each network is built in months (i.e. \code{netwindow=12} would correspond to yearly networks)
#'@param overlap The amount of overlap between netwindows in months (i.e. \code{overlap=2} would result in a second network window starting 2 months before the end of the first). When \code{overlap=0}, there is no overlap between successive network windows
#'@param same.time (TRUE/FALSE) Whether swaps should be restricted to only occur trapping events on the same date or not
#'@param time.restrict Provided as a number of months. Imposes time restrictions on when swaps can take place so that locations can only be swapped with those a fixed time before or after being captured
#'@param spat.restrict Provided on the same scale as the coordinates in the input dataset. Imposes space restrictions on when swaps can take place so that locations can only be swapped with those captued within a fixed distance
#'@param same.id (TRUE/FALSE) Whether swaps should be restricted to only be between captures of the same individual
#'@param n.swaps The number of swaps between each random network being extracted (e.g. \code{n.swaps = 10} would equate to 10 swaps taking place between each random network being saved)
#'@param n.rand The number of randomised networks to be generated
#'@param burnin (TRUE/FALSE) Whether burnin is required
#'@param n.burnin The number of swaps to discard as burn-in before the first random network is created. The total number of swaps conducted is thus n.burnin+n.swaps*n.rand
#'@param warn.thresh The number of times no matches are found (i.e. constraints on randomisations are too restrictive) before the function is stopped and an error message returned
#'@param nextonly (TRUE/FALSE). Determines whether a network edge is only created to the next capture of an individual or all captures within the intwindow. Defaults to FALSE
#'@param iter (TRUE/FALSE) Whether iterative randomisations are being used. If TRUE then D.rand is also returned
#'
#'@return If \code{iter=TRUE} then a list of length 3 is returned with elements corresponding to:
#'\enumerate{
#'    \item The randomised dataset (for feeding back into the next permutation)
#'    \item Randomised adjacency matrix list: a list of with the same number of elements at the number of network windows, with each element containing an array of the randomised adjacency matrices
#'    \item a matrix identifying whether a location was present (i.e. had at least one individual captured there) in each network window
#'}
#'If \code{iter=FALSE} then a list of length 2 is returned with elements corresponding to:
#'\enumerate{
#'    \item Randomised adjacency matrix list: a list of with the same number of elements at the number of network windows, with each element containing an array of the randomised adjacency matrices
#'    \item A matrix identifying whether a location was present (i.e. had at least one individual captured there) in each network window. The edge list is not provided due to to the memory that providing this and the list of matrix arrays would require.
#'}
#'@examples
#'\dontrun{
#'# load example data
#'data(cmrData)
#'
#'# set parameters
#'mindate<-"2010-01-01"
#'maxdate<-"2015-01-01"
#'intwindow<-60
#'netwindow<-12
#'overlap<-0
#'
#'# create networks
#'movenetdat<-MoveNetCreate(data=cmrData,
#'intwindow=intwindow,
#'mindate=mindate,
#'maxdate=maxdate,
#'netwindow=netwindow,
#'overlap=overlap,
#'nextonly=TRUE)
#'
#'# set additional parameters
#'same.time=FALSE
#'time.restrict=6
#'same.id=FALSE
#'n.swaps=10
#'n.rand=100
#'n.burnin=100
#'warn.thresh=100
#'
#'# perform permutations
#'Rs<-DatastreamPermSpat(data=cmrData,
#'intwindow,
#'mindate,
#'maxdate,
#'netwindow,
#'overlap,
#'spacewindow,
#'same.time,
#'time.restrict,
#'spat.restrict="n"
#'same.id,
#'n.swaps,
#'n.rand,
#'burnin=TRUE,
#'n.burnin,
#'warn.thresh,
#'nextonly=TRUE,
#'iter=FALSE)
#'}

#'@export

DatastreamPermSpat<-function(data,intwindow,mindate,maxdate,netwindow,overlap,nextonly=FALSE,same.time,time.restrict,spat.restrict,same.id,n.swaps,n.rand,burnin,n.burnin,warn.thresh,iter){

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

  locdat<-stats::aggregate(D2[,2:4],by=list(D2$loc),unique)[,2:4]
  locmat<-as.matrix(stats::dist(locdat[,2:3]))
  if(spat.restrict!="n"){
    locmat2<-locmat<spat.restrict
    rownames(locmat2)<-colnames(locmat2)<-locs
  } else{locmat2<-locmat}


  n.caps<-length(D2[,1])

  NET<-list()

  for(ts in 1:Ws){

    NET[[ts]]<-array(0,dim=c(n.locs,n.locs,n.rand))
    colnames(NET[[ts]])<-locs
    rownames(NET[[ts]])<-locs
  }

  NODE.EXIST<-matrix(0,nrow=n.locs,ncol=Ws)

  rands.out<-list()

  # set up progress bar
  pb <- progress::progress_bar$new(total = Ws*n.rand, clear = FALSE)
  pb$tick(0)

  #Less than ends
  for (ts in 1:Ws){

    #create edgelist (temporary)
    E1<-rep(locs,each=n.locs-1)
    E1<-factor(E1,levels=levels(locs))
    E2<-locs[2:length(locs)]
    for(i in 2:n.locs){
      E2<-c(as.character(E2),as.character(locs[-i]))
    }
    E2<-factor(E2,levels=levels(locs))
    EDGES<-array(0,dim=c(((n.locs-1)*(n.locs)),3,n.rand))
    EDGES[,1,]<-E1
    EDGES[,2,]<-E2

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

    rands<-cmrPermSpat(D=D3,locmat=locmat2,same.time=same.time,time.restrict=time.restrict,spat.restrict=spat.restrict,same.id=same.id,n.swaps=n.swaps,n.rand=n.rand,burnin=burnin,n.burnin=n.burnin,warn.thresh=warn.thresh)

    if(iter==TRUE){rands.out[[ts]]<-as.data.frame(rands[[1]][,1:5])}

    for(r in 1:length(rands)){

      pb$tick()

      for (i in 1:n.Caps2){
        D3<-rands[[r]]
        range<-seq(D3$Jdays[i],D3$Jdays[i]+X)
        timematch<-which(D3$Jdays%in%range==TRUE)
        idmatch<-which(D3$id%in%D3$id[i]==TRUE)
        twomatch<-timematch[which(timematch%in%idmatch==TRUE)]
        MATCH<-twomatch[-which(twomatch==i)]

        if(nextonly==TRUE){MATCH<-MATCH[which.min(D3$Jdays[MATCH]-D3$Jdays[i])]}

        for (j in 1:length(MATCH)){
          EDGES[which(EDGES[,1,r]%in%as.numeric(D3$loc[i])==TRUE&EDGES[,2,r]%in%as.numeric(D3$loc[MATCH[j]])==TRUE),3,r]<-EDGES[which(EDGES[,1,r]%in%as.numeric(D3$loc[i])==TRUE&EDGES[,2,r]%in%as.numeric(D3$loc[MATCH[j]])==TRUE),3,r]+1
        }

      } #end loop over n.caps

      ##and now turn the edge list into an association matrix as well to put network in double format
      NET.rows<-as.numeric(factor(rownames(NET[[ts]]),levels=levels(D$loc)))

      EDGES.tmp<-EDGES[which(EDGES[,3,r]>0),,r]

      for (i in 1:length(EDGES.tmp[,3])){
        NET[[ts]][which(NET.rows%in%EDGES.tmp[i,1]==TRUE),which(NET.rows%in%EDGES.tmp[i,2]==TRUE),r]<-NET[[ts]][which(NET.rows%in%EDGES.tmp[i,1]==TRUE),which(NET.rows%in%EDGES.tmp[i,2]==TRUE),r]+EDGES.tmp[i,3]
      }

    } #end r loop over randomisations

    #end loop over ts/Ws
  }

  NODE.EXIST<-data.frame(locs,NODE.EXIST)


  if(iter==FALSE){results<-list(NET,NODE.EXIST)}
  if(iter==TRUE){
    rands.out2<-do.call("rbind", rands.out)
    results<-list(rands.out2,NET,NODE.EXIST)
  }

  return(results)

}

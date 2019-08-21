#'DatastreamPermSoc
#'
#'This function creates randomised social networks for each network window using datastream permutations with user-defined restrictions (to constrain swaps according to temporal or spatial windows)
#'
#'@param data A 5 column dataframe with columns for the ID of the captured individual, the location of its capture (a name or number), the x coordinate of its capture location, the y coordinate of the capture location, and the date of capture
#'@param intwindow The maximum period of time (in days) between two co-captures (i.e. if \code{intwindow = 10} then two individuals captured 10 days apart could be considered co-captured but two individuals captured 11 days apart would not be)
#'@param mindate The start date (format = \code{"YYYY-MM-DD"}) of the study (i.e. when you want to build networks from)
#'@param maxdate The end date (format = \code{"YYYY-MM-DD"}) of the study (i.e. when you want to build networks until). Please provide as the day after the last day of the study.
#'@param netwindow The period of time (in months) over which each network is built(i.e. \code{netwindow=12} would correspond to yearly networks)
#'@param overlap The amount of overlap between networks in months (i.e. \code{overlap=2} would result in a second network window starting 2 months before the end of the first). When \code{overlap=0}, there is no overlap between successive network windows
#'@param spacewindow The maximum distance between locations that can be classed as a co-capture (calculated using the coordinate system provided in the input data). Best used when multiple capture locations occur very close together
#'@param same.time (TRUE/FALSE) Whether swaps should be restricted to only occur betwen individuals trapped on the same date or not
#'@param time.restrict Provided as a number of months. Imposes time restrictions on when swaps can take place so that individuals can only be swapped with those a fixed time before or after being captured
#'@param same.spat (TRUE/FALSE) Whether swaps should be restricted to only occur between individuals trapped at the same location
#'@param spat.restrict Provided on the same scale as the coordinates in the input dataset. Imposes space restrictions on when swaps can take place so that individuals can only be swapped with those captued within a fixed distance
#'@param n.swaps The number of swaps between each random network being extracted (e.g. \code{n.swaps = 10} would equate to 10 swaps taking place between each random network being saved)
#'@param n.rand The number of randomised networks to be generated
#'@param burnin (TRUE/FALSE) Whether burnin is required
#'@param n.burnin The number of swaps to discard as burn-in before the first random network is created. The total number of swaps conducted is thus n.burnin+n.swaps*n.rand
#'@param iter (TRUE/FALSE) Whether iterative randomisations are being used. If TRUE then D.rand is also returned
#'@return If \code{iter=TRUE} then a list of length 3 with elements corresponding to: \cr
#'\enumerate{
#'    \item The randomised dataset (for feeding back into the next permutation)
#'    \item Randomised adjacency matrix list: a list of with the same number of elements at the number of network windows, with each element containing an array of the randomised adjacency matrices
#'    \item A matrix identifying whether an individual was present in each network window.
#'}
#'
#'If \code{iter=FALSE} then a list of length 2 with elements corresponding to: \cr
#'\enumerate{
#'    \item Randomised adjacency matrix list: a list of with the same number of elements at the number of network windows, with each element containing an array of the randomised adjacency matrices
#'    \item A matrix identifying whether an individual was present in each network window. The edge list is not provided due to to the memory that providing this and the list of matrix arrays would require.
#'}
#'@examples
#' # load in data
#' data(cmr_dat)
#'
#' # set parameters
#' mindate<-"2010-01-01"
#' maxdate<-"2015-01-01"
#' intwindow<-60
#' netwindow<-12
#' overlap<-0
#' spacewindow<-0
#'
#' # create network
#' netdat<-DynamicNetCreate(data=cmr_dat,intwindow=intwindow,mindate=mindate,maxdate=maxdate,netwindow=netwindow,overlap=overlap,spacewindow=spacewindow)
#'
#'same.time=FALSE
#'time.restrict=6
#'same.spat=FALSE
#'spat.restrict="n"
#'n.swaps=10
#'n.rand=100
#'n.burnin=100
#'
#'Rs <- DatastreamPermSoc(data=cmr_dat,intwindow,mindate,maxdate,netwindow,overlap,spacewindow,same.time,time.restrict,same.spat,spat.restrict,n.swaps,n.rand,burnin=TRUE,n.burnin,iter=FALSE)

#'@export

DatastreamPermSoc<-function(data, intwindow, mindate, maxdate, netwindow, overlap, spacewindow, same.time, time.restrict, same.spat, spat.restrict, n.swaps, n.rand, burnin, n.burnin, iter){

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

  locdat<-stats::aggregate(D2[,2:4],by=list(D2$loc),unique)[,2:4]
  locmat<-as.matrix(stats::dist(locdat[,2:3]))
  locmat2<-locmat<spat.restrict
  rownames(locmat2)<-colnames(locmat2)<-locs

  n.caps<-length(D2[,1])

  NET<-list()

  for(ts in 1:Ws){
    NET[[ts]]<-array(0,dim=c(n.ids,n.ids,n.rand))
    colnames(NET[[ts]])<-ids
    rownames(NET[[ts]])<-ids
  }

  NODE.EXIST<-matrix(0,nrow=n.ids,ncol=Ws)

  rands.out<-list()

  #Less than ends
  for (ts in 1:Ws){

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
    EDGES<-array(0,dim=c(((n.ids-1)*(n.ids))/2,3,n.rand))
    EDGES[,1,]<-E1
    EDGES[,2,]<-E2

    D3<-D2[which(D2$Jdays>=starts[ts]&D2$Jdays<ends[ts]),]
    D3$id<-factor(D3$id,levels=levels(D2$id))
    n.Caps2<-length(D3$id)

    #loop through IDs and record whether each one was recorded
    #in this time period with a binary response
    for (i in 1:n.ids){

      ifelse(ids[i]%in%D3$id>0,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts]+1,NODE.EXIST[i,ts]<-NODE.EXIST[i,ts])
      #print(paste(ts,"-",i,"-tickB"))
    }

    rands<-CMRnet::cmrPermSoc(D=D3,locmat=locmat2,same.time=same.time,time.restrict=time.restrict,same.spat=same.spat,spat.restrict=spat.restrict,n.swaps=n.swaps,n.rand=n.rand,burnin=burnin,n.burnin=n.burnin)

    if(iter==TRUE){rands.out[[ts]]<-as.data.frame(rands[[1]][,1:5])}

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
          EDGES[which(EDGES[,1,r]%in%as.numeric(D3$id[i])==TRUE&EDGES[,2,r]%in%as.numeric(D3$id[MATCH[j]])==TRUE),3,r]<-EDGES[which(EDGES[,1,r]%in%as.numeric(D3$id[i])==TRUE&EDGES[,2,r]%in%as.numeric(D3$id[MATCH[j]])==TRUE),3,r]+1
        }

      }

      ##and now turn the edge list into an association matrix as well to put network in double format
      NET.rows<-as.numeric(factor(rownames(NET[[ts]]),levels=levels(D$id)))

      for (i in 1:length(EDGES[,3,r])){
        NET[[ts]][which(NET.rows%in%EDGES[i,1,r]==TRUE),which(NET.rows%in%EDGES[i,2,r]==TRUE),r]<-NET[[ts]][which(NET.rows%in%EDGES[i,1,r]==TRUE),which(NET.rows%in%EDGES[i,2,r]==TRUE),r]+EDGES[i,3,r]
        NET[[ts]][which(NET.rows%in%EDGES[i,2,r]==TRUE),which(NET.rows%in%EDGES[i,1,r]==TRUE),ts]<-NET[[ts]][which(NET.rows%in%EDGES[i,2,r]==TRUE),which(NET.rows%in%EDGES[i,1,r]==TRUE),r]+EDGES[i,3,r]
      }

    } #end r loop over randomisations

    #end loop over ts/Ws
  }

  NODE.EXIST<-data.frame(ids,NODE.EXIST)

  rands.out2<-do.call("rbind", rands.out)

  if(iter==FALSE){results<-list(NET,NODE.EXIST)}
  if(iter==TRUE){
    rands.out2<-do.call("rbind", rands.out)
    results<-list(rands.out2,NET,NODE.EXIST)
  }

  return(results)

}

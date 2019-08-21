#-------------------#
# testing CMRnet ####
#-------------------#

# install profvis if not present
if('profvis' %in% installed.packages() == FALSE){install.packages('profvis')}

# load profvis
library(profvis)

# load in data
data(cmr_dat)
data(cmr_dat2)

# set parameters
mindate <- "2010-01-01"
maxdate <- "2015-01-01"
intwindow <- 60
netwindow <- 12
overlap <- 0
spacewindow <- 0

# test all network create functions ####

# Dynamic networks
netdat <- DynamicNetCreate(
    data = cmr_dat,
    intwindow = intwindow,
    mindate = mindate,
    maxdate = maxdate,
    netwindow = netwindow,
    overlap = overlap,
    spacewindow = spacewindow
  )

# run profiling on DynamicNetCreate
profvis({
  netdat <- DynamicNetCreate(
    data = cmr_dat,
    intwindow = intwindow,
    mindate = mindate,
    maxdate = maxdate,
    netwindow = netwindow,
    overlap = overlap,
    spacewindow = spacewindow
  )
})

# movement networks
movenetdat <- MoveNetCreate(
  data = cmr_dat,
  intwindow = intwindow,
  mindate = mindate,
  maxdate = maxdate,
  netwindow = netwindow,
  overlap = overlap,
  nextonly = TRUE
)

# Multiple movement networks
multimovenetdat <-
  MultiMoveNetCreate(
    data = cmr_dat2,
    intwindow = intwindow,
    mindate = mindate,
    maxdate = maxdate,
    netwindow = netwindow,
    overlap = overlap,
    nextonly = TRUE
  )

# set new params
same.time = FALSE
time.restrict = 6
same.spat = FALSE
spat.restrict = "n"
n.swaps = 10
n.rand = 10
n.burnin = 2

# run permutations
Rs <- DatastreamPermSoc(
    data = cmr_dat,
    intwindow,
    mindate,
    maxdate,
    netwindow,
    overlap,
    spacewindow,
    same.time,
    time.restrict,
    same.spat,
    spat.restrict,
    n.swaps,
    n.rand,
    burnin = TRUE,
    n.burnin,
    iter = FALSE
  )

Rs <-
  DatastreamPermSpat(
    data = cmr_dat,
    intwindow,
    mindate,
    maxdate,
    netwindow,
    overlap,
    spacewindow,
    same.time,
    time.restrict,
    same.id,
    n.swaps,
    n.rand,
    burnin = TRUE,
    n.burnin,
    warn.thresh,
    nextonly = TRUE,
    iter = FALSE
  )

B <- cmrNodeswap(multimovenetdat, n.rand = 1000, multi = TRUE)

# get the error:
# Error in NET[which(NET.rows %in% EDGES[i, 2, r] == TRUE), which(NET.rows %in%  : incorrect number of subscripts

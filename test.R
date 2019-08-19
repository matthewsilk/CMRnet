# test of DatastreamPermSoc

# load in data
data(cmr_dat)

# set parameters
mindate <- "2010-01-01"
maxdate <- "2015-01-01"
intwindow <- 60
netwindow <- 12
overlap <- 0
spacewindow <- 0

# create network
netdat <- DynamicNetCreate(
    data = cmr_dat,
    intwindow = intwindow,
    mindate = mindate,
    maxdate = maxdate,
    netwindow = netwindow,
    overlap = overlap,
    spacewindow = spacewindow
  )

# [1] "stopped0days early"

# set new params
same.time = FALSE
time.restrict = 6
same.spat = FALSE
spat.restrict = "n"
n.swaps = 10
n.rand = 100
n.burnin = 100

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

# get the error:
# Error in NET[which(NET.rows %in% EDGES[i, 2, r] == TRUE), which(NET.rows %in%  : incorrect number of subscripts

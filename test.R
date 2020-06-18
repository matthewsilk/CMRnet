#-------------------#
# testing CMRnet ####
#-------------------#

# install profvis if not present
if('profvis' %in% installed.packages() == FALSE){install.packages('profvis')}

devtools::install_github('matthewsilk/CMRnet')

3# load profvis
library(profvis)
library(CMRnet)

# load in data
data(cmrData)
data(cmrData2)

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
    data = cmrData,
    intwindow = intwindow,
    mindate = mindate,
    maxdate = maxdate,
    netwindow = netwindow,
    overlap = overlap,
    spacewindow = spacewindow
  )

# movement networks
movenetdat <- MoveNetCreate(
  data = cmrData,
  intwindow = intwindow,
  mindate = mindate,
  maxdate = maxdate,
  netwindow = netwindow,
  overlap = overlap,
  nextonly = TRUE
)

# Multiplex movement networks
multimovenetdat <-
  MultiMoveNetCreate(
    data = cmrData2,
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
n.swaps = 5
n.rand = 2
n.burnin = 0
warn.thresh = 100

# run permutations
Rs <- DatastreamPermSoc(
    data = cmrData,
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
    burnin = FALSE,
    n.burnin,
    iter = FALSE,
    buffer=0
  )

same.time=FALSE
time.restrict=6
same.id=FALSE
n.swaps=10
n.rand=100
n.burnin=100
warn.thresh=100

Rs <- DatastreamPermSpat(
  data = cmrData,
  intwindow,
  mindate,
  maxdate,
  netwindow,
  overlap,
  nextonly = FALSE,
  same.time,
  time.restrict,
  same.id,
  n.swaps,
  n.rand,
  burnin = TRUE,
  n.burnin,
  warn.thresh,
  iter = FALSE
)

# check cmrNodeswap with both multi = TRUE and multi = FALSE
A <- cmrNodeswap(movenetdat, n.rand = 1000, multi = FALSE)
B <- cmrNodeswap(multimovenetdat, n.rand = 1000, multi = TRUE)

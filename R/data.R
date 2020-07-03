#' cmrData
#'
#'cmrData can be used to test social and movement network functions of CMRnet. It is an example of how dataframes should be organised to feed into CMRnet
#'
#' @format A data frame with 1191 rows and 5 columns variables:
#' * \code{id} - provides the identity of the individual captured
#' * \code{loc} - the named capture location
#' * \code{x} - geographic location of a capture event (on the horizontal axis)
#' * \code{y} - geographic location of a capture event (on the vertical axis)
#' * \code{date} - the capture date
#' @md
#' @docType data
#' @usage data("cmrData")
"cmrData"

#' cmrData2
#'
#' A dataset that tests \code{\link{MultiMoveNetCreate}}. \code{cmrData2} is an example of how dataframes should be created to feed into CMRnet when there are multiple layered movement networks (i.e. a separate network for males and females).
#' @format A data frame with 1191 rows and 5 columns variables:
#' * `id` - provides the identity of the individual captured
#' * `loc` - the named capture location
#' * `x` - geographic location of a capture event (on the horizontal axis)
#' * `y` - geographic location of a capture event (on the vertical axis)
#' * `date` - the capture date
#' * `layer` - provides information used to define the layers of the multiplex network
#' @docType data
#' @md
#' @usage data("cmrData2")
"cmrData2"

#' cmrData3
#'
#'cmrData3 to test \code{\link{DynamicNetCreateHi}} and \code{\link{MoveNetCreateHi}}. These functions help create short-term social and movement networks.\code{cmrData3} is an example of how dataframes should be created to feed into CMRnet when capture times are relevant as well as capture dates.
#' @format A data frame with 1191 rows and 5 columns variables:
#' * `id` - provides the identity of the individual captured
#' * `loc` - the named capture location
#' * `x` - geographic location of a capture event (on the horizontal axis)
#' * `y` - geographic location of a capture event (on the vertical axis)
#' * `date` - the capture date in POSIXct format
#' @docType data
#' @md
#' @usage data("cmrData3")
"cmrData3"


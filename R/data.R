#' cmrData
#' cmrData
#'
#'cmrData can be used to test social and movement ntwork functions of CMRnet. It contains five columns: the ID column provides the identity of the individual captured, loc is the named capture location, x and y provide the geographic location of a capture event, and date provide the capture date.
#'
#'
#' @format A data frame with five variables: \code{id}, \code{loc},
#'   \code{x}, \code{y} and \code{date}. \code{cmrData} is an example of how dataframe es should be created to feed into CMRnet.
"cmrData"

#' cmrData2
#'
#'cmrData2 to test the multiplex movement network function.It contains six columns: the ID column provides the identity of the individual captured, loc is the named capture location, x and y provide the geographic location of a capture event, date provides the capture date and layer provide information used to define the layers of the multiplex network.
#'
#' @format A data frame with six variables: \code{id}, \code{loc},
#'   \code{x}, \code{y}, \code{date} and \code{layer}. \code{cmrData2} is an example of how dataframes should be created to feed into CMRnet when there are multiple layered movement networks (i.e. a separate network for males and females).
"cmrData2"

#' cmrData3
#'
#'cmrData3 to test short-term social and movement network functions of CMRnet.It contains five columns: the ID column provides the identity of the individual captured, loc is the named capture location, x and y provide the geographic location of a capture event, and date provide the capture date and time.
#'
#' @format A data frame with five variables: \code{id}, \code{loc},
#'   \code{x}, \code{y} and \code{date}. \code{cmrData3} is an example of how dataframes should be created to feed into CMRnet when capture times are relevant as well as capture dates.
"cmrData3"


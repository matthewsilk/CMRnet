
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CMRnet

<!-- badges: start -->

<!-- badges: end -->

CMRnet can be used to construct social and movement networks from
spatio-temporally referenced capture-mark recapture data.

## Installation

You can install the most up to date version of CMRnet from GitHub.

``` r
# install CMRnet
remotes::install_github("matthewsilk/CMRnet")
```

## Package Description

### Data input

CMRnet can be used to construct social and movement networks from
spatio-temporally referenced CMR data. The simplest form consists of a
dataframe with columns for individual ID, capture location, the XY
coordinates of the capture location and the date of capture. The second
form includes an additional “layer” column that can be used to construct
multiplex movement networks separated according to phenotypic
characteristics of individuals in the population. The advantage of the
latter approach is that it makes it possible to consider the role of
different phenotypes (such as sex), life-history stages or even
among-individual variation in structuring the movement network.
Alongside the input data, the both network construction functions
require a number of other arguments in order to construct networks. For
both types of network it is necessary to provide a start date (using the
mindate= argument) and an equivalent end date (using the maxdate=
argument) to define the period of time that social networks should be
constructed over. It is also necessary to set the length of time (in
months ) that each network should be constructed over using the
netwindow= argument. Finally, the overlap= argument (also in months)
determines whether there is any overlap between the periods of time that
networks are constructed over. For example, with netwindow=12 and
overlap=6 the first network window would run for the first year of the
study, and the second network window would run from the beginning of
month 7 to the end of month 18. For social network construction there is
an additional intwindow= argument in the network construction function.
This defines the amount of time that can elapse between two individuals
being captured at a location before an edge is no longer added to the
network. For example, if intwindow=20 then two individuals captured at
the same location 19 days apart will be connected in the social network
but two individuals captured at the same location 21 days apart would
not.

-----

### Social network construction

CMRnet constructs social networks using the assumption that two
individuals captured at the same location on the same or similar dates
are likely to have had the opportunity to interact or come into contact
with each other (Fig. 1). Therefore, it should be seen as an
approximation of the true social network of the population; the approach
may miss edges if individuals interact away from capture locations, and
could also generate false edges if individuals that are co-captured do
not actually interact in a meaningful way (which will depend on the
question being asked; (Carter, Lee, & Marshall, 2015)). This means it is
especially valuable as an approach when either a) social interaction
data at higher resolutions are not available or prohibitively
impractical or costly to obtain, or b) the social system of the study
population is appropriate for an approach that ties social relationships
to spatio-temporal co-occurrence over the timeframe selected. The
approach is thus likely to be most applicable in situations where
individuals live in relatively stable groups centred around a refuge or
food resource. However, it will also be more widely applicable. For
example, the approach could be used to construct the social networks of
flocks of colour-ringed birds if the interaction window was set to be as
short as possible so that edges were only added if individuals visited
the same site at the same time. The social networks are initially
constructed as weighted directed networks with edge weights either
representing the number of times that the focal individual is
co-captured with the other individual in a dyad (a count). Co-capture is
defined as being captured at the same (or a nearby) location within a
fixed period of time. This time period is set using the intwindow=
argument in the function and extends either side of each capture period
(i.e. intwindow=20 would result in interactions up to 20 days before and
after a capture event as co-captures). It is then possible to generate
an undirected network using this co-capture information by calculating
either the minimum, mean or maximum recorded edge weight for each dyad.
An example of how this is achieved is provided by Case Study 1. Once
constructed, social networks are outputted in both edge list and
adjacency matrix format alongside a dataframe recording the network
windows that each individual was captured in (binary values for captured
and not captured).

-----

### Movement network construction

CMRnet constructs movement networks using the movement of individuals
(the edges) to connect between different capture locations (the nodes)
(Fig. 1). Movement networks are constructed as weighted, directed
networks. Many options in the function are identical to those used to
construct social networks. One important difference is that the
interaction window only runs subsequent to the capture event (i.e. when
intwindow=20 an edge is only added to the movement network when an
individual is captured at a different location within 20 days of the
focal capture event). The key differences in using the function are that
there is no option for to set the spacewindow= argument (i.e. the
locations in the input dataset are used as the nodes in the final
network), and there is an additional nextonly= argument. When this is
set to be true (the default) then, if there are multiple captures of an
individual at different locations within the interaction window, an edge
is only drawn to the location of the next capture (and not to any
subsequent captures that are still within the interaction window). The
package also contains an additional function – multi.move.net.create() –
that can be used to construct multiplex movement networks using an
additional sixth column to the input dataset. This function constructs a
separate network for each unique value in the “layer” column for each
network window. For example, if this column contained information on the
sex of individuals then the function would output two separate networks
for each time period, while if it repeated the ID column then it would
output nt separate networks for each time period (where nt is the size
of the population in that time period). Taking this multiplex approach
may facilitate comparisons between the spatial networks of different
individuals, phenotypes or demographic classes, and so may be preferable
to constructing a single movement network for the whole population.

-----

### Network analysis tools

The package also contains functions that can be used to construct
permuted networks from CMR data. The function cmr.nodeswap() is a shell
for applying the function rmperm() from the package sna (Butts, 2014) to
outputs from CMRnet network construction functions. This can be used to
conduct simultaneous row/column permutations on the already constructed
adjacency matrices for each network window. The package also contains
functions to construct data-stream permutations of the network. The
function datastream.perm.soc() generates randomised networks using an
algorithm that swaps individuals between capture events (see (Bejder,
Fletcher, & Bräger, 1998; Farine, 2017)). After a period of burn-in
(defined using the n.burnin= argument) as the randomised dataset becomes
increasingly dissimilar to the observed dataset, the algorithm
calculates social networks after a set number of swaps (defined using
the n.swaps= argument ). The algorithm proceeds until the pre-defined
(n.rand=) number of permuted networks have been generated. The total
number of swaps therefore equals n.burnin+n.swaps×n.rand. Swaps can be
restricted so that they can only occur within in a particular time
window or geographic area, and examples of how this is coded are
provided in the supplementary material. Datastream permutations for
movement networks are conducted using the datastream.perm.spat()
function, and the algorithm used is very similar apart from that
locations rather than individual identities are swapped between capture
events. This makes it possible to constrain permutations so that swaps
only occur between different captures of the same individual, as well as
within particular time windows. Examples of this approach are provided
in the supplementary material.

-----

### References

  - Bejder, L., Fletcher, D., & Bräger, S. (1998). A method for testing
    association patterns of social animals. Animal Behaviour, 56(3),
    719–725.
  - Butts, C. T. (2014). Package “sna.”
  - Carter, A. J., Lee, A. E. G., & Marshall, H. H. (2015). Research
    questions should drive edge definitions in social network studies.
  - Farine, D. R. (2017). A guide to null models for animal social
    network analysis. Methods in Ecology and Evolution, 8(10),
    1309–1320.

\name{KrigInv-package}
\alias{KrigInv}
\docType{package}

\title{ Kriging-Based Inversion for Deterministic and Noisy Computer Experiments }

\description{
Kriging-based sequential algorithms, meant to identify the excursion set of a real valued function. The algorithms can also identify one or several level sets.
}

\details{
\tabular{ll}{
Package: \tab KrigInv\cr
Type: \tab Package\cr
Version: \tab 1.4.1\cr
Date: \tab 2018-09-04\cr
License: \tab GPL version 3\cr
LazyLoad: \tab yes\cr
}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)

Victor Picheny (INRA, Toulouse, France)

David Ginsbourger (IDIAP Martigny and University of Bern, Switzerland)

Dario Azzimonti (IDSIA, Switzerland)

with contributions from Yann Richet (IRSN, France)

Maintainer: Clement Chevalier (clement.chevalier@unine.ch)
}

\note{
Important functions are \code{\link{EGI}} and \code{\link{EGIparallel}}. The last 1.4 version allows to handle multiple thresholds T, stored in an array and implements conservative excursion set strategies.

A first prototype of this package was originally developed by D. Ginsbourger in the frame of a collaboration with IRSN (Institut de Radioprotection et de Surete Nucleaire), acting through Yann Richet. The three main authors thank IRSN for sponsoring open source research, and allowing them to spread the present package and publish it on CRAN. They also would like to warmly thank Yann Richet for numerous discussions concerning this package, and more!
}




\references{
Azzimonti, D., Ginsbourger, D., Chevalier, C., Bect, J., and Richet, Y. (2018). \emph{Adaptive design of experiments for conservative estimation of excursion sets}. Under revision. Preprint at \href{https://hal.archives-ouvertes.fr/hal-01379642}{hal-01379642}

Chevalier C., Picheny V., Ginsbourger D. (2014), \emph{Kriginv: An efficient and user-friendly implementation of batch sequential inversion strategies based on kriging} (2014) Computational Statistics & Data Analysis, vol. 71, pp 1021-1034

Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465

Chevalier C., Ginsbourger D. (2014), \emph{Corrected Kriging update formulae for batch-sequential data assimilation}, in Pardo-Iguzquiza, E., et al. (Eds.) Mathematics of Planet Earth, pp 119-122

Chevalier C. (2013) \emph{Fast uncertainty reduction strategies relying on Gaussian process models} Ph.D Thesis, University of Bern

Picheny V., Ginsbourger D., Roustant O., Haftka R.T., (2010) \emph{Adaptive designs of experiments for accurate approximation of a target region}, J. Mech. Des. vol. 132(7)

Picheny V. (2009) \emph{Improving accuracy and compensating for uncertainty in surrogate modeling}, Ph.D. thesis, University of Florida and Ecole Nationale Superieure des Mines de Saint-Etienne

Bichon B.J., Eldred M.S., Swiler L.P., Mahadevan S., McFarland J.M. (2008) \emph{Efficient global reliability analysis for nonlinear implicit performance functions}, AIAA Journal 46(10), pp 2459-2468

Ranjan P., Bingham D., Michailidis G. (2008) \emph{Sequential experiment design for contour estimation from complex computer codes} Technometrics 50(4), pp 527-541
}

\keyword{package}

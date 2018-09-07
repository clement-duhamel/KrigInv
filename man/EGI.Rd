\name{EGI}
\Rdversion{1.4}
\alias{EGI}

\title{
Efficient Global Inversion: sequential inversion algorithm based on Kriging. 
}
\description{
Sequential sampling based on the optimization of a kriging-based criterion, with model update after each iteration. The criterias aim at identifying an excursion set or one/many level sets. At each iteration batchsize = 1 new locations are evaluated. 
Different criteria are available for selecting experiments. The pointwise criteria are \code{"bichon"}, \code{"ranjan"}, \code{"tmse"}, \code{"tsee"}  and are fast to compute. In addition, integral criteria require numerical integration and can potentially deliver more than one new location per iteration. Available integral criteria are \code{"imse"}, \code{"timse"}, \code{"sur"}, \code{"jn"}, \code{"vorob"}. The use of the integral criteria is implemented in the \code{\link{EGIparallel}} function.
}
\usage{
EGI(T, model, method = NULL, method.param = NULL,
fun, iter, lower, upper, new.noise.var = 0,
optimcontrol = NULL, kmcontrol = NULL, integcontrol = NULL, ...)
}

\arguments{
  \item{T}{
Array containing one or several thresholds. The criteria which can be used with multiple thresholds are \code{"tmse"}, \code{"timse"}, \code{"sur"}, \code{"jn"}.
}
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{method}{
Criterion used for choosing observations. 
}
  \item{method.param}{
Optional tolerance value (scalar) for methods \code{"ranjan"}, \code{"bichon"}, \cr
\code{"tmse"} and \code{"timse"}. If not provided, default value is used (1 for ranjan and bichon, 0 for tmse and timse).
}
  \item{fun}{
Objective function.
}
  \item{iter}{
Number of iterations (i.e. number of additional sampling points).
}
  \item{lower}{
Vector containing the lower bounds of the design space.
}
  \item{upper}{
Vector containing the upper bounds of the design space.
}
  \item{new.noise.var}{
Optional scalar value of the noise variance of the new observations.
}
  \item{optimcontrol}{
Optional list of control parameters for the optimization of the sampling criterion. The field \code{method} defines which optimization method is used: it can be either \code{"genoud"} (default) for an optimisation using the genoud algorithm, or \code{"discrete"} for an optimisation over a specified discrete set. 
If the field \code{method} is set to \code{"genoud"}, one can set some parameters of this algorithm: 
\code{pop.size}  (default : 50*d),  \code{max.generations} (default : 10*d), \cr
\code{wait.generations} (2),  \code{BFGSburnin} (2) and the mutations \code{P1}, \code{P2}, up to \code{P9} (see \code{\link[rgenoud]{genoud}}). Numbers into brackets are the default values.
If the field \code{method} is set to \code{"discrete"}, one can set the field \code{optim.points}: p * d matrix corresponding to the p points where the criterion will be evaluated. If nothing is specified, 100*d points are chosen randomly. 
}
  \item{kmcontrol}{
Optional list representing the control variables for the re-estimation of the kriging model once new points are sampled. 
The items are the same as in \code{\link[DiceKriging]{km}}.
}
  \item{integcontrol}{
Optional list specifying the procedure to build the integration points and weights, relevant only for the sampling criteria based on numerical integration: \cr
(\code{"imse"}, \code{"timse"}, \code{"sur"} or \code{"jn"}). Many options are possible.
A) If nothing is specified, 100*d points are chosen using the Sobol sequence. 
B) One can directly set the field \code{integration.points} (a p * d matrix) for prespecified integration points. In this case these integration points and the corresponding vector \code{integration.weights} will be used for all the iterations of the algorithm. 
C) If the field \code{integration.points} is not set then the integration points are renewed at each iteration. 
In that case one can control the number of integration points \code{n.points} (default: 100*d) and a specific distribution \code{distrib}. Possible values for \code{distrib} are: 
\code{"sobol"}, \code{"MC"}, \code{"timse"}, \code{"imse"}, \code{"sur"} and \code{"jn"} (default: \code{"sobol"}). 
C.1) The choice \code{"sobol"} corresponds to integration points chosen with the Sobol sequence in dimension d (uniform weight). 
C.2) The choice \code{"MC"} corresponds to points chosen randomly, uniformly on the domain. 
C.3) The choices \code{"timse"}, \code{"imse"}, \code{"sur"} and \code{"jn"} correspond to importance sampling distributions (unequal weights). It is strongly recommended to use the importance sampling distribution corresponding to the chosen sampling criterion. 
When important sampling procedures are chosen, \code{n.points} points are chosen using importance sampling among a discrete set of \code{n.candidates} points (default: \code{n.points*10}) which are distributed according to a distribution \cr
 \code{init.distrib} (default: \code{"sobol"}). Possible values for \code{init.distrib} are the space filling distributions \code{"sobol"} and \code{"MC"} or an user defined distribution \code{"spec"}. 
The \code{"sobol"} and \code{"MC"} choices correspond to quasi random and random points in the domain. 
If the \code{"spec"} value is chosen the user must fill in manually the field \code{init.distrib.spec} to specify himself a n.candidates * d matrix of points in dimension d.
}
  \item{\dots}{
Other arguments of the target function \code{fun}.
}
}
\details{
The function used to build the integration points and weights (based on the options specified in \code{integcontrol}) is the function \code{\link{integration_design}}
}
\value{A list with components:
\item{par}{The added observations (ite * d matrix)}
\item{value}{The value of the function \code{fun} at the added observations (vector of size "ite")}
\item{nsteps}{The number of added observations (=ite).}
\item{lastmodel}{The current (last) kriging model of \code{\link[DiceKriging]{km}} class.}
\item{lastvalue}{The value of the criterion at the last added point.}
\item{allvalues}{If an optimization on a discrete set of points is chosen, the value of the criterion at all these points, for the last iteration.}
}
\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465

Picheny V., Ginsbourger D., Roustant O., Haftka R.T., (2010) \emph{Adaptive designs of experiments for accurate approximation of a target region}, J. Mech. Des. vol. 132(7)

Chevalier C. (2013) \emph{Fast uncertainty reduction strategies relying on Gaussian process models} Ph.D Thesis, University of Bern

Bichon B.J., Eldred M.S., Swiler L.P., Mahadevan S., McFarland J.M. (2008) \emph{Efficient global reliability analysis for nonlinear implicit performance functions}, AIAA Journal 46(10), pp 2459-2468

Ranjan P., Bingham D., Michailidis G. (2008) \emph{Sequential experiment design for contour estimation from complex computer codes} Technometrics 50(4), pp 527-541
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)

Victor Picheny (INRA, Toulouse, France)

David Ginsbourger (IDIAP Martigny and University of Bern, Switzerland)
}

\seealso{\code{\link{EGIparallel}},\code{\link{max_infill_criterion}}}

\examples{
#EGI

set.seed(9)
N <- 20 #number of observations
T <- 80 #threshold
testfun <- branin
lower <- c(0,0)
upper <- c(1,1)

#a 20 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

optimcontrol <- list(method="genoud",pop.size=50)
integcontrol <- list(distrib="sur",n.points=50)
iter <- 1

\dontrun{
obj1 <- EGI(T=T,model=model,method="sur",fun=testfun,iter=iter,
           lower=lower,upper=upper,optimcontrol=optimcontrol,
           integcontrol=integcontrol)

obj2 <- EGI(T=T,model=model,method="ranjan",fun=testfun,iter=iter,
           lower=lower,upper=upper,optimcontrol=optimcontrol)


par(mfrow=c(1,3))
print_uncertainty_2d(model=model,T=T,main="probability of excursion",
type="pn",new.points=0,cex.points=2)

print_uncertainty_2d(model=obj1$lastmodel,T=T,
main="updated probability of excursion, sur sampling",
type="pn",new.points=iter,col.points.end="red",cex.points=2)

print_uncertainty_2d(model=obj2$lastmodel,T=T,
main="updated probability of excursion, ranjan sampling",
type="pn",new.points=iter,col.points.end="red",cex.points=2)
}
##############
#same example with noisy initial observations and noisy new observations
branin.noise <- function(x) return(branin(x)+rnorm(n=1,sd=30))

set.seed(9)
N <- 20;T <- 80
testfun <- branin.noise
lower <- c(0,0);upper <- c(1,1)

design <- data.frame( matrix(runif(2*N),ncol=2) )
response.noise <- apply(design,1,testfun)
response.noise - response


model.noise <- km(formula=~., design = design, response = response.noise,
covtype="matern3_2",noise.var=rep(30*30,times=N))

optimcontrol <- list(method="genoud",pop.size=50)
integcontrol <- list(distrib="sur",n.points=50)
iter <- 1

\dontrun{
obj1 <- EGI(T=T,model=model.noise,method="sur",fun=testfun,iter=iter,
           lower=lower,upper=upper,optimcontrol=optimcontrol,
           integcontrol=integcontrol,new.noise.var=30*30)

obj2 <- EGI(T=T,model=model.noise,method="ranjan",fun=testfun,iter=iter,
           lower=lower,upper=upper,optimcontrol=optimcontrol,
            new.noise.var=30*30)

par(mfrow=c(1,3))
print_uncertainty_2d(model=model.noise,T=T,
main="probability of excursion, noisy obs.",
type="pn",new.points=0,cex.points=2)

print_uncertainty_2d(model=obj1$lastmodel,T=T,
main="probability of excursion, sur sampling, noisy obs.",
type="pn",new.points=iter,col.points.end="red",cex.points=2)

print_uncertainty_2d(model=obj2$lastmodel,T=T,
main="probability of excursion, ranjan sampling, noisy obs.",
type="pn",new.points=iter,col.points.end="red",cex.points=2)
}
}

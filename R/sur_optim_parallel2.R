sur_optim_parallel2 <- function(x,other.points, integration.points,integration.weights=NULL,
		                              intpoints.oldmean,intpoints.oldsd,precalc.data,
		                              model, T, new.noise.var=NULL,batchsize,current.sur,ai_precalc=NULL){
	
	x.complete <- c(x,other.points)
	return(sur_optim_parallel(
		x = x.complete, integration.points = integration.points, integration.weights = integration.weights,
		intpoints.oldmean = intpoints.oldmean,intpoints.oldsd = intpoints.oldsd,precalc.data = precalc.data,
		model = model,T = T,new.noise.var = new.noise.var,batchsize=batchsize,current.sur=current.sur,
		ai_precalc=ai_precalc
		)
	)
	
}

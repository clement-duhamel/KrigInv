tsee_optim <-
function(x, model, T){
	
	y <- t(x)
  if((nrow(y)==1) && (ncol(y)==model@d)){
    z <- y #only one point
  }else{
    #many points
    if(ncol(x) == model@d) z <- x
    if(ncol(x) != model@d) z <- y
  }

	krig <- predict_nobias_km(object=model,newdata=as.data.frame(z),type="UK",se.compute=TRUE)
	mk <- krig$mean; sk <- krig$sd
	t <- (T-mk)/sk; ski_dnorm_t <- sk * dnorm(t)
	C <-  ((T-mk) * pnorm(t) + ski_dnorm_t) * ((mk-T) * pnorm(-t) + ski_dnorm_t)
	C[is.nan(C)] <- 0
	return(C)
}


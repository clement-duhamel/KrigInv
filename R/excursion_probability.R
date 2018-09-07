excursion_probability <- function(mn,sn,T){
  
  if(length(T)==1){
    # one thresholds
    pn <- pnorm( (mn - T)/sn )
  }else{
    # multiple thresholds
    T <- sort(T)
    k <- length(T)
    
    pn <- rep(0,times=length(mn))
    for(i in 1:k){
      Ti <- T[i]
      pn <- pn + (-1)^(i+1) * pnorm( (mn - Ti)/sn )
    }
  }
  return(pn)
  
}
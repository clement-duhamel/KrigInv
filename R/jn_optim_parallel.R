jn_optim_parallel <- function(x, integration.points,integration.weights=NULL,
                               intpoints.oldmean,intpoints.oldsd,precalc.data,
                               model, T, new.noise.var=NULL,batchsize,current.sur,ai_precalc){
  
  if(!is.null(new.noise.var)){
    if(new.noise.var == 0) {
      new.noise.var <- NULL
    }
  }
  
  #x is a vector of size d * batchsize
  d <- model@d
  n <- model@n
  X.new <- matrix(x,nrow=d)
  mindist <- Inf
  
  tp1 <- c(as.numeric(t(model@X)),x)
  for (i in 1:batchsize){
    #distance between the i^th point and all other points (in the DOE or in the batch)
    xx <- X.new[,i]
    tp2<-matrix(tp1-as.numeric(xx),ncol=d,byrow=TRUE)^2
    mysums <- sqrt(rowSums(tp2))
    mysums[n+i] <- Inf #because this one is always equal to zero
    mindist <- min(mindist,mysums)		
  }
  
  if (!identical(colnames(integration.points), colnames(model@X))) colnames(integration.points) <- colnames(model@X)
  
  if ((mindist > 1e-5) || (!is.null(new.noise.var))){
    X.new <- t(X.new)	
    krig  <- predict_nobias_km(object=model, newdata=as.data.frame(X.new), 
                type="UK",se.compute=TRUE, cov.compute=TRUE) 
    
    mk <- krig$mean ; sk <- krig$sd ; newXvar <- sk*sk
    F.newdata <- krig$F.newdata ; c.newdata <- krig$c;Sigma.r <- krig$cov
    
    kn = computeQuickKrigcov(model,integration.points,X.new,precalc.data, F.newdata , c.newdata) 
    
    krig2  <- predict_update_km_parallel (newXmean=mk,newXvar=newXvar,newXvalue=mk, 
                Sigma.r=Sigma.r,newdata.oldmean=intpoints.oldmean,newdata.oldsd=intpoints.oldsd,kn=kn)
    if(!is.null(krig2$error)) return(current.sur)
    
    sk.new <- krig2$sd	
    c <- (intpoints.oldsd*intpoints.oldsd)/(sk.new*sk.new)
    c[c==Inf]<- 1000; c[is.nan(c)] <- 1000
    ai_precalc <- t(ai_precalc)/sk.new
    ai_precalc[ai_precalc==Inf]<- 1000 ;ai_precalc[ai_precalc== -Inf] <- -1000;ai_precalc[is.nan(ai_precalc)] <- 1000
    ai_precalc <- ai_precalc/sqrt(c) #a.new
    b <- krig2$lambda / sk.new
    integration.case <- NULL
    tmp <- 0
    nT <- length(T)
    
    if(length(integration.points) == (model@d * 2 * length(integration.weights))){
      # integration points generated with the "jn" distribution
      integration.case <- 1
      M <- nrow(integration.points) / 2
        
      a1 <- ai_precalc[1:M,]
      a2 <- ai_precalc[(M+1):(2*M),]
      c1 <- c[1:M]
      c2 <- c[(M+1):(2*M)]
      bz1 <- b[(1:M),]
      bz2 <- b[(M+1):(2*M),]
    }else{
      #general case with M^2 integration points 
      integration.case <- 2
      M <- nrow(integration.points)
      indices <- expand.grid(c(1:M),c(1:M))
      col1 <- indices[,1] ; col2 <- indices[,2]
      
      a1 <- ai_precalc[col1,]
      a2 <- ai_precalc[col2,]
      c1 <- c[col1]
      c2 <- c[col2]
      bz1 <- b[col1,]
      bz2 <- b[col2,]
    }
    
    if(nT == 1){a1 <- matrix(a1,ncol=1);a2 <- matrix(a2,ncol=1)}
    if(batchsize==1) bz2 <- matrix(bz2,ncol=1)
    d_z1_z2 <- colSums(t(bz1)*tcrossprod(Sigma.r,bz2))
    arg3 <- d_z1_z2 / sqrt(c1*c2)
    arg3 <- pmin(arg3,1);arg3 <- pmax(arg3,-1)
    
    if(any(is.na(arg3))){return(current.sur)}
    
    for(i in 1:nT){
      ai.new <- a1[,i]
      for(j in 1:nT){
        aj.new <- a2[,j]
        tmp <- tmp + (-1)^(i+j)*pbivnorm(ai.new,aj.new,arg3)
      }
    }
    
    if(integration.case == 1){
      crit <- - sum(tmp*integration.weights)
    }else{
      if (is.null(integration.weights)) {
        crit <- - mean(tmp)
      }else{ 
        new.weights <- integration.weights[col1]*integration.weights[col2]
        crit <- - sum(tmp*new.weights)
      }
    }
  }else crit <- current.sur		
  return(crit)
}
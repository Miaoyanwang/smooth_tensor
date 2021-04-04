
Objective = function (x, mu.array, Cs, Ds, Es, lambda = 0, method="L0") {
  if(method=="L0") return(sum((x - mu.array[Cs, Ds, Es, drop=FALSE])^2,na.rm = T)+2*lambda*sum(mu.array !=0 ,na.rm = T))
  if(method=="L1") return(sum((x - mu.array[Cs, Ds, Es, drop=FALSE])^2,na.rm = T)+2*lambda*sum(abs(mu.array),na.rm = T))
  stop("No such kind of method:", method,".\n")
}

makediagNA = function(n){
  
  tnsr = array(1,dim = c(n,n,n))
  tnsr = array( apply( tnsr, 3, function(x) {x[ row(x) == col(x) ] <-NA; x} ), dim(tnsr) )
  tnsr = array( apply( tnsr, 2, function(x) {x[ row(x) == col(x) ] <-NA; x} ), dim(tnsr) )
  tnsr = array( apply( tnsr, 1, function(x) {x[ row(x) == col(x) ] <-NA; x} ), dim(tnsr) )
  
  
  return(tnsr)
}


Soft = function (a, b, method="L0"){
  if (b < 0) 
    stop("Can soft-threshold by a nonnegative quantity only.")
  if(method == "L0") return(sign(a) * ifelse(abs(a)>b, abs(a), 0))
  if(method == "L1") return(sign(a) * pmax(0, abs(a) - b))
  stop("No such kind of method:", method,".\n")
}


tensor_unfold = function(tensor,dim=1){
  if (dim == 1) unfold = aperm(tensor,c(3,2,1))
  if (dim == 2) unfold = aperm(tensor,c(1,3,2))
  if (dim == 3) unfold = tensor
  unfold = apply(unfold,3,c)
  if (is.vector(unfold)) return(as.matrix(unfold)) else return(t(unfold))
}


ReNumber = function (Cs) 
{
  newCs <- rep(NA, length(Cs))
  uniq <- sort(unique(Cs))
  for (i in 1:length(uniq)) {
    newCs[Cs == uniq[i]] <- i
  }
  return(newCs)
}

## cut function
cut = function(tnsr){
  cond = (dim(tnsr)[1]==dim(tnsr)[2])&(dim(tnsr)[1]==dim(tnsr)[3])
  if (!(cond)){
    stop("* tensor must have same dimension for each mode")
  }
  # n = dim(tnsr)[1]
  # for(i in 1:n){
  #   for(j in 1:n){
  #     for(k in 1:n)
  #       if((i-j)*(j-k)*(k-i)==0){
  #         tnsr[i,j,k] = 0
  #       }
  #   }
  # }
  tnsr = array( apply( tnsr, 3, function(x) {x[ row(x) == col(x) ] <-0; x} ), dim(tnsr) )
  tnsr = array( apply( tnsr, 2, function(x) {x[ row(x) == col(x) ] <-0; x} ), dim(tnsr) )
  tnsr = array( apply( tnsr, 1, function(x) {x[ row(x) == col(x) ] <-0; x} ), dim(tnsr) )
  
  
  return(tnsr)
}






#give x as an array
UpdateMus_tensor = function (x, Cs, Ds, Es, lambda=0, method="L0") {
  uniqCs = sort(unique(Cs))
  uniqDs = sort(unique(Ds))
  uniqEs = sort(unique(Es))

  mus = array(NA, c(length(uniqCs), length(uniqDs), length(uniqEs)))
  if(method=="L1"){
    for (k in uniqCs){
      for (r in uniqDs){
        for (l in uniqEs){
          if (lambda == 0) mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l],na.rm = T)
          if (lambda > 0) mus[k,r,l] = Soft(mean(x[Cs==k,Ds==r,Es==l],na.rm = T),lambda/(2*sum(Cs==k)*sum(Ds==r)*sum(Es==l)),method=method)
          if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
        }
      }
    }
  }## added 
  if(method=="L0"){
    for (k in uniqCs){
      for (r in uniqDs){
        for (l in uniqEs){
          if (lambda == 0) mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l],na.rm = T)
          if (lambda > 0) mus[k,r,l] = Soft(mean(x[Cs==k,Ds==r,Es==l],na.rm = T),sqrt(lambda/(sum(Cs==k)*sum(Ds==r)*sum(Es==l))),method=method)
          ### modified
          ## Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/sum(Cs==k)*sum(Ds==r)*sum(Es==l),method=method)
          if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
        }
      }
    }
  }## added
  return(mus)
}

UpdateClusters_tensor = function (x, mus, curCs, curDs) {
  Cs.new <- rep(NA, length(curCs))
  uniq <- 1:max(curCs)
  mus.expandcolumns <- mus[, curDs, drop = FALSE]
  #if (dim(mus.expandcolumns)[1]==1) mus
  for (i in 1:nrow(x)) {
    dist2.clust <- NULL
    for (k in 1:length(uniq)) {
      #see which cluster is closest to one sample
      dist2.clust <- c(dist2.clust, sum((x[i, , drop = FALSE] - 
                                           mus.expandcolumns[k, , drop = FALSE])^2,na.rm = T))
    }
    wh <- which(dist2.clust == min(dist2.clust))
    Cs.new[i] <- wh[1]
  }
  return(Cs.new)
}




classify2 = function(x,k,r,l,sym = F,diagP = T,lambda=0,max.iter=1000,threshold = 1e-15,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,nstart=80,method="L0",center=FALSE){
  n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]
  if(center == TRUE) {
    mustemp <- mean(x)
    x <- x-mustemp
  }
  if(is.null(Cs.init)){
    if(k==1) Cs = rep(1,n) else {Cs  = kmeans(tensor_unfold(x,1),k,nstart = nstart)$cluster}
  } else {
    Cs = ReNumber(Cs.init)
    k = length(unique(Cs))
  }
  if(is.null(Ds.init)){
    if(r==1) Ds = rep(1,p) else {Ds  = kmeans(tensor_unfold(x,2),r,nstart = nstart)$cluster}
  } else {
    Ds = ReNumber(Ds.init)
    r = length(unique(Ds))
  }
  if(is.null(Es.init)){
    if(l==1) Es = rep(1,q) else {Es  = kmeans(tensor_unfold(x,3),l,nstart = nstart)$cluster}
  } else {
    Es = ReNumber(Es.init)
    l = length(unique(Es))
  }
  if(diagP==F){
    x = x*makediagNA(n)
  }
  objs <- 1e+15
  improvement <- 1e+10
  i <- 1
  mu.array = UpdateMus_tensor(x,Cs,Ds,Es,lambda,method=method)
  while((improvement > threshold) & (i <= max.iter)){
    Cs = UpdateClusters_tensor(tensor_unfold(x),tensor_unfold(mu.array),Cs,(rep(Ds,each=q)-1)*l+rep(Es,times=p))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Cs <- ReNumber(Cs)
    k = length(unique(Cs))
    mu.array = UpdateMus_tensor(x,Cs,Ds,Es,lambda,method=method)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Ds = UpdateClusters_tensor(tensor_unfold(x,2),tensor_unfold(mu.array,2),Ds,(rep(Es,each=n)-1)*k+rep(Cs,times=q))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Ds <- ReNumber(Ds)
    r = length(unique(Ds))
    mu.array = UpdateMus_tensor(x,Cs,Ds,Es,lambda,method=method)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Es = UpdateClusters_tensor(tensor_unfold(x,3),tensor_unfold(mu.array,3),Es,(rep(Ds,each=n)-1)*k+rep(Cs,times=p))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Es <- ReNumber(Es)
    l = length(unique(Es))
    mu.array = UpdateMus_tensor(x,Cs,Ds,Es,lambda, method=method)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    improvement <- abs(objs[length(objs)] - objs[length(objs) - 
                                                   6])/abs(objs[length(objs) - 6])
    i <- i + 1
    if(trace) cat("step",i,",improvement=",improvement,".\n")
    if (is.na(improvement)) break
  }
  if (sym ==T){
    symobj = NULL
    
    temp = UpdateMus_tensor(x,Cs,Cs,Cs,lambda, method=method) 
    symobj = c(symobj,Objective(x, temp, Cs, Cs, Cs, lambda = lambda, method=method))
    temp = UpdateMus_tensor(x,Ds,Ds,Ds,lambda, method=method)
    symobj = c(symobj,Objective(x, temp, Ds, Ds, Ds, lambda = lambda, method=method))
    temp = UpdateMus_tensor(x,Es,Es,Es,lambda, method=method)
    symobj = c(symobj,Objective(x, temp, Es, Es, Es, lambda = lambda, method=method))
    Cs = Ds = Es = list(Cs,Ds,Es)[[which.min(symobj)]]
    mu.array = UpdateMus_tensor(x,Cs,Ds,Es,lambda, method=method)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
  }
  if (i > max.iter) {
    warning("The algorithm has not converged by the specified maximum number of iteration.\n")
  }
  if(center==TRUE){
    mu.array <- mu.array + mustemp
  }
  mu.array[abs(mu.array)<=1e-6] = 0
  return(list("judgeX"=mu.array[Cs,Ds,Es, drop=FALSE],"Cs"=Cs,"Ds"=Ds,"Es"=Es,"objs"=objs[length(objs)], "mus"=mu.array))
}


tbmClustering = function(x,k,r,l,lambda=0,sym = F,diagP = T,max.iter=100,threshold = 1e-3,sim.times=1,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,method="L0"){
  #x=test;lambda=1e-3;max.iter=200;threshold = 5e-3;sim.times=10
  if (sim.times == 1){
    result = classify2(x,k,r,l,sym,diagP,lambda=lambda,max.iter = max.iter,threshold = threshold,Cs.init = Cs.init,Ds.init = Ds.init,Es.init = Es.init,method=method)
    if(diagP == F){
      result$judgeX = cut(result$judgeX)
    }
    return(result)
  } 
  if (.Platform$OS.type == "windows") {
    result = lapply(rep(list(x),sim.times), classify2, k,r,l,sym,diagP,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init,method=method)
    objs = unlist(lapply(result, function(result){result$objs}))
  } else {
    result = mclapply(rep(list(x),sim.times), classify2, k,r,l,sym,diagP,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init,nstart = sample(1:1000,1),method=method,mc.cores = n.cores)
    objs = unlist(lapply(result, function(result){result$objs}))
  }
  result = result[[which(objs == min(objs))[1]]]
  if(diagP == F){
    result$judgeX = cut(result$judgeX)
  }
  return(result)
}




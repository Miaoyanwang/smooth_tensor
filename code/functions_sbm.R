#library("parallel")
#library("rTensor")
ReNumber2 = function(Cs){
  #Cs=c(2,2)
  newCs <- rep(NA, length(Cs))
  uniq <- sort(unique(Cs))
  num = c()
  for (i in 1:length(uniq)) {
    num[i] = sum(Cs==i)
  }
  newCs[order(Cs)] = rep(order(num),num)
  return(newCs)
}

reorderClusters = function(x,Cs,Ds,Es){
  #x=truthX;Cs=truthCs;Ds=truthDs;Es=truthEs
  Cs = ReNumber2(Cs);Ds = ReNumber2(Ds);Es = ReNumber2(Es)
  return(list("x"=x[order(Cs),order(Ds),order(Es), drop=FALSE],"Cs"=sort(Cs),
              "Ds"=sort(Ds),"Es"=sort(Es)))
}

tensor_unfold = function(tensor,dim=1){
  if (dim == 1) unfold = aperm(tensor,c(3,2,1))
  if (dim == 2) unfold = aperm(tensor,c(1,3,2))
  if (dim == 3) unfold = tensor
  unfold = apply(unfold,3,c)
  if (is.vector(unfold)) return(as.matrix(unfold)) else return(t(unfold))
}

design_row = function(sp,C1,D1,E1,k,r,l){
  labels = array(rep(0,k*r*l),c(k,r,l))
  #print(C1[sp[1]])
  labels[C1[sp[1]],D1[sp[2]],E1[sp[3]]] = 1
  return(c(labels))
}

Objective = function (x, mu.array, Cs, Ds, Es, lambda = 0, method="L0") {
  if(method=="L0") return(sum((x - mu.array[Cs, Ds, Es, drop=FALSE])^2)+2*lambda*sum(mu.array !=0 ))
  if(method=="L1") return(sum((x - mu.array[Cs, Ds, Es, drop=FALSE])^2)+2*lambda*sum(abs(mu.array)))
  stop("No such kind of method:", method,".\n")
}

#dim should be a vector
tensor_index = function(index,dims){
  index = index-1
  Cs = (index %% dims[1])+1
  Ds = (index %% (dims[1]*dims[2]))%/%dims[1] +1
  Es = (index %/% (dims[1]*dims[2]))+1
  return(c(Cs,Ds,Es))
}


Soft = function (a, b, method="L0"){
  if (b < 0) 
    stop("Can soft-threshold by a nonnegative quantity only.")
  if(method == "L0") return(sign(a) * ifelse(abs(a)>b, abs(a), 0))
  if(method == "L1") return(sign(a) * pmax(0, abs(a) - b))
  stop("No such kind of method:", method,".\n")
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
        if (lambda == 0) mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l])
        if (lambda > 0) mus[k,r,l] = Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/(2*sum(Cs==k)*sum(Ds==r)*sum(Es==l)),method=method)
        if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
      }
    }
  }
  }## added 
  if(method=="L0"){
      for (k in uniqCs){
          for (r in uniqDs){
              for (l in uniqEs){
                  if (lambda == 0) mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l])
                  if (lambda > 0) mus[k,r,l] = Soft(mean(x[Cs==k,Ds==r,Es==l]),sqrt(lambda/(sum(Cs==k)*sum(Ds==r)*sum(Es==l))),method=method)
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
                                           mus.expandcolumns[k, , drop = FALSE])^2))
    }
    wh <- which(dist2.clust == min(dist2.clust))
    Cs.new[i] <- wh[1]
  }
  return(Cs.new)
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

label_for_krl = function(abc,k,r,l,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,sim.times,lambda,xmiss,method="L0",crossvalidation = TRUE){
  a = abc[1]; b = abc[2]; c = abc[3]
  if(is.null(Cs.init) == FALSE) Cs.init = Cs.init[,a]
  if(is.null(Ds.init) == FALSE) Ds.init = Ds.init[,b]
  if(is.null(Es.init) == FALSE) Es.init = Es.init[,c]
  if (crossvalidation == TRUE){
    return(list(tbmClustering(xmiss, k[a], r[b], l[c], lambda = lambda, 
                       Cs.init = Cs.init, Ds.init = Ds.init, 
                       Es.init = Es.init, sim.times=sim.times, method=method)$judgeX))} else {
                         return(tbmClustering(xmiss, k[a], r[b], l[c], lambda = lambda, 
                                       Cs.init = Cs.init, Ds.init = Ds.init, 
                                       Es.init = Es.init, sim.times=sim.times, method=method))
                       }
}

positionfun=function(d){
  x=rep(1,d[2])%x%rep(1,d[3])%x%c(1:d[1])
  y=rep(1,d[3])%x%(c(1:d[2])%x%rep(1,d[1]))
  z=c(1:d[3])%x%rep(1,d[1])%x%rep(1,d[2])
  position=cbind(x,y,z)
  return(list("position"=position))
}


mse = function(bires, data){
    npq = dim(bires$judgeX)
    n = npq[1]; p = npq[2]; q = npq[3]
    return(sum((bires$judgeX-data$truthX)^2)/n/p/q)
} 



# label_for_cp = function(multiplicative=1,x,k,r,l){
#   cp_result = cp(as.tensor(x),num_components = multiplicative)
#   lambda = cp_result$lambdas
#   fitted=attributes(cp_result$est)$data
#   
# 
#   n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]
#   mus = array(rep(0,n*p*q),c(n,p,q))
#   Cs = kmeans(cp_result$U[[1]],k)$cluster
#   Ds = kmeans(cp_result$U[[2]],r)$cluster
#   Es = kmeans(cp_result$U[[3]],l)$cluster
#   
#   for(i in 1:k){
#       for(j in 1:r){
#           for(k in 1:l){
#               mus[Cs==i,Ds==j,Es==k]=mean(x[Cs==1,Ds==j,Es==k])
#           }
#           }
#       }
#   return(list(judgeX=fitted,s=multiplicative,Cs=Cs,Ds=Ds,Es=Es,blockmean=mus,mus=fitted))
# }


cluster2block=function(mu,Cs,Ds,Es){
    d=c(length(Cs),length(Ds),length(Es))
    r=c(length(unique(Cs)),length(unique(Ds)),length(unique(Es)))
    block=array(0,dim=d)
    for(i in 1:r[1]){
        for(j in 1:r[2]){
            for(k in 1:r[3]){
                if(mu[i,j,k]==0){
                    block[which(Cs==i),which(Ds==j),which(Es==k)]="0 0 0"
                }
                else{
               block[which(Cs==i),which(Ds==j),which(Es==k)]=paste(unique(Cs)[i],unique(Ds)[j],unique(Es)[k]) 
                }
            }
        }
    }
    return(block)
    }


tensor_unfold4 = function(tensor,dim=1){
  return(t(apply(tensor,dim,c)))
}


Objective4 = function (x, mu.array, Cs, Ds, Es, Fs, lambda = 0, method="L0") {
  if(method=="L0") return(sum((x - mu.array[Cs, Ds, Es, Fs, drop=FALSE])^2)+2*lambda*sum(mu.array !=0 ))
  if(method=="L1") return(sum((x - mu.array[Cs, Ds, Es, Fs, drop=FALSE])^2)+2*lambda*sum(abs(mu.array)))
  stop("No such kind of method:", method,".\n")
}

UpdateMus_tensor4 = function (x, Cs, Ds, Es, Fs, lambda=0, method="L0") {
  uniqCs = sort(unique(Cs))
  uniqDs = sort(unique(Ds))
  uniqEs = sort(unique(Es))
  uniqFs = sort(unique(Fs))
  mus = array(NA, c(length(uniqCs), length(uniqDs), length(uniqEs), length(uniqFs)))
  if(method=="L1"){
    for (k in uniqCs){
      for (r in uniqDs){
        for (l in uniqEs){
          for(m in uniqFs){
            if (lambda == 0) mus[k,r,l,m] = mean(x[Cs==k,Ds==r,Es==l,Fs=m])
            if (lambda > 0) mus[k,r,l,m] = Soft(mean(x[Cs==k,Ds==r,Es==l,Fs=m]),lambda/(2*sum(Cs==k)*sum(Ds==r)*sum(Es==l)*sum(Fs==m)),method=method)
            if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
          }
        }
      }
    }
  }## added 
  if(method=="L0"){
    for (k in uniqCs){
      for (r in uniqDs){
        for (l in uniqEs){
          for (m in uniqFs) {
            if (lambda == 0) mus[k,r,l,m] = mean(x[Cs==k,Ds==r,Es==l,Fs==m])
            if (lambda > 0) mus[k,r,l,m] = Soft(mean(x[Cs==k,Ds==r,Es==l,Fs==m]),sqrt(lambda/(sum(Cs==k)*sum(Ds==r)*sum(Es==l)*sum(Fs==m))),method=method)
            ### modified
            ## Soft(mean(x[Cs==k,Ds==r,Es==l]),lambda/sum(Cs==k)*sum(Ds==r)*sum(Es==l),method=method)
            if (lambda < 0) stop("Cannot have a negative tuning parameter value.")
          }
        }
      }
    }
  }## added
  return(mus)
}

# Perform tensor clustering via TBM method
# 
# This function performs sparse clustering on a three-dimensional tensor via TBM method.
# @param x a three-dimensional array
# @param k \eqn{d_1}: the clusters number of mode 1
# @param r \eqn{d_2}: the clusters number of mode 2
# @param l \eqn{d_3}: the clusters number of mode 3
# @param lambda a positive numeric value. The coefficient of the regularized term.
# @param max.iter a positive integer. The Maximum times of iteration.
# @param threshold a positive small numeric value which determines whether the algorithm converges or not.
# @param trace logic value. If true, it would print the iteration situation.
# @param Cs.init vector or NULL. Initial clsuter result of mode 1.
# @param Ds.init vector or NULL. Initial clsuter result of mode 2.
# @param Es.init vector or NULL. Initial clsuter result of mode 3.
# @param nstart positive interger. The same as the "nstart" in kmeans().
# @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
# @param center logic value that indicates whether run "x = x-mean(x)" before performing clustering.
# @return a list   
# \code{judgeX} estimated underlying mean signal.   
# 
#                \code{Cs} clustering result of mode 1.  
#                
#                \code{Ds} clustering result of mode 2.  
#                
#                \code{Es} clustering result of mode 3.  
#                
#                \code{mus} estimated underlying mean signal of each cluster.  
# 

classify2 = function(x,k,r,l,sym = F,lambda=0,max.iter=1000,threshold = 1e-15,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,nstart=20,method="L0",center=FALSE){
  n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]
  if(center == TRUE) {
    mustemp <- mean(x)
    x <- x-mustemp
  }
  if(is.null(Cs.init)){
    if(k==1) Cs = rep(1,n) else {Cs  = kmeans(tensor_unfold(x,1),k,nstart = nstart)$cluster}
  } else {
    Cs = Cs.init
  }
  if(is.null(Ds.init)){
    if(r==1) Ds = rep(1,p) else {Ds  = kmeans(tensor_unfold(x,2),r,nstart = nstart)$cluster}
  } else {
    Ds = Ds.init
  }
  if(is.null(Es.init)){
    if(l==1) Es = rep(1,q) else {Es  = kmeans(tensor_unfold(x,3),l,nstart = nstart)$cluster}
  } else {
    Es = Es.init
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
    #mu.array = UpdateMus_tensor(x,Cs,Ds,Es,lambda,method=method)
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Ds = UpdateClusters_tensor(tensor_unfold(x,2),tensor_unfold(mu.array,2),Ds,(rep(Es,each=n)-1)*k+rep(Cs,times=q))
    objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es, lambda = lambda, method=method))
    Ds <- ReNumber(Ds)
    r = length(unique(Ds))
    #mu.array = UpdateMus_tensor(x,Cs,Ds,Es,lambda,method=method)
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
    Cs = Ds = ES = list(Cs,Ds,Es)[[which.min(symobj)]]
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

# Perform tensor clustering via TBM method
# 
# This function performs sparse clustering on a three-dimensional tensor via TBM method.
# @param x a four-dimensional array
# @param k \eqn{d_1}: the clusters number of mode 1
# @param r \eqn{d_2}: the clusters number of mode 2
# @param l \eqn{d_3}: the clusters number of mode 3
# @param m \eqn{d_4}: the clusters number of mode 4
# @param lambda a positive numeric value. The coefficient of the regularized term.
# @param max.iter a positive integer. The Maximum times of iteration.
# @param threshold a positive small numeric value which determines whether the algorithm converges or not.
# @param trace logic value. If true, it would print the iteration situation.
# @param Cs.init vector or NULL. Initial clsuter result of mode 1.
# @param Ds.init vector or NULL. Initial clsuter result of mode 2.
# @param Es.init vector or NULL. Initial clsuter result of mode 3.
# @param Fs.init vector or NULL. Initial clsuter result of mode 4.
# @param nstart positive interger. The same as the "nstart" in kmeans().
# @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
# @param center logic value that indicates whether run "x = x-mean(x)" before performing clustering.
# @return a list   
# \code{judgeX} estimated underlying mean signal.   
# 
#                \code{Cs} clustering result of mode 1.  
#                
#                \code{Ds} clustering result of mode 2.  
#                
#                \code{Es} clustering result of mode 3.
#                
#                \code{Fs} clustering result of mode 4. 
#                
#                \code{mus} estimated underlying mean signal of each cluster.  



classify4 = function(x,k,r,l,m,lambda=0,max.iter=1000,threshold = 1e-15,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,Fs.init=NULL,nstart=20,method="L0",center=FALSE){
  #x=try$x;lambda=0;max.iter=1000;threshold=1e-15;trace=FALSE;Cs.init=NULL;Ds.init=NULL;Es.init=NULL;nstart=20;method="L0";center=FALSE
  n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]; s = dim(x)[4]
  if(center == TRUE) {
    mustemp <- mean(x)
    x <- x-mustemp
  }
  if(is.null(Cs.init)){
    if(k==1) Cs = rep(1,n) else {Cs  = kmeans(tensor_unfold4(x,1),k,nstart = nstart)$cluster}
  } else {
    Cs = Cs.init
  }
  if(is.null(Ds.init)){
    if(r==1) Ds = rep(1,p) else {Ds  = kmeans(tensor_unfold4(x,2),r,nstart = nstart)$cluster}
  } else {
    Ds = Ds.init
  }
  if(is.null(Es.init)){
    if(l==1) Es = rep(1,q) else {Es  = kmeans(tensor_unfold4(x,3),l,nstart = nstart)$cluster}
  } else {
    Es = Es.init
  }
  if(is.null(Es.init)){
    if(m==1) Es = rep(1,s) else {Fs  = kmeans(tensor_unfold4(x,4),m,nstart = nstart)$cluster}
  } else {
    Fs = Fs.init
  }
  objs <- 1e+15
  improvement <- 1e+10
  i <- 1
  mu.array = UpdateMus_tensor4(x,Cs,Ds,Es,Fs,lambda,method=method)
  while((improvement > threshold) & (i <= max.iter)){
    #Cs = UpdateClusters_tensor(tensor_unfold4(x),tensor_unfold4(mu.array),Cs,(rep(Ds,each=q)-1)*l+rep(Es,times=p))
    Cs = UpdateClusters_tensor(tensor_unfold4(x),tensor_unfold4(mu.array),Cs,
                               rep(Ds,times=q*s)+r*(rep(rep(Es,each=p),times=s)-1)+r*l*(rep(Fs,each=p*q)-1))
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    Cs <- ReNumber(Cs)
    k = length(unique(Cs))
    
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    #Ds = UpdateClusters_tensor(tensor_unfold4(x,2),tensor_unfold4(mu.array,2),Ds,(rep(Es,each=n)-1)*k+rep(Cs,times=q))
    Ds = UpdateClusters_tensor(tensor_unfold4(x,2),tensor_unfold4(mu.array,2),Ds,
                               rep(Cs,times=q*s)+k*(rep(rep(Es,each=n),times=s)-1)+k*l*(rep(Fs,each=n*q)-1))
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    Ds <- ReNumber(Ds)
    r = length(unique(Ds))
    
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    #Es = UpdateClusters_tensor(tensor_unfold4(x,3),tensor_unfold4(mu.array,3),Es,(rep(Ds,each=n)-1)*k+rep(Cs,times=p))
    Es = UpdateClusters_tensor(tensor_unfold4(x,3),tensor_unfold4(mu.array,3),Es,
                               rep(Cs,times=p*s)+k*(rep(rep(Ds,each=n),times=s)-1)+k*r*(rep(Fs,each=n*p)-1))
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    Es <- ReNumber(Es)
    l = length(unique(Es))
    
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    #Fs = UpdateClusters_tensor(tensor_unfold4(x,3),tensor_unfold4(mu.array,3),Fs,(rep(Ds,each=n)-1)*k+rep(Cs,times=p))
    Fs = UpdateClusters_tensor(tensor_unfold4(x,4),tensor_unfold4(mu.array,4),Fs,
                               rep(Cs,times=p*q)+k*(rep(rep(Ds,each=n),times=q)-1)+k*r*(rep(Es,each=n*p)-1))
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    Fs <- ReNumber(Fs)
    m = length(unique(Fs))
    
    mu.array = UpdateMus_tensor4(x,Cs,Ds,Es,Fs,lambda, method=method)
    objs <- c(objs, Objective4(x, mu.array, Cs, Ds, Es, Fs, lambda = lambda, method=method))
    improvement <- abs(objs[length(objs)] - objs[length(objs) - 
                                                   7])/abs(objs[length(objs) - 7])
    i <- i + 1
    if(trace) cat("step",i,",improvement=",improvement,".\n")
    if (is.na(improvement)) break
  }
  if (i > max.iter) {
    warning("The algorithm has not converged by the specified maximum number of iteration.\n")
  }
  if(center==TRUE){
    mu.array <- mu.array + mustemp
  }
  mu.array[abs(mu.array)<=1e-6] = 0
  return(list("judgeX"=mu.array[Cs,Ds,Es,Fs, drop=FALSE],"Cs"=Cs,"Ds"=Ds,"Es"=Es,"Fs"=Fs,"objs"=objs[length(objs)], "mus"=mu.array))
}

# Calculate the correct rate of selecting \eqn{d_1}, \eqn{d_2}, \eqn{d_3} 
# 
# Calculate the correct rate of the result of sparse_choosekrl().
# @param true Vector. The true number of clusters in each mode (eg. c(3,3,3)).
# @param results List. A list consists of 3*1 matrix and each matrix is the estimated c(d_1,d_2,d_3). (The return value of sparse_choosekrl().)
# @return A numeric value. The correct rate of selecting \eqn{d_1}, \eqn{d_2}, \eqn{d_3} 


Calculate<-function(true,results){
  real<-matrix(true,ncol=3)
  percent<-0
  for(i in 1:length(results)){
    if(nrow(results[[i]])>1){
      for(a in 1:nrow(results[[i]])){
        #  	  cat("iteration",a,fill=TRUE)
        if(sum(results[[i]][a,]==real)==3){
          percent<-percent+(1/nrow(results[[i]]))
        }
      }
    }
    else if(nrow(results[[i]])<2){
      if(sum(results[[i]]==real)==3){
        percent<-percent+1
      }	
    }
  }
  return(percent/length(results))
}

# Summarize the result of estimating the true c(d_1,d_2,d_3) in a simulation.
# 
# Summarize the return objct of sparse_choosekrl().
# @param results List. A list consists of 3-dimensional vectors and each vector is estimated c(d_1,d_2,d_3). (The result returned by sparse_choosekrl().)
# 
# @return \code{meank} mean estimated d_1
#         \code{meanr} mean estimated d_2
#         \code{meanl} mean estimated d_3
#         \code{sdk} the standard deviation of estimated d_1
#         \code{sdr} the standard deviation of estimated d_2
#         \code{sdl} the standard deviation of estimated d_3

Calculatekrl<-function(results){
  k<-rep(NA,length(results))
  r<-rep(NA,length(results))
  l<-rep(NA,length(results))
  #length(result): the number of samples
  for(i in 1:length(results)){
    #i = 1
    if(nrow(results[[i]]>1)){
      tempk<-0
      tempr<-0
      templ<-0
      for(a in 1:nrow(results[[i]])){
        #because the return value of Do function sometimes there are not only one classification result in one iteration
        tempk<-tempk+(results[[i]][a,1]/nrow(results[[i]]))
        tempr<-tempr+(results[[i]][a,2]/nrow(results[[i]]))
        templ<-templ+(results[[i]][a,3]/nrow(results[[i]]))
      }
      k[i]<-tempk
      r[i]<-tempr
      l[i]<-templ
    } else if(nrow(results[[i]]<2)){
      k[i]<-results[[i]][1,1]
      r[i]<-results[[i]][1,2]
      l[i]<-results[[i]][1,3]
    }
  }
  return(list(meank=mean(k),meanr=mean(r),meanl=mean(l),sdek=sd(k)/sqrt(length(k)),sder=sd(r)/sqrt(length(r)),sdel=sd(l)/sqrt(length(l)) ))
}



cer = function(bires,data){
  cerC<-1-adjustedRand(data$truthCs,bires$Cs,randMethod=c("Rand"))
  cerD<-1-adjustedRand(data$truthDs,bires$Ds,randMethod=c("Rand"))
  cerE<-1-adjustedRand(data$truthEs,bires$Es,randMethod=c("Rand"))
  cat("The CER(clustering error rate) is ",cerC,",",cerD,",",cerE,".\n")
  return(list(cerC = cerC, cerD = cerD, cerE = cerE))
}


# Calculate the clustering error rate (CER)
# 
# Calculate clustering error rate (CER) by comparing the clustering result and the true underlying mean signal tensor.
# @param bires List. The return value of label2() or classify2() which consists of judgeX, Cs, Ds, Es, objs, mus.
# @param data List. The return value of get.data().
# 
# @return A list which consists of the CER of all three modes.

# Draw 3d plot of tensor
# 
# Draw 3d plot of a given tensor.
# @param tensor a three-dimensional array
# @return None


#require two packages: rgl and viridis
# plot_tensor=function(tensor){
#   
#   
#   n=prod(dim(tensor))   
#   color_choice=min(round(prod(dim(tensor))/(6*6*6)),100)+1
#   marker=viridis_pal(option = "B")(color_choice)
#   
#   position=positionfun(dim(tensor))$position
#   quan=c(quantile(tensor,(0:color_choice)/color_choice))
#   col=tensor
#   for(i in 1:color_choice){
#     col[(tensor>=quan[i])&(tensor<=quan[i+1])]=marker[i]
#   }
#   col[tensor==quan[i+1]]=marker[i]
#   
#   plot3d(position[,1],position[,2],position[,3],col=col,alpha=0.3,size=5,xlab="",ylab="",zlab="")
# }


# Perform simulation: selecting the best \eqn{d_1}, \eqn{d_2} and \eqn{d_3}
# 
# Simulation: selecting the best \eqn{d_1}, \eqn{d_2} and \eqn{d_3}.
# @param n the dimension of mode 1
# @param p the dimension of mode 2
# @param q the dimension of mode 3
# @param k \eqn{d_1}: the clusters number of mode 1
# @param r \eqn{d_2}: the clusters number of mode 2
# @param l \eqn{d_3}: the clusters number of mode 3
# @param error positive numeric value. The noise when producing the data
# @param sim.times positive integer. Simulation times
# @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
# @param mode two options: "bic": selecting cluster numbers by bic; "crossvalidation":selecting cluster numbers by cross validation.
# @param seed logic value. Whether set seed to each simulation.
# @return A list consists of the estimating result in different iteration.


sim_choosekrl <- function(n,p,q,k,r,l,error=1,sim.times=5,method="L0",mode="bic",seed=TRUE){
  classification<-list()
  if(mode == "bic"){
    for(a in 1:sim.times){
      cat("Starting", a, fill=TRUE)
      if(seed == TRUE) set.seed(a)
      x = getOrder3Tensor(n,p,q,k,r,l,error=error)$x
      classification[[a]]<-chooseClusteringSize(x,k=2:5,r=2:5,l=2:5,method=method)$estimated_krl#estimate clusters
    }
  }
  if(mode == "crossvalidation"){
    for(a in 1:sim.times){
      cat("Starting", a, fill=TRUE)
      if(seed == TRUE) set.seed(a)
      x = getOrder3Tensor(n,p,q,k,r,l,error=error)$x
      classification[[a]]<-sparse_choosekrl(x,k=2:5,r=2:5,l=2:5,method=method)$estimated_krl#estimate clusters
    }
  }
  if(mode!="bic" & mode!="crossvalidation") stop("No such kind of mode:", mode, ".\n")
  return(classification)
}

# Perform simulation: selecting the lambda
# 
# 
# Simulation: selecting the lambda.  Select the lambda with lowest BIC while given a certain range.
# @param n the dimension of mode 1
# @param p the dimension of mode 2
# @param q the dimension of mode 3
# @param k \eqn{d_1}: the clusters number of mode 1
# @param r \eqn{d_2}: the clusters number of mode 2
# @param l \eqn{d_3}: the clusters number of mode 3
# @param sparse the sparse percent of data
# @param iteration iteration times
# @param lambda a vector of possible lambda, for example: lambda = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200)
# @param standarddeviation the standard deviation when producing data
# @param center if True, then x = x- mean(x)
# @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
# @return \code{selectedlambda}: a vector of selecting lambdas with lowest BIC in each iteration.

sim_chooseLambda = function(n,p,q,k,r,l,sparse,iteration,lambda,standarddeviation=4,center = FALSE,method="L0"){
  selectedLambda = 1:iteration
  for(iter in 1:iteration){
    cat("Iteration",iter,fill=TRUE)
    set.seed(iter)
    smp = getOrder3Tensor(n,p,q,k,r,l,error=standarddeviation,sort=FALSE,sparse.percent = sparse)$x
    if(center == TRUE)smp = smp - mean(smp)
    selectedLambda[iter] = chooseLambda(smp,k,r,l,lambda=lambda,method=method)$lambda
  }
  return(selectedLambda)
}

# Perform simulation: calculating the mean CER over iterations
# 
# Perform simulation: calculating the mean CER over iterations.
# @param n the dimension of mode 1
# @param p the dimension of mode 2
# @param q the dimension of mode 3
# @param k \eqn{d_1}: the clusters number of mode 1
# @param r \eqn{d_2}: the clusters number of mode 2
# @param l \eqn{d_3}: the clusters number of mode 3
# @param error the standard deviation when producing data
# @param lambda ...
# @param iteration iteration times
# @param method two options: "L0", "L1". Two methods use different penalties, where "L1" indicating Lasso penalty.
# @return A list consists of mean CER and standard deviation of CER across iterations.

simulation  = function(n,p,q,k,r,l,error,lambda,iteration=1,method="L0"){
  cer = c()
  for (i in 1:iteration){
    set.seed(i)
    data = getOrder3Tensor(n,p,q,k,r,l,error=error)
    test = data$x
    truthCs = data$truthCs
    truthDs = data$truthDs
    truthEs = data$truthEs
    sim = classify2(test,k,r,l,method=method)
    cerC<-1-adjustedRand(truthCs,sim$Cs,randMethod=c("Rand"))
    cerD<-1-adjustedRand(truthDs,sim$Ds,randMethod=c("Rand"))
    cerE<-1-adjustedRand(truthEs,sim$Es,randMethod=c("Rand"))
    cer = c(cer,cerC,cerD,cerE)
  }
  meancer = apply(matrix(cer,ncol=3,byrow=TRUE),MARGIN=2,mean)
  sdcer = apply(matrix(cer,ncol=3,byrow=TRUE),MARGIN=2,sd)
  return(list("meancer"=meancer,"sdcer"=sdcer))
}

# Perform tuning parameter (\eqn{d_1}, \eqn{d_2} and \eqn{d_3}) selection for sparse tensor clustering via cross validation
# 
# Select the best \eqn{d_1}, \eqn{d_2} and \eqn{d_3} to perform clustering. A range of values of d[1], d[2], d[3] is usually considered.
# @param x a three-dimensional array
# @param k the range of \eqn{d_1}: a vector, the possible clusters numbers of mode 1
# @param r the range of \eqn{d_2}: a vector, the possible clusters numbers of mode 2
# @param l the range of \eqn{d_3}: a vector, the possible clusters numbers of mode 3
# @param lambda a numeric value. The coefficient of the regularization term.
# @param percent a numeric value between 0 and 1
# @param trace  trace logic value. If true, it would print the iteration situation.
# @param nstart positive interger. The same as the "nstart" in kmeans().
# @param sim.times the same as tbmClustering(): the times of calling classify2() with different seeds.
# @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.
# @return a list   
# \code{estimated_krl} a 1*3 matrix which is the estimated c(d_1,d_2,d_3).  
# 
#                \code{results.se} the standard error of each possible combination of the cluster numbers.  
#                
#                \code{results.mean} the mean residual of each possible combination of the cluster numbers.  
#                

sparse_choosekrl = function (x,k,r,l,lambda=0,percent=0.2,trace=FALSE,nstart=20,sim.times=1,method="L0") {
  #x=test;l=range.l;lambda=0;percent=0.2;trace=TRUE
  #k=2:4;r=2:4;l=2:4
  if ((1%%percent) != 0) 
    stop("1 must be divisible by the specified percentage")
  if (percent <= 0) 
    stop("percentage cannot be less than or equal to 0")
  if (percent >= 1) 
    stop("percentage cannot be larger or equal to 1")
  
  #Returns suitably lagged and iterated differences.
  if (sum(diff(k) <= 0) > 0 || sum(diff(r) <= 0) > 0 || sum(diff(l) <= 0) > 0) 
    stop("k and r has to be an increasing sequence.  Please sort k and r before using the function")
  n=dim(x)[1];p=dim(x)[2];q=dim(x)[3]
  miss <- sample(1:(n*p*q), n*p*q, replace = FALSE)
  numberoftimes <- 1/percent
  allresults <- array(NA, dim = c(numberoftimes, length(k), 
                                  length(r), length(l)))
  Cs.init <- matrix(NA, nrow = dim(x)[1], ncol = length(k))
  #put the kmeans results into columns
  for (i in 1:length(k)) {
    Cs.init[, i] <- kmeans(tensor_unfold(x), k[i], nstart = nstart)$cluster
  }
  
  Ds.init <- matrix(NA, nrow = dim(x)[2], ncol = length(r))
  for (j in 1:length(r)) {
    Ds.init[, j] <- kmeans(tensor_unfold(x,2), r[j], nstart = nstart)$cluster
  }
  
  Es.init <- matrix(NA, nrow = dim(x)[3], ncol = length(l))
  for (j in 1:length(l)) {
    Es.init[, j] <- kmeans(tensor_unfold(x,3), l[j], nstart = nstart)$cluster
  }
  
  for (i in 1:numberoftimes) {#i = 1 
    if (trace == TRUE) 
      cat("Iteration", i, fill = TRUE)
    xmiss <- x
    missing <- miss[1:round(n*p*q*percent)]
    xmiss[missing] <- NA
    xmiss[missing] <- mean(xmiss, na.rm = TRUE)
    
    
    krl = matrix(c(rep(1:length(k),each=length(r)*length(l)),
                   rep(1:length(r),times=length(k)*length(l)),
                   rep(rep(1:length(l),each=length(r)),times=length(k))),byrow=TRUE,
                 nrow=3)
    
    
    if (.Platform$OS.type == "windows") {
      res = apply(krl,MARGIN=2,label_for_krl,k,r,l,Cs.init,Ds.init,Es.init,sim.times=sim.times,lambda=lambda,xmiss=xmiss,method=method)
    } else {
      krl_list = as.list(as.data.frame(krl))
      res = mclapply(krl_list, label_for_krl,k,r,l,Cs.init,Ds.init,Es.init,sim.times=sim.times,lambda=lambda,xmiss=xmiss,method=method,mc.cores=n.cores)
    }
    
    
    for (a in 1:dim(krl)[2]){
      #print(krl[,a])
      #print(sum((x[missing]-res[[a]][[1]][missing])^2))
      allresults[i,krl[,a][1],krl[,a][2],krl[,a][3]] = sum((x[missing]-res[[a]][[1]][missing])^2)
    }
    miss <- miss[-1:-(dim(x)[1] * dim(x)[2] * dim(x)[3]/numberoftimes)]
  }
  results.se <- apply(allresults, MARGIN=c(2, 3, 4), sd)/sqrt(numberoftimes)
  results.mean <- apply(allresults, c(2, 3, 4), mean)
  #comparing every mean with the mean with higher k,r,l
  IndicatorArray <- 1 * (results.mean[1:(length(k) - 1), 1:(length(r) - 
                                                              1), 1:(length(l) - 1)] <= results.mean[2:length(k), 2:length(r), 2:length(l)] + results.se[2:length(k), 2:length(r), 2:length(l)])
  if (max(IndicatorArray) == 0) {
    warning("The krl has reached the upper boundary. Please enlarger the range.")
    return(list(estimated_krl = c(max(k),max(r),max(l))))}
  
  #outer(outer(matrix(1:4,2),matrix(5:8,2)),matrix(1:4,2))
  ModeIndex <- outer(outer(k[-length(k)], r[-length(r)],  "*"), l[-length(l)], "*")
  smallestIndicatorTrue <- min(ModeIndex[IndicatorArray == TRUE])
  out <- which(IndicatorArray == TRUE & ModeIndex == 
                 smallestIndicatorTrue, arr.ind = TRUE)
  out <- array(c(k[out[,1]], r[out[,2]], l[out[,3]]), dim=dim(out))
  tempmode1 <- NULL
  tempmode2 <- NULL
  tempmode3 <- NULL
  for (i in 1:length(k)) {
    tempmode1 <- c(tempmode1, paste("K = ", k[i], sep = ""))
  }
  for (i in 1:length(r)) {
    tempmode2 <- c(tempmode2, paste("R = ", r[i], sep = ""))
  }
  for (i in 1:length(l)) {
    tempmode3 <- c(tempmode3, paste("L = ", l[i], sep = ""))
  }
  
  
  dimnames(results.se)[[1]] = tempmode1
  dimnames(results.se)[[2]] = tempmode2
  dimnames(results.se)[[3]] = tempmode3
  dimnames(results.mean)[[1]] = tempmode1
  dimnames(results.mean)[[2]] = tempmode2
  dimnames(results.mean)[[3]] = tempmode3
  return(list(estimated_krl = out, results.se = results.se, 
              results.mean = results.mean))
}

# Evaluate the accuracy of clustering result
# 
# Given the input tensor, perform tensor clustering and evaluate the accuracy of the clustering result. 
# @param bires the return value of tbmClustering()
# @param data the return value of getOrder3Tensor()
# @param CER logic value. If true, it would return CER also
# @param show logic value. If true, it would print the result
# @return a list:
#     \code{sparsityrate}: estimated sparsity rate;
#     
#     \code{correctzerorate}: correct zero rate;
#     
#     \code{correctonerate}: correct one rate;
#     
#     \code{totalincorrectrate}: total incorrect rate;
#     
#     \code{cerC}: CER on mode 1;
#     
#     \code{cerD}: CER on mode 2;
#     
#     \code{cerE}: CER on mode 3;
#     
#     \code{cerTotal}: overall CER;
#     
#     \code{error}: relative error.
#     

sparse_evaluate = function(bires, data, CER=TRUE, show=TRUE){
  npq=dim(data$x);n=npq[1];p=npq[2];q=npq[3]
  restensor<-(abs(bires$judgeX)>1e-10)*1
  totalzero<-sum(restensor==0)/(n*p*q)
  correctzero<-sum(restensor[which(data$binaryX==0)]==0)/sum(data$binaryX==0)
  correctone<-sum(restensor[which(data$binaryX==1)]==1)/(n*p*q-sum(data$binaryX==0))
  totalincorrect<-sum(abs(restensor-data$binaryX))/(n*p*q)
  error = sum((bires$judgeX-data$truthX)^2)/prod(dim(data$x))
  result = list(sparsityrate=totalzero, correctzerorate=correctzero,correctonerate=correctone,totalincorrectrate=totalincorrect,error=error)
  if (CER == TRUE){
    cerC<-1-adjustedRand(data$truthCs,bires$Cs,randMethod=c("Rand"))
    cerD<-1-adjustedRand(data$truthDs,bires$Ds,randMethod=c("Rand"))
    cerE<-1-adjustedRand(data$truthEs,bires$Es,randMethod=c("Rand"))
    blocktruth=cluster2block(data$mus,data$truthCs,data$truthDs,data$truthEs)
    blockest=cluster2block(bires$mus,bires$Cs,bires$Ds,bires$Es)
    
    cerTotal<-1-adjustedRand(as.numeric(as.factor(blocktruth)),as.numeric(as.factor(blockest)),randMethod=c("Rand"))
    #if (show == TRUE) cat("CerC is", cerC, ", cerD is", cerD, ", cerE is", cerE, ", cerTotal is", cerTotal, ".\n")
    if (show == TRUE) cat("CerC is", cerC, ", cerD is", cerD, ", cerE is", cerE, ".\n")
    result = list(sparsityrate=totalzero, correctzerorate=correctzero,correctonerate=correctone,totalincorrectrate=totalincorrect,cerC=cerC,cerD=cerD,cerE=cerE,cerTotal=cerTotal,error=error)
  }
  if(show == TRUE) cat("Total zero rate is", totalzero, ", correct zero rate is", correctzero, ", correct one rate is", correctone, ", total incorrect rate is", totalincorrect, ", MSE is", error,".\n")
  return(result)
}

# Calculate the BIC of clustering result
# 
# Given the clustering result, calculate the BIC.
# @param x a three-dimensional array
# @param clusterobj the return object of tbmClustering() or classify2()
# @param method two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty, "L1" indicating Lasso penalty.
# @param apply "main": apply in the main formula; "cp": apply in the CPD k-means.
# @return a vector [1]BIC, [2]nonzeromus
#
tensor_calculateBIC = function (x, clusterobj,method="L0",apply="main") 
{
  if(apply!="main" & apply!="cp") stop("parameter apply does not take such a value named", apply, ".\n")
  
  ## Modified; add the degree of freedom due to clustering. 
  npq = dim(x)
  n = npq[1]; p = npq[2]; q = npq[3]
  RSS=log(sum((x-clusterobj$judgeX)^2))*(n*p*q)
  if (apply == "main"){
    reducedCs=ifelse(dim(clusterobj$mus)[1]>1,dim(clusterobj$mus)[1],0)
    reducedDs=ifelse(dim(clusterobj$mus)[2]>1,dim(clusterobj$mus)[2],0)
    reducedEs=ifelse(dim(clusterobj$mus)[3]>1,dim(clusterobj$mus)[3],0)
    
    #df_clustering=min(reducedCs*log(n)+reducedDs*log(p)+reducedEs*log(q),sum(clusterobj$mus != 0)*log(n*p*q))
    df_clustering=reducedCs*log(n)+reducedDs*log(p)+reducedEs*log(q)
    
    df=log(n*p*q)*(sum(clusterobj$mus != 0)+df_clustering)
    #df=df_clustering
    return(RSS+df)
  }
  if (apply == "cp"){
    df = log(n*p*q)/n/p/q*(n+p+q-2)*clusterobj$s
    return(RSS+df)
  }
  
  
}

adjustedRand = function (cl1, cl2, randMethod = c("Rand", "HA", "MA", "FM", 
                                   "Jaccard")) 
{
  if (!is.vector(cl1)) {
    stop("cl1 is not a vector!\n")
  }
  if (!is.vector(cl2)) {
    stop("cl2 is not a vector!\n")
  }
  if (length(cl1) != length(cl2)) {
    stop("two vectors have different lengths!\n")
  }
  len <- length(randMethod)
  if (len == 0) {
    stop("The argument 'randMethod' is empty!\n")
  }
  cl1u <- unique(cl1)
  m1 <- length(cl1u)
  cl2u <- unique(cl2)
  m2 <- length(cl2u)
  n <- length(cl1)
  randVec <- rep(0, len)
  names(randVec) <- randMethod
  for (i in 1:len) {
    randMethod[i] <- match.arg(arg = randMethod[i], choices = c("Rand", 
                                                                "HA", "MA", "FM", "Jaccard"))
    flag <- match(randMethod[i], c("Rand", "HA", "MA", "FM", 
                                   "Jaccard"))
    c.res <- .C("adjustedRand", as.integer(cl1), as.integer(cl1u), 
                as.integer(cl2), as.integer(cl2u), as.integer(m1), 
                as.integer(m2), as.integer(n), as.integer(flag), 
                r = as.double(0))
    randVec[i] <- c.res$r
  }
  return(randVec)
}




tbmClustering = function(x,k,r,l,lambda=0,sym = F,max.iter=1000,threshold = 1e-10,sim.times=1,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,method="L0"){
  #x=test;lambda=1e-3;max.iter=200;threshold = 5e-3;sim.times=10
    if (sim.times == 1) return(classify2(x,k,r,l,sym,lambda=lambda,max.iter = max.iter,threshold = threshold,Cs.init = Cs.init,Ds.init = Ds.init,Es.init = Es.init,method=method))
    if (.Platform$OS.type == "windows") {
      result = lapply(rep(list(x),sim.times), classify2, k,r,l,sym,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init,method=method)
      objs = unlist(lapply(result, function(result){result$objs}))
    } else {
      result = mclapply(rep(list(x),sim.times), classify2, k,r,l,sym,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init,nstart = sample(1:1000,1),method=method,mc.cores = n.cores)
      objs = unlist(lapply(result, function(result){result$objs}))
    }
    result = result[[which(objs == min(objs))[1]]]

  return(result)
}
  
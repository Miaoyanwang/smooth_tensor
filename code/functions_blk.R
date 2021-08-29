# New smooth tensor estimation algorithms

# Block Average function
block = function(A,kvec){
  ouput = list()
  ftor = array(1:prod(kvec),dim = kvec)

  g = list()
  for(i in 1:length(dim(A))){
    g[[i]] = sort(rep(1:kvec[i],ceiling(dim(A)[i]/kvec[i]))[1:dim(A)[i]])
  }

  if(length(dim(A))==2){

    fmat = ftor[g[[1]],g[[2]]]
    Avg_mat = array(aggregate(c(A),by = list(as.factor(fmat)),mean)[,-1],dim = kvec)
    Blk_mat = Avg_mat[g[[1]],g[[2]]]

  }else if(length(dim(A))==3){

    fmat = ftor[g[[1]],g[[2]],g[[3]]]
    Avg_mat = array(aggregate(c(A),by = list(as.factor(fmat)),mean)[,-1],dim = kvec)
    Blk_mat = Avg_mat[g[[1]],g[[2]],g[[3]]]

  }else{
    stop("* Input array should be matrix or 3-order tensor")
  }

  return(list(Avg_mat = Avg_mat,Blk_mat= Blk_mat))
}
#Borda count estimation
Borda = function(A,kvec){
  d = dim(A)[1]
  if(length(dim(A))==2){
    #sorting step
    o1 = order(sapply(1:d, function(x) sum(A[x,])))
    As = A[o1,]

    #block approximation
    Blk = block(As,kvec)

    #sorting back
    Theta = Blk$Blk_mat[o1,]

  }else if(length(dim(A))==3){
    #sorting
    o1 = order(sapply(1:d, function(x) sum(A[x,,])))
    As = A[o1,,]

    #block approximation
    Blk = block(As,kvec)

    #sorting back
    Theta = Blk$Blk_mat[o1,,]
  }else{
    stop("* Input array should be matrix or 3-order tensor")
  }
  return(Theta)
}


## fit blockwise polynomial tensor
## polynomial degree l; block size k;

polytensor=function(tensor, l, k){
  d=dim(tensor)[1]
  est=array(dim=c(d,d,d))
  z=rep(1:k,rep(floor(d/k),k))
  z=c(z,rep(k,d-length(z)))
  for(i in 1:k){
    for(j in 1:k){
      for(q in 1:k){
        subtensor=tensor[which(z==i),which(z==j),which(z==q)]
        X1=c(slice.index(subtensor,1))
        X2=c(slice.index(subtensor,2))
        X3=c(slice.index(subtensor,3))
        if(l==0)
          fit=lm(c(subtensor)~1)
        else
          fit=lm(c(subtensor)~polym(X1,X2,X3,degree=l,raw=TRUE))
        est[which(z==i),which(z==j),which(z==q)]=fit$fitted.value
      }
    }
  }
  return(est)
}


# Borda count estimation
Borda2 = function(A,l,k){
    d = dim(A)[1]
    #sorting
    o1 = order(sapply(1:d, function(x) sum(A[x,,])))
    As = A[o1,,]
    
    #polynomial block approximation
    est = polytensor(As,l,k)
    
    #sorting back
    Theta = est[o1,,]
    
  return(Theta)
}



# Spectral method with threshold.
Spectral = function(A,row_idx,col_idx,threshold = NULL){
  A_unfolded = unfold(as.tensor(A),row_idx = row_idx,col_idx = col_idx)@data
  if(is.null(threshold)){
    threshold = sqrt(max(dim(A_unfolded))) 
  }
  Decomp = svd(A_unfolded)
  
  s = max(length(which(ifelse(Decomp$d>threshold,Decomp$d,0)>0)),1)
  D = diag(Decomp$d,s)
  Theta = Decomp$u[,1:s,drop = F] %*% D %*% t(Decomp$v[,1:s,drop = F])
  Theta = fold(Theta,row_idx = row_idx,col_idx = col_idx,dim(A))@data
  return(Theta)
}

# Least square estimation
library(rTensor)
tensor_unfold = function(tensor,dim=1){
  if (dim == 1) unfold = aperm(tensor,c(3,2,1))
  if (dim == 2) unfold = aperm(tensor,c(1,3,2))
  if (dim == 3) unfold = tensor
  unfold = apply(unfold,3,c)
  if (is.vector(unfold)) return(as.matrix(unfold)) else return(t(unfold))
}

ReNumber = function (Cs,sort = T){
  newCs <- rep(NA, length(Cs))
  if(sort==T){
    uniq <- unique(Cs)
  }else{
    uniq <- sort(unique(Cs)) 
  }
  for (i in 1:length(uniq)) {
    newCs[Cs == uniq[i]] <- i
  }
  return(newCs)
}

UpdateMus_tensor = function (x, Cs, Ds, Es) {
  Cs = ReNumber(Cs,sort = F); Ds = ReNumber(Ds,sort = F); Es = ReNumber(Es,sort = F)
  
  mus = array(NA, c(length(unique(Cs)), length(unique(Ds)), length(unique(Es))))
  d = dim(mus)
  for (k in 1:d[1]){
    for (r in 1:d[2]){
      for (l in 1:d[3]){
        mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l],na.rm = T)
      }
    }
  }
  return(mus)
}



HSC = function(x,k,l,r,nstart = 40,sym = F){
  result = list()
  u1 = svd(tensor_unfold(x,1))$u[,1:k]
  u2 = svd(tensor_unfold(x,2))$u[,1:l]
  u3 = svd(tensor_unfold(x,3))$u[,1:r]
  hu1 = svd(tensor_unfold(ttl(as.tensor(x),list(t(u2),t(u3)),ms = c(2,3))@data,1))$u[,1:k]
  hu2 = svd(tensor_unfold(ttl(as.tensor(x),list(t(u1),t(u3)),ms = c(1,3))@data,2))$u[,1:l]
  hu3 = svd(tensor_unfold(ttl(as.tensor(x),list(t(u1),t(u2)),ms = c(1,2))@data,3))$u[,1:r]
  Y1 = hu1%*%t(hu1)%*%tensor_unfold(ttl(as.tensor(x),list(t(hu2),t(hu3)),ms = c(2,3))@data,1)
  Y2 = hu2%*%t(hu2)%*%tensor_unfold(ttl(as.tensor(x),list(t(hu1),t(hu3)),ms = c(1,3))@data,2)
  Y3 = hu3%*%t(hu3)%*%tensor_unfold(ttl(as.tensor(x),list(t(hu1),t(hu2)),ms = c(1,2))@data,3)
  if(sym ==F){
    result$Cs  = kmeans(Y1,k,nstart = nstart)$cluster
    result$Ds  = kmeans(Y2,r,nstart = nstart)$cluster
    result$Es  = kmeans(Y3,l,nstart = nstart)$cluster  
  }else{
    result$Ds =result$Es =result$Cs = kmeans(Y1,k,nstart = nstart)$cluster
  }
  return(result)
}




Bal = function(A,k,max_iter = 100,threshold = 1,rep = 5){
  n = dim(A)[1]
  zlist = list()
  for(r in 1:rep){
    z = kmeans(tensor_unfold(A,1),k)$cluster
    E = matrix(nrow = n,ncol = k)
    xi = vector(length=k)
    iter = 0; improvement = n
    while((iter<=max_iter)&(improvement > threshold)){
      iter = iter+1
      eta = as.numeric(table(z))
      for(i in 1:n){
        for(a in 1:k){
          E[i,a] = sum(A[i,which(z==a),])
        }
      }
      for(a in 1:k){
        xi[a] = eta[a]*(n-eta[a])+eta[a]*(eta[a]-1)
      }
      nz = apply(E%*%diag(1/xi),1,which.max);nz
      improvement =  sum(abs(nz-z)>0);improvement
      z = nz
    }
    zlist[[r]]  = ReNumber(z)
    z = zlist[[which.max(unlist(lapply(zlist,function(x) length(unique(x)))))]]
  }
  return(z)
}


LSE = function(A,k,mode = 2,max_iter = 100,threshold = 1,rep = 5){
  n = dim(A)[1]
  if(mode==1){
    #Spectral membership estimation
    z = kmeans(tensor_unfold(A,1),k,nstart = 100)$cluster
  }else if (mode==2){
    # Balasubramanian estimation
    z = Bal(A,k,max_iter,threshold, rep)
  }else if (mode==3){
    # HSC membership estimation
    z = HSC(A,k,k,k,sym= T)$Cs
  }
  z = ReNumber(z)
  
  mu.array = UpdateMus_tensor(A,z,z,z)
  # this is for the technical error
  
  Theta = mu.array[z,z,z, drop=FALSE]
  return(Theta)
}





######### Simulation functions ############################################

f1 = function(a){
  return(a[1]*a[2]*a[3])
}
f2 = function(a){
  return(mean(a))
}
f3 = function(a){
  return(1/(1+exp(-3*sum(a^2))))
}

fb3 = function(a){
  return((a[1]^2+a[2]^2+a[3]^2)/(exp(cos(1/(a[1]^2+a[2]^2+a[3]^2)))))
}
f4 = function(a){
  return(log(1+max(a)))
}
f5=function(a){
  return(min(a)/exp(-max(a)-sqrt(a[1])-sqrt(a[2])-sqrt(a[3])))
}
fb5=function(a){
  return(exp(-min(a)-sqrt(a[1])-sqrt(a[2])-sqrt(a[3])))
}

simulation = function(d, mode = 1,sigma = 0.5,signal_level=10){
  tensor=array(dim=c(d,d,d))
  X1=c(slice.index(tensor,1))/d
  X2=c(slice.index(tensor,2))/d
  X3=c(slice.index(tensor,3))/d
  if(mode==1){
    signal = array(apply(cbind(X1,X2,X3),1,f1),dim=rep(d,3))
  }else if(mode==2){
    signal = array(apply(cbind(X1,X2,X3),1,f2),dim=rep(d,3))
  }else if(mode==3){
    signal = array(apply(cbind(X1,X2,X3),1,f3),dim=rep(d,3))
  }else if(mode==4){
    signal = array(apply(cbind(X1,X2,X3),1,f4),dim=rep(d,3))
  }else if(mode==5){
    signal = array(apply(cbind(X1,X2,X3),1,f5),dim=rep(d,3))
  }
  ### edited by Miaoyan
  signal=signal_level*signal/sqrt(mean(signal^2)) ## normalize signal tensor to have averaged magnitude = signal_level
  observe = signal+rnorm(d^3,0,sigma)

  return(list(signal= signal,observe=observe))
}



Observe_A = function(P){
  n = dim(P)[1]; m = length(dim(P))
  U=array(rnorm(length(P),0,1),dim(P)) 
  tempA=1*(U<qnorm(P,0,1))
  A = array(0,dim(P))
  for (i in 1:n){
    for (j in i:n){
      for(k in j:n){
        A[i,j,k] = A[i,k,j] = A[j,i,k]  = A[j,k,i] = A[k,i,j] = A[k,j,i]= tempA[i,j,k]
      }
    }
  }
  return(A)
}


simulation_bin = function(d,mode = 1){
  tensor=array(dim=c(d,d,d))
  X1=c(slice.index(tensor,1))/d
  X2=c(slice.index(tensor,2))/d
  X3=c(slice.index(tensor,3))/d
  if(mode==1){
    signal = array(apply(cbind(X1,X2,X3),1,f1),dim=rep(d,3))
  }else if(mode==2){
    signal = array(apply(cbind(X1,X2,X3),1,f2),dim=rep(d,3))
  }else if(mode==3){
    signal = pmin(array(apply(cbind(X1,X2,X3),1,f3),dim=rep(d,3)),1)
  }else if(mode==4){
    signal = array(apply(cbind(X1,X2,X3),1,f4),dim=rep(d,3))
  }else if(mode==5){
    signal = array(apply(cbind(X1,X2,X3),1,fb5),dim=rep(d,3))
  }
  
  observe = Observe_A(signal)
  return(list(signal = signal, observe = observe))
}



summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}












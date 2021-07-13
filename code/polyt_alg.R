##################### Polynomial time algorithms #####################
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


# Sort and smooth method
SASmooth = function(A,kvec){
  d = dim(A)[1]
  if(length(dim(A))==2){
    #sorting
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



# Spectral method with threshold.
Spectral = function(A,row_idx,col_idx){
  A_unfolded = unfold(as.tensor(A),row_idx = row_idx,col_idx = col_idx)@data
  
  threshold = sqrt(max(dim(A_unfolded)))
  Decomp = svd(A_unfolded)
  
  s = max(length(which(ifelse(Decomp$d>threshold,Decomp$d,0)>0)),1)
  D = diag(Decomp$d,s)
  Theta = Decomp$u[,1:s,drop = F] %*% D %*% t(Decomp$v[,1:s,drop = F])
  Theta = fold(Theta,row_idx = row_idx,col_idx = col_idx,dim(A))@data
  return(Theta)
}



# High-order spectral method without threshold.
Hspectral = function(A,rk){
  m = length(dim(A))
  if(m==2){
    k = rk[1]
    u1 = svd(A)$u[,1:k]
    u2 = svd(t(A))$u[,1:k]
    hu1 = svd(A%*%u2)$u[,1:k]
    hu2 = svd(t(A)%*%u1)$u[,1:k]
    Theta = hu1%*%t(hu1)%*%A%*%hu2%*%t(hu2)
  }else if(m==3){
    k = rk[1]; l = rk[2]; r = rk[3]
    u1 = svd(tensor_unfold(A,1))$u[,1:k]
    u2 = svd(tensor_unfold(A,2))$u[,1:l]
    u3 = svd(tensor_unfold(A,3))$u[,1:r]
    hu1 = svd(tensor_unfold(ttl(as.tensor(A),list(t(u2),t(u3)),ms = c(2,3))@data,1))$u[,1:k]
    hu2 = svd(tensor_unfold(ttl(as.tensor(A),list(t(u1),t(u3)),ms = c(1,3))@data,2))$u[,1:l]
    hu3 = svd(tensor_unfold(ttl(as.tensor(A),list(t(u1),t(u2)),ms = c(1,2))@data,3))$u[,1:r]
    
    Theta = ttl(as.tensor(A),list(hu1%*%t(hu1),hu2%*%t(hu2),hu3%*%t(hu3)),ms = c(1,2,3))@data
  }else{
    stop("* Input array should be matrix or 3-order tensor")
  }
  
  
  return(Theta)
}


## High-order spectral method combining rank and threshold techniques
Hspectral2 = function(A){
  m = length(dim(A)); 
  if(m==2){
    d = max(dim(A))
    k1 = sqrt(d)
    k2 = sum(svd(A)$d>= sqrt(d))
    k = min(k1,k2)
    
    u1 = svd(A)$u[,1:k]
    u2 = svd(t(A))$u[,1:k]
    hu1 = svd(A%*%u2)$u[,1:k]
    hu2 = svd(t(A)%*%u1)$u[,1:k]
    Theta = hu1%*%t(hu1)%*%A%*%hu2%*%t(hu2)
  }else if(m==3){
    d = dim(A)
    k1 = sqrt(d[1]); l1 = sqrt(d[2]); r1 = sqrt(d[3])
    k2 = sum(svd(tensor_unfold(A,1))$d>= d[1]^(3/4))
    l2 = sum(svd(tensor_unfold(A,2))$d>= d[2]^(3/4))
    r2 = sum(svd(tensor_unfold(A,3))$d>= d[3]^(3/4))
    k = min(k1,k2); l = min(l1,l2); r = min(r1,r2)
    
    u1 = svd(tensor_unfold(A,1))$u[,1:k]
    u2 = svd(tensor_unfold(A,2))$u[,1:l]
    u3 = svd(tensor_unfold(A,3))$u[,1:r]
    hu1 = svd(tensor_unfold(ttl(as.tensor(A),list(t(u2),t(u3)),ms = c(2,3))@data,1))$u[,1:k]
    hu2 = svd(tensor_unfold(ttl(as.tensor(A),list(t(u1),t(u3)),ms = c(1,3))@data,2))$u[,1:l]
    hu3 = svd(tensor_unfold(ttl(as.tensor(A),list(t(u1),t(u2)),ms = c(1,2))@data,3))$u[,1:r]
    
    Theta = ttl(as.tensor(A),list(hu1%*%t(hu1),hu2%*%t(hu2),hu3%*%t(hu3)),ms = c(1,2,3))@data
  }else{
    stop("* Input array should be matrix or 3-order tensor")
  }
  
  
  return(Theta)
}

##################### simple test #####################

hgm = hgmodel.smooth(10)
A = hgm$A
P = hgm$P


## Block approximation method
# Calculating optimal block size
m = length(dim(A))
k = round(d^{m/(m+2)})
mean((P-SASmooth(A,c(k,k,k)))^2)

## Spectral method
mean((P-Spectral(A,1,c(2,3)))^2)


## Hspectral method
mean((P-Hspectral(A,c(k,k,k)))^2)

source("functions_hg.R") 
source("functions_sbm.R")
ini=HSC(A,k,k,k,sym=T)
res =tbmClustering(A,k,k,k,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = T, diagP = F)
mean((res$judgeX-P)^2)




##################### Smooth graphon #####################

Dat = as.data.frame(matrix(nrow = 10,ncol = 4))
colnames(Dat) = c("d","SAS","Spectral","SBM")
Dat[,1] = 10*(1:10)

for(d in 1:10){
  hgm = hgmodel.smooth(10*d)
  A = hgm$A
  P = hgm$P
  ## Block approximation method
  # Calculating optimal block size
  m = length(dim(A))
  k = round((10*d)^{m/(m+2)})
  Dat[d,2] = mean((P-SASmooth(A,c(k,k,k)))^2)
  tic()
  ini=HSC(A,k,k,k,sym=T)
  res =tbmClustering(A,k,k,k,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = T, diagP = F)
  toc()
  ## Spectral method
  Dat[d,3] = mean((P-Spectral(A,1,c(2,3)))^2)
  

  ini=HSC(A,k,k,k,sym=T)
  res =tbmClustering(A,k,k,k,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = T, diagP = F)
  Dat[d,4] = mean((res$judgeX-P)^2)
  
}


dat = as.data.frame(matrix(nrow = 30,ncol = 3))
colnames(dat) = c("dim","method","MSE")
dat[,1] = rep(Dat[,1],3)
dat[,2] = rep(c("SAS","Spectral","SBM"),each = 10)
dat[,3] = c(Dat[,2],Dat[,3],Dat[,4])

library(ggplot2)
g1 = ggplot(data = dat, aes(x = dim,y =MSE,color = method ))+geom_point()+geom_line()+
  ggtitle("Smooth graphon")




##################### Piece-wise constant graphon  #####################
Dat = as.data.frame(matrix(nrow = 10,ncol = 4))
colnames(Dat) = c("d","SAS","Spectral","SBM")
Dat[,1] = 10*(1:10)

for(d in 1:10){
  W = pbtensor(d*2,sym = T,smooth = T)
  ## results seem invariant to cluster permutations within each mode.
  #shuffle=sample(1:g,g,replace=F)##
  #W=W[shuffle,shuffle,shuffle]
  ## image(W[,,1])
  hgmod = hgmodel.block(W,d*10,diagP = F,type="Gaussian")
  A = hgmod$A
  P = hgmod$P
  ## Block approximation method
  # Calculating optimal block size
  m = length(dim(A))
  k = round((10*d)^{m/(m+2)})
  Dat[d,2] = mean((P-SASmooth(A,c(k,k,k)))^2)
  
  ## Spectral method
  Dat[d,3] = mean((P-Spectral(A,1,c(2,3)))^2)
  
  
  ini=HSC(A,k,k,k,sym=T)
  res =tbmClustering(A,k,k,k,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = T, diagP = F)
  Dat[d,4] = mean((res$judgeX-P)^2)
  
}


dat = as.data.frame(matrix(nrow = 30,ncol = 3))
colnames(dat) = c("dim","method","MSE")
dat[,1] = rep(Dat[,1],3)
dat[,2] = rep(c("SAS","Spectral","SBM"),each = 10)
dat[,3] = c(Dat[,2],Dat[,3],Dat[,4])


g2 = ggplot(data = dat, aes(x = dim,y =MSE,color = method ))+geom_point()+geom_line()+
  ggtitle("Piecewise-constant graphon")






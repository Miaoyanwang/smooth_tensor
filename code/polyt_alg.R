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
  
  s = length(which(ifelse(Decomp$d>threshold,Decomp$d,0)>0))
  D = diag(Decomp$d,s)
  Theta = Decomp$u[,1:s,drop = F] %*% D %*% t(Decomp$v[,1:s,drop = F])
  Theta = fold(Theta,row_idx = row_idx,col_idx = col_idx,dim(A))@data
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






## produce probability tensor
pbtensor = function(K){
  W = array(0,c(K,K,K))
  for (i in 1:K){
    for (j in 1:K){
      for(k in 1:K){
        W[i,j,k] = runif(1)
      }
    }
  }
  return(W)
}

## cut function
cut = function(tnsr){
  cond = (dim(tnsr)[1]==dim(tnsr)[2])&(dim(tnsr)[1]==dim(tnsr)[3])
  if (!(cond)){
    stop("* tensor must have same dimension for each mode")
  }
  n = dim(tnsr)[1]
  for(i in 1:n){
    for(j in 1:n){
      for(k in 1:n)
        if((i-j)*(j-k)*(k-i)==0){
          tnsr[i,j,k] = 0
        }
    }
  }
  return(tnsr)
}


## given probability tensor, generate adjacenty tensor.
hgmodel.block = function(W,n,order = T){
  ## Check W
  cond1 = ((all(W>=0))&&(all(W<=1)))
  cond2 = (dim(W)[1]==dim(W)[2])&(dim(W)[1]==dim(W)[3])
  if (!(cond1&&cond2)){
    stop("* hgmodel.block : W is not a valid stochastic block model.")
  }
  K = dim(W)[1]
  if (K>n){
    stop("* hgmodel.block : the number of nodes should be >= the number of blocks.")
  }
  if(order == T){
    v = sort(runif(n))
  }else{
    v = runif(n)
  }
  
  u = round((K-1)*v)+1
  
  ## Construct Probability tensor and  a random hyper graph
  P = array(0,c(n,n,n)); A = array(0,c(n,n,n))
  for (i in 1:(n-2)){
    for (j in (i+1):(n-1)){
      for(k in (j+1):n){
        p = W[u[i],u[j],u[k]]
        P[i,j,k] = P[i,k,j] = P[j,i,k] = P[j,k,i] = P[k,i,j] = P[k,j,i] = p
        A[i,j,k] = A[i,k,j] = A[j,i,k] = A[j,k,i] = A[k,i,j] = A[k,j,i] = rbinom(1,1,p)
      }
    }
  }
  
  ## output
  output = list()
  output$A = A
  output$P = P
  return(output)
}

f1 = function(x,c) return(1/(1+exp(-c*sum(x^2))))

hgmodel.smooth = function(n,c,option = 1){
  if (option ==1){
    xi = runif(n)
    P = array(0,c(n,n,n)); A = array(0,c(n,n,n))
    
    for (i in 1:(n-2)){
      for (j in (i+1):(n-1)){
        for(k in (j+1):n){
          x = c(xi[i],xi[j],xi[k])
          P[i,j,k] = P[i,k,j] = P[j,i,k] = P[j,k,i] = P[k,i,j] = P[k,j,i] = f1(x,c)
          A[i,j,k] = A[i,k,j] = A[j,i,k] = A[j,k,i] = A[k,i,j] = A[k,j,i] = rbinom(1,1,f1(x,c))
        }
      }
    }
  }
  
  
  
  ## output
  output = list()
  output$A = A
  output$P = P
  return(output)
}


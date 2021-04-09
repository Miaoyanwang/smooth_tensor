## produce probability tensor
pbtensor = function(K,sym = T,smooth = T){
  ####################################
  ## Update by Miaoyan. Replace loop by array operation
  ####################################
  if(smooth ==T){
    temp=sort(runif(K*(K-1)*(K-2)/6+K*(K-1)+K,0,1)) ## sorted
    s=0
    if(sym==T){
      W = array(0,c(K,K,K))
      for (i in 1:K){
        for (j in i:K){
          for(k in j:K){
            s=s+1
            W[i,j,k] = W[i,k,j] = W[j,i,k]  = W[j,k,i] = W[k,i,j] = W[k,j,i]=temp[s]
          }
        }
      }
    }else{
      W=array(sort(runif(K^3,0,1)),c(K,K,K)) 
    }
  }else{
    if(sym==T){
      W = array(0,c(K,K,K))
      for (i in 1:K){
        for (j in i:K){
          for(k in j:K){
            W[i,j,k] = W[i,k,j] = W[j,i,k]  = W[j,k,i] = W[k,i,j] = W[k,j,i]=runif(1)
          }
        }
      }
    }else{
      W=array(runif(K^3,0,1),c(K,K,K)) 
    }    
  }

  return(W)
}


symmetrize=function(W){
  W_sym=W
  perm_group=rbind(c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
  for(s in 1:nrow(perm_group)){
    W_sym=W_sym+aperm(W,perm=perm_group[s,])
  }
  return(W_sym/6)
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





## given probability tensor, generate adjacenty tensor.
hgmodel.block = function(W,n,order = T,diagP = T,type="Bernoulli"){
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
  
  ####################################
  ## Update by Miaoyan. Loop is slow. should be avoided if possible.
  ####################################
  ## Construct Probability tensor and  a random hyper graph
  # P = array(0,c(n,n,n)); A = array(0,c(n,n,n))
  # for (i in 1:(n-2)){
  # for (j in (i+1):(n-1)){
  #   for(k in (j+1):n){
  #     p = W[u[i],u[j],u[k]]
  #     P[i,j,k] = P[i,k,j] = P[j,i,k] = P[j,k,i] = P[k,i,j] = P[k,j,i] = p
  #       A[i,j,k] = A[i,k,j] = A[j,i,k] = A[j,k,i] = A[k,i,j] = A[k,j,i] = rbinom(1,1,p)
  #    }
  #  }
  # }
  
  ## replace loop by array operation
  P=symmetrize(W[u,u,u]) ## symmetric probability tensor
  if(diagP==F){
    P = cut(P)  
  }
  
  
  U=array(rnorm(n^3,0,1),c(n,n,n)) ## i.i.d. Gaussian tensor
  U=symmetrize(U) ## symmetric Gaussian tensor. i.i.d. in subtensor
  if(type=="Bernoulli"){
  A=1*(U<qnorm(P,0,1/sqrt(6))) ## an equivalent way of generating Bernoulli entries with prob P. Both U and P are symmetric, so A is symmetric
  }else if(type=="Gaussian"){
      A=P+0.1*U## gaussian noise
  }

  
  ## output
  output = list()
  output$A = A
  output$P = P
  output$Cs = u
  return(output)
}

f1 = function(x,c) return(1/(1+exp(-(c*sum(x^2)))))

hgmodel.smooth = function(n,c,option = 1, diagP = T){
  if (option ==1){
    xi = runif(n)
    P = array(apply(expand.grid(xi, xi, xi), 1, function(x) f1(x,c)), dim=c(n,n,n))
    if(diagP==F){
      P = cut(P)
    }
    
    U=array(rnorm(n^3,0,1),c(n,n,n)) ## i.i.d. Gaussian tensor
    U=symmetrize(U)   
    A=1*(U<qnorm(P,0,1/sqrt(6)))
    # for (i in 1:(n-2)){
    #   for (j in (i+1):(n-1)){
    #     for(k in (j+1):n){
    #       x = c(xi[i],xi[j],xi[k])
    #       P[i,j,k] = P[i,k,j] = P[j,i,k] = P[j,k,i] = P[k,i,j] = P[k,j,i] = f1(x,c)
    #       A[i,j,k] = A[i,k,j] = A[j,i,k] = A[j,k,i] = A[k,i,j] = A[k,j,i] = rbinom(1,1,f1(x,c))
    #     }
    #   }
    # }
  }

  
  ## output
  output = list()
  output$A = A
  output$P = P
  return(output)
}




#### signal + noise model #########################

signaltensor = function(K,A,sym = F,smooth = T){

  if(smooth == T){
    if(sym==T){
      temp=sort(runif(K*(K-1)*(K-2)/6+K*(K-1)+K,-A,A))
      s = 0
      W = array(0,c(K,K,K))
      for (i in 1:K){
        for (j in i:K){
          for(k in j:K){
            s = s+1
            W[i,j,k] = W[i,k,j] = W[j,i,k]  = W[j,k,i] = W[k,i,j] = W[k,j,i] =temp[s]
          }
        }
      }
    }else{
      W=array(sort(runif(K^3,-A,A)),c(K,K,K)) 
    }
  }else{
    if(sym==T){
      W = array(0,c(K,K,K))
      for (i in 1:K){
        for (j in i:K){
          for(k in j:K){
            W[i,j,k] = W[i,k,j] = W[j,i,k]  = W[j,k,i] = W[k,i,j] = W[k,j,i] =runif(1,-A,A)
          }
        }
      }
    }else{
      W=array(runif(K^3,-A,A),c(K,K,K)) 
    }
  }
  
  return(W)
}



snmod.block = function(W,n,order = T){
  K = dim(W)[1]
  
  if(order == T){
    v1 = sort(runif(n));v2 = sort(runif(n));v3 = sort(runif(n))
  }else{
    v1 = runif(n);  v2 = runif(n); v3 = runif(n)
  }
  
  
  u1 = round((K-1)*v1)+1;u2 = round((K-1)*v2)+1;u3 = round((K-1)*v3)+1
  
  
  noise=array(rnorm(n^3,0,1),c(n,n,n)) ## i.i.d. Gaussian noise tensor
  obsTensor = W[u1,u2,u3]+noise
  ## output
  output = list()
  output$noise = noise
  output$obsTensor = obsTensor
  output$Sig = obsTensor-noise
  output$Cs = u1;output$Ds = u2; output$Es = u3
  
  return(output)
}








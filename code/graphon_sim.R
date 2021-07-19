Observe_A = function(P){
  n = dim(P)[1]; m = length(dim(P))
  U=array(rnorm(length(P),0,1),dim(P)) 
  tempA=1*(U<qnorm(P,0,1))
  A = array(0,dim(P))
  if(m==2){
    tempA[lower.tri(tempA)] = t(tempA)[lower.tri(tempA)]
    A = tempA
  }else if(m==3){
    for (i in 1:n){
      for (j in i:n){
        for(k in j:n){
          A[i,j,k] = A[i,k,j] = A[j,i,k]  = A[j,k,i] = A[k,i,j] = A[k,j,i]= tempA[i,j,k]
        }
      }
    }
  }else{
    stop("* mode should be 2 or 3")
  }
  return(A)
}




# Graphon 1.
graphon1 = function(n,m = 2){
  K = floor(log(n))
  if(m==2){
    P = array(rep(0.3/(K+1),n^2),dim = c(n,n) )
    for(k in 1:K){
      P[ceiling(n*(k-1)/K):floor(n*k/K),ceiling(n*(k-1)/K):floor(n*k/K)] = k/(K+1)
    }
    A = Observe_A(P)
    
  }else if(m==3){
    P = array(rep(0.3/(K+1),n^3),dim = c(n,n,n) )
    for(k in 1:K){
      P[ceiling(n*(k-1)/K):floor(n*k/K),ceiling(n*(k-1)/K):floor(n*k/K),
        ceiling(n*(k-1)/K):floor(n*k/K)] = k/(K+1)
    }
    A = Observe_A(P)

  }else{
    stop("* mode should be 2 or 3")
  }
  return(list(A = A,P = P))
}



# Graphon 2.
graphon2 = function(n,m = 2){
  if(m==2){
    P = array(0,dim = c(n,n))
    for(i in 1:n){
      for(j in i:n){
        P[i,j] = P[j,i] = sin(5*pi*((i+j)/n-1)+1)/2+0.5
      }
    }
    A = Observe_A(P)
  }else if(m==3){
    P = array(0,dim = c(n,n,n) )
    for(i in 1:n){
      for(j in 1:n){
        for(k in 1:n){
          P[i,j,k] = P[i,k,j] = P[j,i,k]  = P[j,k,i] = P[k,i,j] = P[k,j,i]= sin(5*pi*((i+j+k)/n-1)+1)/2+0.5
        }
      }
    }
    A = Observe_A(P)
    
  }else{
    stop("* mode should be 2 or 3")
  }
  return(list(A = A,P = P))
}


# Graphon 3.
graphon3 = function(n,m = 2){
  if(m==2){
    P = array(0,dim = c(n,n))
    for(i in 1:n){
      for(j in i:n){
        P[i,j] = P[j,i] = min(((i/n)^2+(j/n)^2)/(3**cos(1/((i/n)^2+(j/n)^2))),1)
      }
    }
    A = Observe_A(P)
  }else if(m==3){
    P = array(0,dim = c(n,n,n) )
    for(i in 1:n){
      for(j in 1:n){
        for(k in 1:n){
          P[i,j,k] = P[i,k,j] = P[j,i,k]  = P[j,k,i] = P[k,i,j] = P[k,j,i]= min(((i/n)^2+(j/n)^2+(k/n)^2)/(3**cos(1/((i/n)^2+(j/n)^2+(k/n)^2))),1)
        }
      }
    }
    A = Observe_A(P)
    
  }else{
    stop("* mode should be 2 or 3")
  }
  return(list(A = A,P = P))
}



# Graphon 4.
graphon4 = function(n,m = 2){
  if(m==2){
    P = array(0,dim = c(n,n))
    for(i in 1:n){
      for(j in i:n){
        P[i,j] = P[j,i] = max(i,j)/n
      }
    }
    A = Observe_A(P)
  }else if(m==3){
    P = array(0,dim = c(n,n,n) )
    for(i in 1:n){
      for(j in 1:n){
        for(k in 1:n){
          P[i,j,k] = P[i,k,j] = P[j,i,k]  = P[j,k,i] = P[k,i,j] = P[k,j,i]= max(i,j,k)/n
        }
      }
    }
    A = Observe_A(P)
    
  }else{
    stop("* mode should be 2 or 3")
  }
  return(list(A = A,P = P))
}


n = 50
m = 2
gp1 = graphon1(n,m)
gp2 = graphon2(n,m)
gp3 = graphon3(n,m)
gp4 = graphon4(n,m)




A1 = gp1$A; A2 = gp2$A; A3 = gp3$A; A4 = gp4$A
P1 = gp1$P; P2 = gp2$P; P3 = gp3$P; P4 = gp4$P
par(mar=c(5.1, 4.1, 4.1, 4.1)) 
library(RColorBrewer)
plot(P1,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=0.11*(0:9))
plot(P2,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=0.11*(0:9))
plot(P3,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=0.11*(0:9))
plot(P4,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=0.11*(0:9))


k = round((n)^{m/(m+2)})

P11 = SASmooth(A1,c(k,k))
P12 = Spectral(A1,1,2)
P13 = Hspectral(A1,c(k,k))
P21 = SASmooth(A2,c(k,k))
P22 = Spectral(A2,1,2)
P23 = Hspectral(A2,c(k,k))
P31 = SASmooth(A3,c(k,k))
P32 = Spectral(A3,1,2)
P33 = Hspectral(A3,c(k,k))
P41 = SASmooth(A4,c(k,k))
P42 = Spectral(A4,1,2)
P43 = Hspectral(A4,c(k,k))

trunc = function(x){
  return(ifelse(x>1,1,ifelse(x<0,0,x)))
}

################## P1: graphon1 analysis ###########################################
mean((P1-P11)^2)
mean((P1-P12)^2)
mean((P1-trunc(P12))^2)
mean((P1-trunc(P13))^2)
mean((P1-P13)^2)

mean((P2-P21)^2)
mean((P2-P22)^2)
mean((P2-trunc(P22))^2)
mean((P2-trunc(P23))^2)
mean((P2-P23)^2)

mean((P3-P31)^2)
mean((P3-P32)^2)
mean((P3-trunc(P32))^2)
mean((P3-trunc(P33))^2)
mean((P3-P33)^2)



par(mfrow = c(1,1))
pP11 = matrix(0,nrow = n,ncol =n) ;pP11[lower.tri(pP11)] = P11[lower.tri(P11)]; pP11[upper.tri(pP11)] = P1[upper.tri(P1)]
pP12 = matrix(0,nrow = n,ncol =n) ;pP12[lower.tri(pP12)] = P12[lower.tri(P12)]; pP12[upper.tri(pP12)] = P1[upper.tri(P1)]
pP13 = matrix(0,nrow = n,ncol =n) ;pP13[lower.tri(pP13)] = P13[lower.tri(P13)]; pP13[upper.tri(pP13)] = P1[upper.tri(P1)]

plot(pP11,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")
plot(pP12,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")
plot(pP13,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")

################## P2: graphon2 analysis ###########################################
mean((P2-P21)^2)
mean((P2-P22)^2)
mean((P2-P23)^2)

pP21 = matrix(0,nrow = n,ncol =n) ;pP21[lower.tri(pP21)] = P21[lower.tri(P21)]; pP21[upper.tri(pP21)] = P2[upper.tri(P2)]
pP22 = matrix(0,nrow = n,ncol =n) ;pP22[lower.tri(pP22)] = P22[lower.tri(P22)]; pP22[upper.tri(pP22)] = P2[upper.tri(P2)]
pP23 = matrix(0,nrow = n,ncol =n) ;pP23[lower.tri(pP23)] = P23[lower.tri(P23)]; pP23[upper.tri(pP23)] = P2[upper.tri(P2)]

plot(pP21,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")
plot(pP22,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")
plot(pP23,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")


################## P3: graphon3 analysis ###########################################
mean((P3-P31)^2)
mean((P3-P32)^2)
mean((P3-P33)^2)

pP31 = matrix(0,nrow = n,ncol =n) ;pP31[lower.tri(pP31)] = P31[lower.tri(P31)]; pP31[upper.tri(pP31)] = P3[upper.tri(P3)]
pP32 = matrix(0,nrow = n,ncol =n) ;pP32[lower.tri(pP32)] = P32[lower.tri(P32)]; pP32[upper.tri(pP32)] = P3[upper.tri(P3)]
pP33 = matrix(0,nrow = n,ncol =n) ;pP33[lower.tri(pP33)] = P33[lower.tri(P33)]; pP33[upper.tri(pP33)] = P3[upper.tri(P3)]

plot(pP31,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")
plot(pP32,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")
plot(pP33,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")


################## P4: graphon4 analysis ###########################################
mean((P4-P41)^2)
mean((P4-P42)^2)
mean((P4-P43)^2)

pP41 = matrix(0,nrow = n,ncol =n) ;pP41[lower.tri(pP41)] = P41[lower.tri(P41)]; pP41[upper.tri(pP41)] = P4[upper.tri(P4)]
pP42 = matrix(0,nrow = n,ncol =n) ;pP42[lower.tri(pP42)] = P42[lower.tri(P42)]; pP42[upper.tri(pP42)] = P4[upper.tri(P4)]
pP43 = matrix(0,nrow = n,ncol =n) ;pP43[lower.tri(pP43)] = P43[lower.tri(P43)]; pP43[upper.tri(pP43)] = P4[upper.tri(P4)]

plot(pP41,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")
plot(pP42,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")
plot(pP43,border = NA,col =brewer.pal(n =9, name = "Blues"),breaks=c(0.11*(0:8),1),na.col="gray")


################## hypergrapon analysis (updated) ############################################
Dat1 = Dat2 = Dat3 = Dat4 = as.data.frame(matrix(nrow = 10,ncol = 6))
colnames(Dat1) = colnames(Dat2) = colnames(Dat3) =colnames(Dat4) = c("d","SAS","Spectral","HSpectral","SBM","HSpectral2")
Dat1[,1]= Dat2[,1] = Dat3[,1] = Dat4[,1] = 10*(1:10)

for(d in 1:10){
  n = 10*d
  m = 3
  k = round((n)^{m/(m+2)})
  gp1 = graphon1(n,m)
  gp2 = graphon2(n,m)
  gp3 = graphon3(n,m)
  gp4 = graphon4(n,m)
  
  # SAS method
  Dat1[d,2] = mean((gp1$P-SASmooth(gp1$A,c(k,k,k)))^2)
  Dat2[d,2] = mean((gp2$P-SASmooth(gp2$A,c(k,k,k)))^2)
  Dat3[d,2] = mean((gp3$P-SASmooth(gp3$A,c(k,k,k)))^2)
  Dat4[d,2] = mean((gp4$P-SASmooth(gp4$A,c(k,k,k)))^2)
  ## Spectral method
  Dat1[d,3] = mean((gp1$P-Spectral(gp1$A,1,c(2,3)))^2)
  Dat2[d,3] = mean((gp2$P-Spectral(gp2$A,1,c(2,3)))^2)
  Dat3[d,3] = mean((gp3$P-Spectral(gp3$A,1,c(2,3)))^2)
  Dat4[d,3] = mean((gp4$P-Spectral(gp4$A,1,c(2,3)))^2)
  
  ## Hspectral method
  k = sqrt(n)
  Dat1[d,4] = mean((gp1$P-Hspectral(gp1$A,c(k,k,k)))^2)
  Dat2[d,4] = mean((gp2$P-Hspectral(gp2$A,c(k,k,k)))^2)
  Dat3[d,4] = mean((gp3$P-Hspectral(gp3$A,c(k,k,k)))^2)
  Dat4[d,4] = mean((gp4$P-Hspectral(gp4$A,c(k,k,k)))^2)
  
  
  ## SBM method
  ini=HSC(gp1$A,k,k,k,sym=T)
  res =tbmClustering(gp1$A,k,k,k,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = T, diagP = T)
  Dat1[d,5] = mean((res$judgeX-gp1$P)^2)
  
  ini=HSC(gp2$A,k,k,k,sym=T)
  res =tbmClustering(gp2$A,k,k,k,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = T, diagP = T)
  Dat2[d,5] = mean((res$judgeX-gp2$P)^2)
  
  ini=HSC(gp3$A,k,k,k,sym=T)
  res =tbmClustering(gp3$A,k,k,k,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = T, diagP = T)
  Dat3[d,5] = mean((res$judgeX-gp3$P)^2)
  
  ini=HSC(gp4$A,k,k,k,sym=T)
  res =tbmClustering(gp4$A,k,k,k,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = T, diagP = T)
  Dat4[d,5] = mean((res$judgeX-gp4$P)^2)
  
  Dat1[d,6] = mean((gp1$P-Hspectral2(gp1$A))^2)
  Dat2[d,6] = mean((gp2$P-Hspectral2(gp2$A))^2)
  Dat3[d,6] = mean((gp3$P-Hspectral2(gp3$A))^2)
  Dat4[d,6] = mean((gp4$P-Hspectral2(gp4$A))^2)
  
  print(paste(d,"-th iteration ended",sep = ""))
}




dat1 = dat2 = dat3 = dat4 = as.data.frame(matrix(nrow = 50,ncol = 4))
colnames(dat1) = colnames(dat2) = colnames(dat3) = colnames(dat4) = c("dim","method","MSE","graphon")
dat1[,1] = rep(Dat1[,1],5)
dat1[,2] = rep(c("SAS","Spectral","Hspectral","SBM","Hspectral2"),each = 10)
dat1[,3] = c(Dat1[,2],Dat1[,3],Dat1[,4],Dat1[,5],Dat1[,6])
dat1[,4] = "Graphon1"

dat2[,1] = rep(Dat2[,1],5)
dat2[,2] = rep(c("SAS","Spectral","Hspectral","SBM","Hspectral2"),each = 10)
dat2[,3] = c(Dat2[,2],Dat2[,3],Dat2[,4],Dat2[,5],Dat2[,6])
dat2[,4] = "Graphon2"

dat3[,1] = rep(Dat3[,1],5)
dat3[,2] = rep(c("SAS","Spectral","Hspectral","SBM","Hspectral2"),each = 10)
dat3[,3] = c(Dat3[,2],Dat3[,3],Dat3[,4],Dat3[,5],Dat3[,6])
dat3[,4] = "Graphon3"

dat4[,1] = rep(Dat4[,1],5)
dat4[,2] = rep(c("SAS","Spectral","Hspectral","SBM","Hspectral2"),each = 10)
dat4[,3] = c(Dat4[,2],Dat4[,3],Dat4[,4],Dat4[,5],Dat4[,6])
dat4[,4] = "Graphon4"


dat = rbind(dat1,dat2,dat3,dat4)

dat[,2] = as.factor(dat[,2])
dat[,4] = as.factor(dat[,4])


ggplot(data = dat,aes(x = dim, y = MSE, color = method))+geom_point(size = .5)+geom_line()+facet_wrap(~graphon,ncol = 2)
ggplot(data = dat[dat$method=="Hspectral"|dat$method=="Hspectral2",],aes(x = dim, y = MSE, color = method))+geom_point()+geom_line()+facet_wrap(~graphon,ncol = 2)


############################ grapon analysis ############################################

Dat1 = Dat2 = Dat3 = Dat4 = as.data.frame(matrix(nrow = 10,ncol = 4))
colnames(Dat1) = colnames(Dat2) = colnames(Dat3) =colnames(Dat4) = c("d","SAS","Spectral","HSpectral")
Dat1[,1]= Dat2[,1] = Dat3[,1] = Dat4[,1] = 10*(1:10)

for(d in 1:10){
  n = 10*d
  m = 2
  k = round((n)^{m/(m+2)})
  gp1 = graphon1(n,m)
  gp2 = graphon2(n,m)
  gp3 = graphon3(n,m)
  gp4 = graphon4(n,m)
  
  # SAS method
  Dat1[d,2] = mean((gp1$P-SASmooth(gp1$A,c(k,k)))^2)
  Dat2[d,2] = mean((gp2$P-SASmooth(gp2$A,c(k,k)))^2)
  Dat3[d,2] = mean((gp3$P-SASmooth(gp3$A,c(k,k)))^2)
  Dat4[d,2] = mean((gp4$P-SASmooth(gp4$A,c(k,k)))^2)
  ## Spectral method
  Dat1[d,3] = mean((gp1$P-Spectral(gp1$A,1,2))^2)
  Dat2[d,3] = mean((gp2$P-Spectral(gp2$A,1,2))^2)
  Dat3[d,3] = mean((gp3$P-Spectral(gp3$A,1,2))^2)
  Dat4[d,3] = mean((gp4$P-Spectral(gp4$A,1,2))^2)
  
  ## Hspectral method
  k = sqrt(n)
  Dat1[d,4] = mean((gp1$P-Hspectral(gp1$A,c(k,k)))^2)
  Dat2[d,4] = mean((gp2$P-Hspectral(gp2$A,c(k,k)))^2)
  Dat3[d,4] = mean((gp3$P-Hspectral(gp3$A,c(k,k)))^2)
  Dat4[d,4] = mean((gp4$P-Hspectral(gp4$A,c(k,k)))^2)
  
  print(paste(d,"-th iteration ended",sep = ""))
}




dat1 = dat2 = dat3 = dat4 = as.data.frame(matrix(nrow = 30,ncol = 4))
colnames(dat1) = colnames(dat2) = colnames(dat3) = colnames(dat4) = c("dim","method","MSE","graphon")
dat1[,1] = rep(Dat1[,1],3)
dat1[,2] = rep(c("SAS","Spectral","Hspectral"),each = 10)
dat1[,3] = c(Dat1[,2],Dat1[,3],Dat1[,4])
dat1[,4] = "Graphon1"

dat2[,1] = rep(Dat2[,1],3)
dat2[,2] = rep(c("SAS","Spectral","Hspectral"),each = 10)
dat2[,3] = c(Dat2[,2],Dat2[,3],Dat2[,4])
dat2[,4] = "Graphon2"

dat3[,1] = rep(Dat3[,1],3)
dat3[,2] = rep(c("SAS","Spectral","Hspectral"),each = 10)
dat3[,3] = c(Dat3[,2],Dat3[,3],Dat3[,4])
dat3[,4] = "Graphon3"

dat4[,1] = rep(Dat4[,1],3)
dat4[,2] = rep(c("SAS","Spectral","Hspectral"),each = 10)
dat4[,3] = c(Dat4[,2],Dat4[,3],Dat4[,4])
dat4[,4] = "Graphon4"


dat = rbind(dat1,dat2,dat3,dat4)

dat[,2] = as.factor(dat[,2])
dat[,4] = as.factor(dat[,4])


ggplot(data = dat,aes(x = dim, y = MSE, color = method))+geom_point()+geom_line()+facet_wrap(~graphon,ncol = 2)



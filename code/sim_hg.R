source("functions_hg.R")
source("functions_sbm.R")
library("RSKC")
set.seed(1)
g = 20;
result=result2=matrix(nrow=5,ncol=5);

n_list=c(50,100,150,200,250);


for(i in 1:5){
    for(rep in 1:5){
W = pbtensor(g,sym = T,smooth = T)
## results seem invariant to cluster permutations within each mode.
#shuffle=sample(1:g,g,replace=F)##
#W=W[shuffle,shuffle,shuffle]
## image(W[,,1])
hgmod = hgmodel.block(W,n_list[i],diagP = F,type="Gaussian")
k = r = l = g
A = hgmod$A
P = hgmod$P
Cs = hgmod$Cs

ini=HSC(A,k,r,l,sym=T)
res =tbmClustering(A,k,r,l,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = T, diagP = F)
result[i,rep]=mean((res$judgeX-P)^2)
result2[i,rep]=CER(hgmod$Cs,res$Cs)
}
}

plot(apply(result,1,mean))


source("functions_hg.R")
source("functions_sbm2.R")
library("RSKC")
set.seed(1)
g = 20;
result=result2=matrix(nrow=5,ncol=1);

n_list=c(50,100,150,200,250);

for(i in 1:5){
    
    W = signaltensor(g,10,sym = T,smooth = T)
    
    snmod = snmod.block(W,n_list[i],sym = F,diagP = T)
    k = r = l = g
    Sig = snmod$Sig
    obsTensor = snmod$obsTensor
    Cs = snmod$Cs
    
    ini=HSC(obsTensor,k,r,l,sym=F)
    res =tbmClustering(obsTensor,k,r,l,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = F, diagP = T)
    result[i]=mean((res$judgeX-Sig)^2)
    result2[i]=CER(snmod$Cs,res$Cs)
    
}

plot(apply(result,1,mean))

source("functions_hg.R")
source("functions_sbm.R")
set.seed(2)
g = 20;
result=result2=matrix(nrow=5,ncol=5);

n_list=c(50,100,150,200,250);
W = pbtensor(g,sym = T)

for(i in 1:5){
    for(rep in 1:5){
hgmod = hgmodel.block(W,n_list[i],diagP = F)
k = r = l = g
A = hgmod$A
P = hgmod$P
Cs = hgmod$Cs

ini=HSC(A,k,r,l)
res =tbmClustering(A,k,r,l,Cs.init=ini$Cs,Ds.init=ini$Cs,Es.init=ini$Cs,sym = T, diagP = F)
result[i,rep]=mean((res$judgeX-P)^2)

ini2=HSC(P,k,r,l)
res2 =tbmClustering(P,k,r,l,Cs.init=ini2$Cs,Ds.init=ini2$Cs,Es.init=ini2$Cs,sym = T, diagP = F)
result2[i,rep]=mean((res2$judgeX-P)^2)

}
}

plot(apply(result,1,mean)-apply(result2,1,mean))

args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])
} else {
  stop()
}
source("functions_sbm.R")
source("functions_hg.R")




index = matrix(nrow = 80,ncol = 2)
s = 0
for(i in 1:4){
  for(j in 1:20){
    s = s+1
    index[s,] = c(5*i,20*(j+1))
  }
}

g = index[BATCH,1]
n = index[BATCH,2]


err1 = NULL; err2 = NULL
for(r in 1:5){
      W = pbtensor(g)
      hgmod = hgmodel.block(W,n)
      A = hgmod$A
      P = hgmod$P
      res = tbmClustering(A,g,g,g,sym=T)
      err1 = c(err1,mean(abs(cut(res$judgeX)-P)^2))
      err2 = c(err2,sum(abs(cut(res$judgeX)-P)^2)/length(which(P!=0)))
}
MSE1 = mean(err1); MSE2 = mean(err2)
sd1 = sd(err1); sd2 = sd(err2)

save(MSE1,MSE2,sd1,sd2,file = paste("sbm",g,"_",n,".RData",sep = ""))

# smooth hypergraphon model
hgmod = hgmodel.smooth(n,1)
A = hgmod$A
P = hgmod$P

res = tbmClustering(A,g,g,g)
mean(abs(cut(res$judgeX)-P)^2)

err = NULL
for(r in 1:5){
  hgmod = hgmodel.smooth(n,1)
  A = hgmod$A
  P = hgmod$P
  
  res = tbmClustering(A,g,g,g)
  err = c(err,mean(abs(cut(res$judgeX)-P)^2))
}
MSE = mean(err)
sd = sd(err)
save(MSE,sd,file = paste("smooth",g,"_",n,".RData",sep = ""))



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


err = NULL; 
for(r in 1:5){
      W = pbtensor(g)
      W2 = pbtensor(g,sym = F)
      hgmod = hgmodel.block(W,n,diagP = F)

      A = hgmod$A
      P = hgmod$P
      res = tbmClustering(A,g,g,g,sym=T,diagP = F)

      
      
      err = c(err,mean(abs(res$judgeX-P)^2))

}

#A is generated from sym
#A2 is generated from non sym



MSE = mean(err)
sd = sd(err)
save(MSE,sd,file = paste("sbm_verify",g,"_",n,".RData",sep = ""))


# 
# ##################################################################################
# # smooth hypergraphon model
# hgmod = hgmodel.smooth(n,1)
# A = hgmod$A
# P = hgmod$P
# 
# res = tbmClustering(A,g,g,g)
# mean(abs(cut(res$judgeX)-P)^2)
# 
# err = NULL
# for(r in 1:5){
#   hgmod = hgmodel.smooth(n,1)
#   A = hgmod$A
#   P = hgmod$P
#   
#   res = tbmClustering(A,g,g,g)
#   err = c(err,mean(abs(cut(res$judgeX)-P)^2))
# }
# MSE = mean(err)
# sd = sd(err)
# save(MSE,sd,file = paste("smooth",g,"_",n,".RData",sep = ""))
# 
# library(plot.matrix)
# par(mar=c(5.1, 4.1, 4.1, 4.1))
# plot(A[,,2])
# plot(P[,,3])
# plot(W[,,1])

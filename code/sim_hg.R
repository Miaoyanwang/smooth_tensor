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



# comparison of two scenarios
err1 = NULL; err2= NULL
rep = 1
for(r in 1:rep){
  W = pbtensor(g)  # symmetric core => picking one non-diagonal entry
  W2 = pbtensor(g,sym = F) # nonsymmetric=> average corresponding non diagnal entries
  hgmod = hgmodel.block(W,n,diagP = F)
  hgmod2 = hgmodel.block(W2,n,diagP = F)
  
    
  A = hgmod$A
  P = hgmod$P
    
  A2 = hgmod2$A
  P2 = hgmod2$P
  
  res = tbmClustering(A,g,g,g,sym=T,diagP = F)
  res2 = tbmClustering(A2,g,g,g,sym=T,diagP = F)
  
  err1 = c(err1,mean(abs(res$judgeX-P)^2))
  err2 = c(err2,mean(abs(res2$judgeX-P2)^2))
}
  
  

#A is generated from sym
#A2 is generated from non ]-sym

mean(abs(res$judgeX-P)^2)
plot(c(res$judgeX[,,1]),c(P[,,1]))

mean(abs(res2$judgeX-P2)^2)
plot(c(res2$judgeX[,,1]),c(P2[,,1]))




# MSE = mean(err)
# sd = sd(err)
# save(MSE,sd,file = paste("sbm_verify",g,"_",n,".RData",sep = ""))


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

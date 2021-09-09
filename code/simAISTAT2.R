args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])
} else {
  stop()
}
source("functions_blk.R")
library(rTensor)



sim = 1:20
dim  = (1:10)*10

ind = matrix(nrow = 200,ncol = 2)
r = 0
for(i in sim){
  for(j in dim){
    r = r+1
    ind[r,] = c(i,j)
  }
}
s = ind[BATCH,1]
d = ind[BATCH,2]

set.seed(s)

result = as.data.frame(matrix(nrow =70 ,ncol = 6))
names(result) =  c("sim","dim","type","model","method","MSE")

result$sim = s
result$dim = d
result$type = rep(c("conti","binary"),each = 35)
result$model = rep(rep(c("model1","model2","model3","model4","model5"),each = 7),2)
result$method = rep(c("Borda0","Borda1","Borda2","Borda3","LSE","Bicluster","Spectral"),10)



MSE = NULL
for(m in 1:5){
  s1 = simulation(d,mode = m,signal_level = 5)
  k = ceiling(d^(3/5))
  MSE = c(MSE,mean((Borda2(s1$observe,0,k)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k)-s1$signal)^2))
  k = ceiling(d^(1/3))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k)-s1$signal)^2))
  k = ceiling(d^(3/11))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k)-s1$signal)^2))
  k = ceiling(d^(3/5))
  MSE = c(MSE,mean((LSE(s1$observe,k,mode = 3)-s1$signal)^2))
  MSE = c(MSE,mean((Bicluster(s1$observe,k)-s1$signal)^2))
  MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2))
  
}
for(m in 1:5){
  s1 = simulation_bin(d,mode = m)
  k = ceiling(d^(3/5))
  MSE = c(MSE,mean((Borda2(s1$observe,0,k)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k)-s1$signal)^2))
  k = ceiling(d^(1/3))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k)-s1$signal)^2))
  k = ceiling(d^(3/11))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k)-s1$signal)^2))
  k = ceiling(d^(3/5))
  MSE = c(MSE,mean((LSE(s1$observe,k,mode = 3)-s1$signal)^2))
  MSE = c(MSE,mean((Bicluster(s1$observe,k)-s1$signal)^2))
  MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2))
}


result$MSE = MSE

save(result,file = paste("nsim_",s,"_dim_",d,".RData",sep =""))


# 
# setwd("/Users/chanwoolee/Desktop/Research projects/Hypergraphon/code/AISTATSIM")
# tot= NULL
# for(BATCH in c(1:30,33:45,47:53,55:61,63:69,71:80,82:86,88:94,96,98:100)){
#   s = ind[BATCH,1]
#   m = ind[BATCH,2]
#   load( paste("model",m,"sim",s,".RData",sep =""))
#   tot = rbind(tot,result)   
# }
# 
# 
# tot$model = as.factor(tot$model)
# tot$method = as.factor(tot$method)
# 
# 
# totsummary = summarySE(tot, measurevar="MSE", groupvars=c("model","method","dim"))
# totsummary
# 
# g1 = ggplot(totsummary[totsummary$model=="model1",], aes(x=dim, y=MSE,colour = method)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
# g2 = ggplot(totsummary[totsummary$model=="model2",], aes(x=dim, y=MSE,colour = method)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
# g3 = ggplot(totsummary[totsummary$model=="model3",], aes(x=dim, y=MSE,colour = method)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
# g4 = ggplot(totsummary[totsummary$model=="model4",], aes(x=dim, y=MSE,colour = method)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
# g5 = ggplot(totsummary[totsummary$model=="model5",], aes(x=dim, y=MSE,colour = method)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
# 
# library(ggpubr)
# ggarrange(g1,g2,g3,g4,g5,nrow = 2,ncol = 3,labels = c("model 1","model 2","model 3","model 4","model 5"))









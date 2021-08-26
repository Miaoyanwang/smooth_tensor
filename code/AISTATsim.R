source("functions_blk.R")
### Tensor visualization #########################################
source("tensor_visualization.R")
s1 = simulation(30,mode = 1)
plot_tensor(s1$signal)
plot_tensor(s1$observe)

s2 = simulation(30,mode = 2)
plot_tensor(s2$signal)
plot_tensor(s2$observe)

s3 = simulation(30,mode = 3)
plot_tensor(s3$signal)
plot_tensor(s3$observe)


s4 = simulation(30,mode = 4)
plot_tensor(s4$signal)
plot_tensor(s4$observe)


s5 = simulation(30,mode = 5)
plot_tensor(s5$signal)
plot_tensor(s5$observe)

l=3; k=floor(d^(1/3))
est=polytensor(data,l,k)
plot(est,tensor)
mean((est-tensor)^2)


##############################################################################
################## comparison between polynomial vs piecewise constant
poly_vs_const = as.data.frame(matrix(nrow = 600,ncol = 5))
names(poly_vs_const) = c("dim","sim","model","deg","MSE")
dim = c(30,60,90)
poly_vs_const$dim = rep(dim,each = 200)
poly_vs_const$sim = rep(rep(1:10,each = 20),3)
poly_vs_const$model = rep(rep(c("model1","model2","model3","model4","model5"),each = 4),30)
poly_vs_const$deg = rep(0:3,150)

MSE = NULL
for(d in dim){
  k = ceiling(d^(1/3))
  for(s in 1:10){
    set.seed(s)
    for(m in 1:5){
      s1 = simulation(d,mode = m)
      MSE = c(MSE,mean((Borda2(s1$observe,0,k)-s1$signal)^2))
      MSE = c(MSE,mean((Borda2(s1$observe,1,k)-s1$signal)^2))
      MSE = c(MSE,mean((Borda2(s1$observe,2,k)-s1$signal)^2))
      MSE = c(MSE,mean((Borda2(s1$observe,3,k)-s1$signal)^2))
      print(paste(d,"_",s,"_",m,"_ended",sep = ""))
    }
    
  }
}
poly_vs_const$MSE = MSE
poly_vs_const$deg = as.factor(poly_vs_const$deg)

datsummary = summarySE(poly_vs_const, measurevar="MSE", groupvars=c("model","deg","dim"))

library(ggplot2)
ggplot(datsummary, aes(x=dim, y=MSE, colour=deg)) + facet_wrap(~model, nrow=1) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)



#save(datsummary50,poly_vs_const50,datsummary100,poly_vs_const100,file = "poly_vs_const.RData")


##### Comparison with other methods #####################################

alternative = as.data.frame(matrix(nrow = 600,ncol = 5))
names(alternative) = c("dim","sim","model","method","MSE")
dim = c(20,40,60,80)
alternative$dim = rep(dim,each = 150)
alternative$sim = rep(rep(1:10,each = 15),4)
alternative$model = rep(rep(c("model1","model2","model3","model4","model5"),each = 3),40)
alternative$method = rep(c("Borda","LSE","Spectral"),200)

MSE = NULL
for(d in dim){
  k = ceiling(d^(1/3))
  for(s in 1:10){
    set.seed(s)
    for(m in 1:5){
      s1 = simulation_bin(d,mode = m)
      MSE = c(MSE,mean((Borda2(s1$observe,3,k)-s1$signal)^2))
      MSE = c(MSE,mean((LSE(s1$observe,k)-s1$signal)^2))
      MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2))
      print(paste(d,"_",s,"_",m,"_ended",sep = ""))
    }
    
  }
}

alternative$MSE = MSE

alternative$model = as.factor(alternative$model)
alternative$method = as.factor(alternative$method)


datsummary = summarySE(alternative, measurevar="MSE", groupvars=c("model","method","dim"))


library(ggplot2)
ggplot(datsummary, aes(x=dim, y=MSE, colour=method)) + facet_wrap(~model, nrow=1) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)







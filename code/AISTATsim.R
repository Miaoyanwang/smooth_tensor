source("functions_blk.R")
### Tensor visualization #########################################
##### signal tensor estimation ###################################
source("tensor_visualization.R")
s1 = simulation(30,mode = 1)
plot_tensor(s1$signal)
plot_tensor(s1$observe)
k = ceiling(30^(1/3))
plot_tensor(Borda2(s1$observe,2,k));mean((Borda2(s1$observe,2,k)-s1$signal)^2)
plot_tensor(Spectral(s1$observe,1,c(2,3)));mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2)
plot_tensor(tucker(as.tensor(s1$observe),c(k,k,k))$est@data);mean((tucker(as.tensor(s1$observe),c(k,k,k))$est@data-s1$signal)^2)



s2 = simulation(30,mode = 2)
plot_tensor(s2$signal)
plot_tensor(s2$observe)
plot_tensor(Borda2(s2$observe,2,k));mean((Borda2(s2$observe,2,k)-s2$signal)^2)
plot_tensor(Spectral(s2$observe,1,c(2,3)));mean((Spectral(s2$observe,1,c(2,3))-s2$signal)^2)
plot_tensor(tucker(as.tensor(s2$observe),c(k,k,k))$est@data);mean((tucker(as.tensor(s2$observe),c(k,k,k))$est@data-s2$signal)^2)


s3 = simulation(30,mode = 3)
plot_tensor(s3$signal)
plot_tensor(s3$observe)
plot_tensor(Borda2(s3$observe,2,k));mean((Borda2(s3$observe,2,k)-s3$signal)^2)
plot_tensor(Spectral(s3$observe,1,c(2,3)));mean((Spectral(s3$observe,1,c(2,3))-s3$signal)^2)
plot_tensor(tucker(as.tensor(s3$observe),c(k,k,k))$est@data);mean((tucker(as.tensor(s3$observe),c(k,k,k))$est@data-s3$signal)^2)


s4 = simulation(30,mode = 4)
plot_tensor(s4$signal)
plot_tensor(s4$observe)
plot_tensor(Borda2(s4$observe,2,k));mean((Borda2(s4$observe,2,k)-s4$signal)^2)
plot_tensor(Spectral(s4$observe,1,c(2,3)));mean((Spectral(s4$observe,1,c(2,3))-s4$signal)^2)
plot_tensor(tucker(as.tensor(s4$observe),c(k,k,k))$est@data);mean((tucker(as.tensor(s4$observe),c(k,k,k))$est@data-s4$signal)^2)



s5 = simulation(30,mode = 5)
plot_tensor(s5$signal)
plot_tensor(s5$observe)
plot_tensor(Borda2(s5$observe,2,k));mean((Borda2(s5$observe,2,k)-s5$signal)^2)
plot_tensor(Spectral(s5$observe,1,c(2,3)));mean((Spectral(s5$observe,1,c(2,3))-s5$signal)^2)
plot_tensor(tucker(as.tensor(s5$observe),c(k,k,k))$est@data);mean((tucker(as.tensor(s5$observe),c(k,k,k))$est@data-s5$signal)^2)


##### Probability tensor estimation visualization ################################
s1 = simulation_bin(30,mode = 1)
plot_tensor(s1$signal)
plot_tensor(s1$observe)
k = ceiling(30^(1/3))
plot_tensor(Borda2(s1$observe,2,k));mean((Borda2(s1$observe,2,k)-s1$signal)^2)
plot_tensor(Spectral(s1$observe,1,c(2,3)));mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2)
k = ceiling(0.6*30^(3/5))
plot_tensor(LSE(s1$observe,k,mode = 1));mean((LSE(s1$observe,k,mode = 1)-s1$signal)^2)
plot_tensor(LSE(s1$observe,k,mode = 2));mean((LSE(s1$observe,k,mode = 2)-s1$signal)^2)
plot_tensor(LSE(s1$observe,k,mode = 3));mean((LSE(s1$observe,k,mode = 3)-s1$signal)^2)



s2 = simulation_bin(30,mode = 2)
plot_tensor(s2$signal)
plot_tensor(s2$observe)
plot_tensor(Borda2(s2$observe,2,k));mean((Borda2(s2$observe,2,k)-s2$signal)^2)
plot_tensor(Spectral(s2$observe,1,c(2,3)));mean((Spectral(s2$observe,1,c(2,3))-s2$signal)^2)
plot_tensor(LSE(s2$observe,k));mean((LSE(s2$observe,k,rep = 100)-s2$signal)^2)


s3 = simulation_bin(30,mode = 3)
plot_tensor(s3$signal)
plot_tensor(s3$observe)
plot_tensor(Borda2(s3$observe,2,k));mean((Borda2(s3$observe,2,k)-s3$signal)^2)
plot_tensor(Spectral(s3$observe,1,c(2,3)));mean((Spectral(s3$observe,1,c(2,3))-s3$signal)^2)
plot_tensor(LSE(s3$observe,k));mean((LSE(s3$observe,k,rep = 10)-s3$signal)^2)

s4 = simulation_bin(30,mode = 4)
plot_tensor(s4$signal)
plot_tensor(s4$observe)
plot_tensor(Borda2(s4$observe,2,k));mean((Borda2(s4$observe,2,k)-s4$signal)^2)
plot_tensor(Spectral(s4$observe,1,c(2,3)));mean((Spectral(s4$observe,1,c(2,3))-s4$signal)^2)
plot_tensor(LSE(s4$observe,k));mean((LSE(s4$observe,k,rep = 10)-s4$signal)^2)


s5 = simulation_bin(30,mode = 5)
plot_tensor(s5$signal)
plot_tensor(s5$observe)
plot_tensor(Borda2(s5$observe,2,k));mean((Borda2(s5$observe,2,k)-s5$signal)^2)
plot_tensor(Spectral(s5$observe,1,c(2,3)));mean((Spectral(s5$observe,1,c(2,3))-s5$signal)^2)
plot_tensor(LSE(s5$observe,k));mean((LSE(s5$observe,k,rep = 10)-s5$signal)^2)



############## Reproducing the hypergraphon estimaiton paper ####################


RMSE1 = RMSE2 = RMSE3 = NULL
for( i in 10*(1:8)){
  s1 = simulation_bin(i,mode = 1)
  k = ceiling(0.6*i^(3/5))
  
  RMSE1 = c(RMSE1,sum((LSE(s1$observe,k,mode = 1)-s1$signal)^2)/sum(s1$signal^2))
  RMSE2 = c(RMSE2,sum((LSE(s1$observe,k,mode = 2)-s1$signal)^2)/sum(s1$signal^2))
  RMSE3 = c(RMSE3,sum((LSE(s1$observe,k,mode = 3)-s1$signal)^2)/sum(s1$signal^2))
  
  print(paste(i,"_ended",sep = ""))
}  
plot(10*(0:8),c(1,RMSE1),xlab = "Number of nodes (n)",ylab = "Normalized Reconstruction Error",type = "b",ylim = c(0,1))
points(10*(0:8),c(1,RMSE2),type = "b",col = "red")
points(10*(0:8),c(1,RMSE3),type = "b",col = "blue")
legend("topright",legend = c("MSC","BAL","HSC"),
       col = c("black","red","blue"),lty = 1)



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


##### Comparison with other methods (Binary) #####################################

alternative = as.data.frame(matrix(nrow = 600,ncol = 5))
names(alternative) = c("dim","sim","model","method","MSE")
dim = c(20,40,60,80)
alternative$dim = rep(dim,each = 150)
alternative$sim = rep(rep(1:10,each = 15),4)
alternative$model = rep(rep(c("model1","model2","model3","model4","model5"),each = 3),40)
alternative$method = rep(c("Borda","LSE","Spectral"),200)

MSE = NULL
for(d in dim){
  for(s in 1:10){
    set.seed(s)
    for(m in 1:5){
      s1 = simulation_bin(d,mode = m)
      k = ceiling(d^(1/3))
      MSE = c(MSE,mean((Borda2(s1$observe,2,k)-s1$signal)^2))
      k = ceiling(d^(3/5))
      MSE = c(MSE,mean((LSE(s1$observe,k,mode = 3)-s1$signal)^2))
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
ggplot(datsummary[datsummary$model =="model5",], aes(x=dim, y=MSE, colour=method)) + facet_wrap(~model, nrow=1) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)


##### Comparison with other methods (continuous) #####################################

alternativec = as.data.frame(matrix(nrow = 1000,ncol = 5))
names(alternativec) = c("dim","sim","model","method","MSE")
dim = c(20,40,60,80)
alternativec$dim = rep(dim,each = 250)
alternativec$sim = rep(rep(1:10,each = 25),4)
alternativec$model = rep(rep(c("model1","model2","model3","model4","model5"),each = 5),40)
alternativec$method = rep(c("Borda","Tucker1","Spectral","LSE","Tucker2"),200)


MSE = NULL
for(d in dim){
  for(s in 1:10){
    set.seed(s)
    for(m in 1:5){
      s1 = simulation(d,mode = m)
      k = ceiling(d^(1/3))
      MSE = c(MSE,mean((Borda2(s1$observe,2,k)-s1$signal)^2))
      MSE = c(MSE,mean((tucker(as.tensor(s1$observe),c(k,k,k))$est@data-s1$signal)^2))
      MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2))
      k = ceiling(d^(3/5))
      MSE = c(MSE,mean((LSE(s1$observe,k,mode = 3)-s1$signal)^2))
      MSE = c(MSE,mean((tucker(as.tensor(s1$observe),c(k,k,k))$est@data-s1$signal)^2))
      print(paste(d,"_",s,"_",m,"_ended",sep = ""))
    }
    
  }
}

alternativec$MSE = MSE

alternativec$model = as.factor(alternativec$model)
alternativec$method = as.factor(alternativec$method)


datsummaryc = summarySE(alternativec, measurevar="MSE", groupvars=c("model","method","dim"))


library(ggplot2)
ggplot(datsummaryc[datsummaryc$model =="model5",], aes(x=dim, y=MSE, colour=method)) + facet_wrap(~model, nrow=1) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)



# signal
# est = tucker(as.tensor(s1$observe),c(k,k,k))
# mean((tucker(as.tensor(s1$observe),c(k,k,k))$est@data-s1$signal)^2)




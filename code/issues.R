############## Replicating the hypergraphon estimaiton paper ####################


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



########### Producing missing imputation simulation ##############################
#### Continuous case #################################
missing = as.data.frame(matrix(nrow = 1600,ncol = 5))
names(missing) = c("fraction","sim","model","method","MSE")
fraction = (3:10)*0.1
missing$fraction = rep(fraction,each = 200)
missing$sim = rep(rep(1:10,each = 10),8)
missing$model = rep(rep(c("model1","model2","model3","model4","model5"),each = 2),160)
missing$method = rep(c("Borda","Spectral"),800)


MSE = NULL
for(rho in fraction){
  for(s in 1:20){
    set.seed(s)
    for(m in 1:5){
      d = 50
      s1 = simulation(d,mode = m,signal_level = 5)
      A = s1$observe
      obs=rbinom(d^3,1,rho)
      A[obs ==0] = NA
      k = ceiling(d^(1/3))
      MSE = c(MSE,mean((Borda2(A,2,k)-s1$signal)^2))
      MSE = c(MSE,mean((Spectral(A,1,c(2,3),lambda = d)-s1$signal)^2))
      print(paste(rho,"_",s,"_",m,"_ended",sep = ""))
    }
    
  }
}

missing$MSE = MSE

missing$model = as.factor(missing$model)
missing$method = as.factor(missing$method)

save(missing,file = "conti_missing.RData")


######## Binary  case #################################
missing2 = as.data.frame(matrix(nrow = 1600,ncol = 5))
names(missing2) = c("fraction","sim","model","method","MSE")
fraction = (3:10)*0.1
missing2$fraction = rep(fraction,each = 200)
missing2$sim = rep(rep(1:10,each = 10),8)
missing2$model = rep(rep(c("model1","model2","model3","model4","model5"),each = 2),160)
missing2$method = rep(c("Borda","Spectral"),800)


MSE = NULL
for(rho in fraction){
  for(s in 1:20){
    set.seed(s)
    for(m in 1:5){
      d = 50
      s1 = simulation_bin(d,mode = m)
      A = s1$observe
      obs=rbinom(d^3,1,rho)
      A[obs ==0] = NA
      k = ceiling(d^(1/3))
      MSE = c(MSE,mean((Borda2(A,2,k)-s1$signal)^2))
      MSE = c(MSE,mean((Spectral(A,1,c(2,3),lambda = 10)-s1$signal)^2))
      print(paste(rho,"_",s,"_",m,"_ended",sep = ""))
    }
    
  }
}

missing2$MSE = MSE

missing2$model = as.factor(missing2$model)
missing2$method = as.factor(missing2$method)


save(missing2,file = "binary_missing.RData")
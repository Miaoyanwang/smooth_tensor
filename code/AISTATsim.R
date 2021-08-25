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
poly_vs_const = matrix(nrow = 400,ncol = 5)
d = 150
k1 = floor(d^(1/3))
k2 = floor(d^(3/5))



poly_vs_const[,1] = trials = rep(1:10,each = 40)
poly_vs_const[,2] = model = rep(rep(c("sim1","sim2","sim3","sim4","sim5"),each = 8),10)
poly_vs_const[,3] = degree = rep(rep(0:3, each = 2),50)
poly_vs_const[,4] =  cluster = rep(c(k1,k2),200)


MSE = NULL
for(sim in 1:10){
  set.seed(sim)
  s1 = simulation(d,mode = 1)
  MSE = c(MSE,mean((Borda(s1$observe,c(k1,k1,k1))-s1$signal)^2))
  MSE = c(MSE,mean((Borda(s1$observe,c(k2,k2,k2))-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k2)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k2)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k2)-s1$signal)^2))
  
  
  
  s1 = simulation(d,mode = 2)
  MSE = c(MSE,mean((Borda(s1$observe,c(k1,k1,k1))-s1$signal)^2))
  MSE = c(MSE,mean((Borda(s1$observe,c(k2,k2,k2))-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k2)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k2)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k2)-s1$signal)^2))
  
  
  s1 = simulation(d,mode = 3)
  MSE = c(MSE,mean((Borda(s1$observe,c(k1,k1,k1))-s1$signal)^2))
  MSE = c(MSE,mean((Borda(s1$observe,c(k2,k2,k2))-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k2)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k2)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k2)-s1$signal)^2))
  
  s1 = simulation(d,mode = 4)
  MSE = c(MSE,mean((Borda(s1$observe,c(k1,k1,k1))-s1$signal)^2))
  MSE = c(MSE,mean((Borda(s1$observe,c(k2,k2,k2))-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k2)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k2)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k2)-s1$signal)^2))
  
  s1 = simulation(d,mode = 5)
  MSE = c(MSE,mean((Borda(s1$observe,c(k1,k1,k1))-s1$signal)^2))
  MSE = c(MSE,mean((Borda(s1$observe,c(k2,k2,k2))-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k2)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k2)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k1)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k2)-s1$signal)^2))
  
  print(paste(sim,"-th sim ended",sep = ""))
}

poly_vs_const[,5] =  as.numeric(MSE)
colnames(poly_vs_const) = c("trials","model","degree","cluster","MSE")
dat = as.data.frame(poly_vs_const)

dat$model = as.factor(dat$model)
dat$cluster = as.factor(dat$cluster)
dat$MSE = as.numeric(dat$MSE)
dat$degree = as.numeric(dat$degree)

datsummary = summarySE(dat, measurevar="MSE", groupvars=c("model","degree","cluster"))

library(ggplot2)
ggplot(datsummary, aes(x=degree, y=MSE, colour=cluster)) + facet_wrap(~model, nrow=1) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)



#save(datsummary50,poly_vs_const50,datsummary100,poly_vs_const100,file = "poly_vs_const.RData")


##### Comparison with other methods #####################################

alternative100
altdatsummary100

save(altdatsummary50,alternative50,altdatsummary100,alternative100,file = "alternative.RData")


alternative = matrix(nrow = 200,ncol = 5)
d = 50
k1 = floor(d^(1/3))
k2 = floor(d^(3/5))



alternative[,1]  = trials = rep(1:10,each = 20)
alternative[,2] = model = rep(rep(c("sim1","sim2","sim3","sim4","sim5"),each = 4),10)
alternative[,3] = degree = rep(rep(c("spectral","LSE"), each = 2),50)
alternative[,4] =  cluster = rep(c(k1,k2),100)

MSE = NULL
for(sim in 1:10){
  set.seed(sim)
  s1 = simulation(d,mode = 1)
  MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2),NA)
  MSE = c(MSE,mean((LSE(s1$observe,k1)-s1$signal)^2))
  MSE = c(MSE,mean((LSE(s1$observe,k2)-s1$signal)^2))
 
  
  
  s1 = simulation(d,mode = 2)
  MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2),NA)
  MSE = c(MSE,mean((LSE(s1$observe,k1)-s1$signal)^2))
  MSE = c(MSE,mean((LSE(s1$observe,k2)-s1$signal)^2))
  
  
  
  s1 = simulation(d,mode = 3)
  MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2),NA)
  MSE = c(MSE,mean((LSE(s1$observe,k1)-s1$signal)^2))
  MSE = c(MSE,mean((LSE(s1$observe,k2)-s1$signal)^2))
  
  
  s1 = simulation(d,mode = 4)
  MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2),NA)
  MSE = c(MSE,mean((LSE(s1$observe,k1)-s1$signal)^2))
  MSE = c(MSE,mean((LSE(s1$observe,k2)-s1$signal)^2))
  
  
  s1 = simulation(d,mode = 5)
  MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2),NA)
  MSE = c(MSE,mean((LSE(s1$observe,k1)-s1$signal)^2))
  MSE = c(MSE,mean((LSE(s1$observe,k2)-s1$signal)^2))
  
  
  print(paste(sim,"-th sim ended",sep = ""))
}


alternative[,5] =  as.numeric(MSE)
colnames(alternative) = c("trials","model","method","cluster","MSE")
dat = as.data.frame(alternative)

dat[,2] = as.factor(dat[,2])
dat[,3] = as.factor(dat[,3])
dat[,4] = as.factor(dat[,4])
dat$MSE = as.numeric(dat$MSE)


datsummary = summarySE(dat, measurevar="MSE", groupvars=c("model","method","cluster"))
alternative50 = na.omit(dat)
altdatsummary50 = na.omit(datsummary)

library(ggplot2)
ggplot(datsummary50, aes(x=degree, y=MSE, colour=cluster)) + facet_wrap(~model, nrow=1) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
ggplot(datsummary100, aes(x=degree, y=MSE, colour=cluster)) + facet_wrap(~model, nrow=1) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)

altdatsummary100$method = rep(c("LSE_k2","LSE_k1","Spectral"),5)
altdatsummary50$method = rep(c("LSE_k2","LSE_k1","Spectral"),5)

ind = c(1,2,7,8,9,10,15,16,17,18,23,24,25,26,31,32,33,34,39,40)
dsummary100 = datsummary100[ind,]
dsummary50 = datsummary50[ind,]
dsummary100

dsummary100$method = rep(c("Bordac_k2","Bordac_k1","Bordap_k2","Bordap_k1"),5)
dsummary50$method = rep(c("Bordac_k2","Bordac_k1","Bordap_k2","Bordap_k1"),5)


tot100 = rbind(altdatsummary100[,c(1,2,5,7)],dsummary100[,c(1,5,7,9)])
tot100$method =as.factor(tot100$method)
tot50 = rbind(altdatsummary50[,c(1,2,5,7)],dsummary50[,c(1,5,7,9)])
tot50$method =as.factor(tot50$method)
ggplot(tot50, aes(x=method, y=MSE)) + facet_wrap(~model, nrow=1) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)


tot50$dim = 50
tot100$dim = 100
tot100
tot = rbind(tot50,tot100)
tot$method = as.factor(tot$method)
ind1 = tot$method == "LSE_k2"|tot$method == "Spectral"|tot$method == "Bordac_k2"|tot$method == "Bordap_k2"
ind2 = tot$method == "LSE_k1"|tot$method == "Spectral"|tot$method == "Bordac_k1"|tot$method == "Bordap_k1"

si5 = tot$model=="sim5"
ggplot(tot[si5,], aes(x=dim, y=MSE)) + facet_wrap(~method, nrow=1) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)




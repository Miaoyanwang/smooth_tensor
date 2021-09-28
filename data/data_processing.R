setwd("../data/")
data=read.table("chicago-crime-comm.tns")
date_map=read.table("mode-1-date.map")
hour_map=read.csv("mode-2-hour.map",head = F)
area_map=read.csv("mode-3-communityarea.map",head = F)
crimetype_map=read.csv("mode-4-crimetype.map",head = F)
# reshape dataframe to tensor
names(data) = c("date","hour","area","type","count")
ndata = aggregate(data[5],data[2:4],sum)


library(reshape2)
tns = acast(ndata,hour~area~type,value.var = "count")
tns[is.na(tns)] = 0
hist(c(tns),breaks = 100)
# Take logarithm for the analysis
ltns = log(1+tns)
hist(c(ltns))
plot_tensor(ltns)
result = Borda2_asym(ltns,2,c(6,4,10))
plot_tensor2(ltns)
rgl.postscript('original.pdf', fmt = 'pdf')

plot_tensor2(ltns[result$order[[1]],result$order[[2]],result$order[[3]]])
result$Theta[result$Theta<0]=0
plot_tensor2(result$Theta[result$order[[1]],result$order[[2]],result$order[[3]]]);mean((result$Theta-ltns)^2)
rgl.postscript('denoisded.pdf', fmt = 'pdf')

result3 = Borda2_asym(ltns,0,c(6,4,10))
plot_tensor2(result3$Theta[result3$order[[1]],result3$order[[2]],result3$order[[3]]]);mean((result3$Theta-ltns)^2)
rgl.postscript('block.pdf', fmt = 'pdf')


save(tns,ltns,hour_map,area_map,crimetype_map,file = "Chicago.RData")


####### F statistics
source("functions_asym.R")
load("Chicago.RData")

result = Borda2_asym(ltns,2,c(6,4,10))
result3 = Borda2_asym(ltns,0,c(6,4,10))
kvec = c(6,4,10)
d = dim(ltns)
df1 = 9*prod(kvec)
df2 = prod(d)-10*prod(kvec)
Fstat = ((sum((result3$Theta-ltns)^2)-sum((result$Theta-ltns)^2))/df1)/(sum((result$Theta-ltns)^2)/df2)
pf(Fstat, df1, df2, lower.tail = TRUE, log.p = FALSE)




# finding correlations

ftor = array(1:prod(c(1,4,10)),dim = c(1,4,10))
fmat = ftor[rep(1,d[1]),z2,z3]
Avg_mat = array(aggregate(c(ltns[o1,o2,o3]),by = list(as.factor(fmat)),function(x) mean(x,na.rm = T))[,-1],dim = c(1,4,10))









plot_tensor2(result$Theta)
result2 = Spectral(ltns,2,c(1,3),lambda = 40)
plot_tensor(result2[result$order[[1]],result$order[[2]],result$order[[3]]]);mean((result2-ltns)^2)

result3 = Borda2_asym(ltns,0,c(7,18,10))
plot_tensor(result3$Theta[result$order[[1]],result3$order[[2]],result3$order[[3]]]);mean((result3$Theta-ltns)^2)

d = dim(ltns)
denoised = result$Theta
denoised[denoised<0]=0
denoised_sorted = denoised[1:d[1],result$order[[2]],result$order[[3]]]

library(plot.matrix)
cM = array(dim = c(d[1],4,d[3]))
for(i in 1:4){
  temp = matrix(0,nrow = d[1],ncol = d[3])
  for(j in which(z2==i)){
    temp = temp + denoised_sorted[,j,] 
  }
  cM[,i,] = temp/length(which(z2==i))
}

par(mar = c(8, 6, 2, 4))

plot(cM[,1,],axis.col=NULL, axis.row=NULL, xlab='', ylab='Hour',main = "Crime Area 1", 
     breaks=0:8,col =  brewer.pal(8, "OrRd"),cex.lab = 1.3,cex.main = 1.5,cex.axis = 1.2)
text(x = 1:32,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c(crimetype_map[o3,]),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 45,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 0.6)
axis(side = 2,
     ## Rotate labels perpendicular to y-axis.
     las = 2,
     ## Adjust y-axis label positions.
     at = c(4,8,12,16,20,24),
     mgp = c(3, 0.75, 0),cex = 1.7)
mtext("Crime Type", side=1, line=5,cex = 1.3)




plot(denoised_sorted[,7,],xlab = "Crime type",ylab = "Hour")
plot(denoised_sorted[,27,],xlab = "Crime type",ylab = "Hour")
plot(denoised_sorted[,57,],xlab = "Crime type",ylab = "Hour")
plot(denoised_sorted[,77,],xlab = "Crime type",ylab = "Hour")



order(sapply(1:d[1], function(x) sum(denoised_sorted[x,,],na.rm = T)))
plot(1:d[1],sapply(1:d[1], function(x) sum(denoised_sorted[x,,],na.rm = T)),type = "l",xlab = "indices",ylab = "")
order(sapply(1:d[2], function(x) sum(denoised_sorted[,x,],na.rm = T)))
plot(1:d[2],sapply(1:d[2], function(x) sum(denoised_sorted[,x,],na.rm = T)),type = "l",xlab = "indices",ylab = "")
order(sapply(1:d[3], function(x) sum(denoised_sorted[,,x],na.rm = T)))
plot(1:d[3],sapply(1:d[3], function(x) sum(denoised_sorted[,,x],na.rm = T)),type = "l",xlab = "indices",ylab = "")



plot(1:d[1],sapply(1:d[1], function(x) sum(ltns[o1,o2,o3][x,,],na.rm = T)),type = "l",xlab = "indices",ylab = "")
plot(1:d[2],sapply(1:d[2], function(x) sum(ltns[o1,o2,o3][,x,],na.rm = T)),type = "l",xlab = "indices",ylab = "")
plot(1:d[3],sapply(1:d[3], function(x) sum(ltns[o1,o2,o3][,,x],na.rm = T)),type = "l",xlab = "indices",ylab = "")



sort(hour_map[o1[which(z1==1)],])
sort(hour_map[o1[which(z1==2)],])
sort(hour_map[o1[which(z1==3)],])
sort(hour_map[o1[which(z1==4)],])
sort(hour_map[o1[which(z1==5)],])
sort(hour_map[o1[which(z1==6)],])

sort(o2[which(z2 == 1)])
sort(o2[which(z2 == 2)])
sort(o2[which(z2 == 3)])
sort(o2[which(z2 == 4)])

crimetype_map[o3[which(z3==1)],]
crimetype_map[o3[which(z3==2)],]
crimetype_map[o3[which(z3==3)],]
crimetype_map[o3[which(z3==4)],]
crimetype_map[o3[which(z3==5)],]
crimetype_map[o3[which(z3==6)],]
crimetype_map[o3[which(z3==7)],]
crimetype_map[o3[which(z3==8)],]
crimetype_map[o3[which(z3==9)],]
crimetype_map[o3[which(z3==10)],]






## hyper parameter setting

# CV indices
test = list()
tindex = sample(1:length(ltns),length(ltns))
q = length(tns)%/%5
for(i in 1:5){
  test[[i]] = tindex[(q*(i-1)+1):(q*i)]
}


test1 = test2 = train1 = train2  = NULL
for(i in 1:5){
  test_index = test[[i]]
  train_index = setdiff(1:length(ltns),test_index)
  train_tensor = ltns
  train_tensor[test_index] = NA
  result = Borda2_asym(train_tensor,2,c(6,4,10))
  train1 = c(train1,mean((result$Theta[train_index]-ltns[train_index])^2))
  test1 = c(test1,mean((result$Theta[test_index]-ltns[test_index])^2))

  
  result2 = Spectral(train_tensor,2,c(1,3),lambda = 50)
  
  train2 = c(train2,mean((result2[train_index]-ltns[train_index])^2))
  test2 = c(test2,mean((result2[test_index]-ltns[test_index])^2))
}


result = Borda2_asym(train_tensor,0,c(6,4,10)) #give us the best result

mean((result$Theta[test_index]-ltns[test_index])^2)


result2 = Spectral(train_tensor,1,c(2,3),lambda = 39)
mean((result2[train_index]-ltns[train_index])^2)
mean((result2[test_index]-ltns[test_index])^2)

result2 = Spectral(train_tensor,2,c(1,3),lambda = 6)
mean((result2[train_index]-ltns[train_index])^2)
mean((result2[test_index]-ltns[test_index])^2)

result2 = Spectral(train_tensor,3,c(1,2),lambda = 5)
mean((result2[train_index]-ltns[train_index])^2)
mean((result2[test_index]-ltns[test_index])^2)




result = Borda2_asym(ltns,2,c(2,3,1))
mean((result$Theta-ltns)^2)
result2 = Spectral(ltns,1,c(2,3),lambda = 16)
mean((result2-ltns)^2)


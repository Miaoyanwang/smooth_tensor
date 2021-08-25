library(rgl)
library("RColorBrewer")
marker = list(color = brewer.pal(9, "RdBu"))$color

### function to plot 3d array
plot_tensor=function(tensor){
  
  n=prod(dim(tensor))   
  #color_choice=min(round(prod(dim(tensor))/(6*6*6)),100)+1
  #marker=viridis_pal(option = "B")(color_choice)
  color_choice=round(prod(dim(tensor))/(20*20*20))+1
  
  position=positionfun(dim(tensor))$position
  quan=c(quantile(tensor,(0:color_choice)/color_choice))
  col=tensor
  for(i in 1:color_choice){
    col[(tensor>=quan[i])&(tensor<quan[i+1])]=marker[i]
  }
  col[tensor==quan[i+1]]=marker[i]
  
  plot3d(position[,1],position[,2],position[,3],col=col,alpha=0.3,size=5,xlab="",ylab="",zlab="")
}


## function to find the array index
positionfun=function(d){
  x=rep(1,d[2])%x%rep(1,d[3])%x%c(1:d[1])
  y=rep(1,d[3])%x%(c(1:d[2])%x%rep(1,d[1]))
  z=c(1:d[3])%x%rep(1,d[1])%x%rep(1,d[2])
  position=cbind(x,y,z)
  return(list("position"=position))
}


d1=40 ## number of rows
d2=40 ## number of columns
d3=50 ## number of 3rd dimension
K1=3  ## number of row clusters 
K2=4  ## number of column clusters 
K3=5  ## number of clusters in the 3rd dimension

cluster1_index=sample(1:K1,d1,replace=T) ## randomly assign d1 rows to K1 clusters
cluster2_index=sample(1:K2,d2,replace=T) 
cluster3_index=sample(1:K3,d3,replace=T)

K1=length(unique(cluster1_index)) ## actual number of clusters in realization; no larger than pre-specified K1
K2=length(unique(cluster2_index))
K3=length(unique(cluster3_index))

group_mean=array(runif(K1*K2*K3,-3,3),dim=c(K1,K2,K3)) ## generate group-specific mean
signal=group_mean[cluster1_index,cluster2_index,cluster3_index] ## generate a stochastic block tensor following the groups-spedific mean
noise=array(rnorm(d1*d2*d3,0,0.1),dim=c(d1,d2,d3))

plot_tensor(signal+noise) ## plot data tensor
signal_reorganize=signal[sort(cluster1_index,index=T)$ix,sort(cluster2_index,index=T)$ix,sort(cluster3_index,index=T)$ix] ## re-organize the tensor entries; put entries in the same groups together.
plot_tensor(signal_reorganize) ## plot the signal tensor with reorganized indices

########################################################################################
n = 40
m = 3
gp1 = graphon1(n,m)
gp2 = graphon2(n,m)
gp3 = graphon3(n,m)
gp4 = graphon4(n,m)

A1 = gp1$A; A2 = gp2$A; A3 = gp3$A; A4 = gp4$A
P1 = gp1$P; P2 = gp2$P; P3 = gp3$P; P4 = gp4$P
N3 = array(rnorm(40^3),dim = c(40,40,40))
plot_tensor(P3)
rgl.postscript('P3.pdf', fmt = 'pdf')
plot_tensor(P3+N3)
rgl.postscript('A3.pdf', fmt = 'pdf')

blk = block(P3+N3,c(9,9,9))
plot_tensor(blk$Blk_mat)
rgl.postscript('B3.pdf', fmt = 'pdf')






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








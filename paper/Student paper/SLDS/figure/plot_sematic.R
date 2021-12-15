rate=function(alpha,m){
    n=length(alpha)
    output=rep(0,n)
    for(i in 1:n){
    if(alpha[i]<=m*(m-1)/2)
    output[i]=2*m*alpha[i]/(m+2*alpha[i])
    else
        output[i]=m-1
        }
    return(output)
}

alpha=seq(from=0,to=7,by=0.1)
plot(alpha,rate(alpha,2),type="l",ylim=c(0,3),ylab="mean squared error rate",las=1,axes=F)
points(alpha,rate(alpha,3),type="l")
points(alpha,rate(alpha,4),type="l")


axis(side = 1)
axis(side = 2)

mseq=seq(from=0,to=4,by=0.1)
aseq=(mseq+1)*mseq/2
points(aseq,mseq,type="l",lty=2)


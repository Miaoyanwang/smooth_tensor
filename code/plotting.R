
bd = function(x) return(x*(x-1)/2)
plot(seq(1,5,length = 100),bd(seq(1,5,length = 100)),type = "l",xlab = "tensor order (m)",ylab =expression(smoothness~(alpha)),axes = F)
axis(1)
axis(2)


polygon(x = c(seq(1,5,length = 100), 1,1), 
        y = c(bd(seq(1,5,length = 100)), 10, 0),
        col = "antiquewhite1")
polygon(x =  c(seq(1,5,length = 100), 5,1),
        y = c(bd(seq(1,5,length = 100)), 0, 0),
        col = "darkslategray3")

points(t(c(1,0)), pch=16)
text(t(c(1,1)), "(1,0)", adj=-.1,cex =.9)

points(t(c(2,1)), pch=16)
text(t(c(2,2)), "(2,1)",cex = .9)

points(t(c(3,3)), pch=16)
text(t(c(3,4)), "(3,3)",cex = .9)

points(t(c(4,6)), pch=16)
text(t(c(3.9,7)), "(4,6)",cex = .9)

points(t(c(5,10)), pch=16)
text(t(c(4.4,9.4)), "(5,10)",cex = .9)

text(t(c(1.5,8)),"A",cex = 1.5)
text(t(c(4.5,2)),"B",cex = 1.5)

library(latex2exp)




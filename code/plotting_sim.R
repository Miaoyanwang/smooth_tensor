source(file = "../functions_blk.R")
load(file = "summary.RData")
library(ggplot2)
# Plotting block size versus MSE according to degree
# Continuous
g11 = ggplot(summary1[summary1$model =="model1"&summary1$type =="conti",], aes(x=groupsize, y=MSE, colour=deg)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.5)+
  labs(x="The block size (k)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));g11
g12 = ggplot(summary1[summary1$model =="model2"&summary1$type =="conti",], aes(x=groupsize, y=MSE, colour=deg)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.5)+
  labs(x="The block size (k)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));g12
g13 = ggplot(summary1[summary1$model =="model3"&summary1$type =="conti",], aes(x=groupsize, y=MSE, colour=deg)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.5)+
  labs(x="The block size (k)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));g13
g14 = ggplot(summary1[summary1$model =="model4"&summary1$type =="conti",], aes(x=groupsize, y=MSE, colour=deg)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.5)+
  labs(x="The block size (k)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));g14
g15 = ggplot(summary1[summary1$model =="model5"&summary1$type =="conti",], aes(x=groupsize, y=MSE, colour=deg)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.5)+
  labs(x="The block size (k)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));g15

library(ggpubr)

# Binary
g21 = ggplot(summary1[summary1$model =="model1"&summary1$type =="binary",], aes(x=groupsize, y=MSE, colour=deg)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.5)+
  labs(x="The block size (k)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));g21
g22 = ggplot(summary1[summary1$model =="model2"&summary1$type =="binary",], aes(x=groupsize, y=MSE, colour=deg)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.5)+
  labs(x="The block size (k)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));g22
g23 = ggplot(summary1[summary1$model =="model3"&summary1$type =="binary",], aes(x=groupsize, y=MSE, colour=deg)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.5)+
  labs(x="The block size (k)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));g23
g24 = ggplot(summary1[summary1$model =="model4"&summary1$type =="binary",], aes(x=groupsize, y=MSE, colour=deg)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.5)+
  labs(x="The block size (k)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));g24
g25 = ggplot(summary1[summary1$model =="model5"&summary1$type =="binary",], aes(x=groupsize, y=MSE, colour=deg)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.5)+
  labs(x="The block size (k)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));g25




## Dimension versus MSE according to different degrees
## Conti
s11 = ggplot(best[best$model =="model1"&best$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s11

s12 = ggplot(best[best$model =="model2"&best$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s12

s13 = ggplot(best[best$model =="model3"&best$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s13

s14 = ggplot(best[best$model =="model4"&best$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s14

s15 = ggplot(best[best$model =="model5"&best$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s15



## Binary
s21 = ggplot(best[best$model =="model1"&best$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s21

s22 = ggplot(best[best$model =="model2"&best$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s22

s23 = ggplot(best[best$model =="model3"&best$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s23

s24 = ggplot(best[best$model =="model4"&best$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s24

s25 = ggplot(best[best$model =="model5"&best$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s25



## Dimension versus MSE according to 3 methods: Borda count, LSE, spectral
## Conti
p11 = ggplot(best2[best2$model =="model1"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count", "LSE","Spectral"));p11

p12 = ggplot(best2[best2$model =="model2"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Borda count", "LSE","Spectral"));p12

p13 = ggplot(best2[best2$model =="model3"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count", "LSE","Spectral"));p13

p14 = ggplot(best2[best2$model =="model4"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Borda count", "LSE","Spectral"));p14

p15 = ggplot(best2[best2$model =="model5"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Borda count", "LSE","Spectral"));p15



## Dimension versus MSE according to 4 methods: Borda count, LSE, BAL, spectral
## Binary
p21 = ggplot(best2[best2$model =="model1"&best2$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("BAL","Borda count", "LSE","Spectral"));p21

p22 = ggplot(best2[best2$model =="model2"&best2$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("BAL","Borda count", "LSE","Spectral"));p22

p23 = ggplot(best2[best2$model =="model3"&best2$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("BAL","Borda count", "LSE","Spectral"));p23

p24 = ggplot(best2[best2$model =="model4"&best2$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("BAL","Borda count", "LSE","Spectral"));p24

p25 = ggplot(best2[best2$model =="model5"&best2$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("BAL","Borda count", "LSE","Spectral"));p25



### Estimating missing entries

save(summary31,summary32,file ="missing.RData")
load(file = "missing.RData")

## Conti
q11 = ggplot(summary31[summary31$model =="model1",], aes(x=fraction, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Observation fraction", y="MSE") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count","Spectral"));q11


q12 = ggplot(summary31[summary31$model =="model2",], aes(x=fraction, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Observation fraction", y="MSE") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count","Spectral"));q12


q13 = ggplot(summary31[summary31$model =="model3",], aes(x=fraction, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Observation fraction", y="MSE") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count","Spectral"));q13

q14 = ggplot(summary31[summary31$model =="model4",], aes(x=fraction, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Observation fraction", y="MSE") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count","Spectral"));q14

q15 = ggplot(summary31[summary31$model =="model5",], aes(x=fraction, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Observation fraction", y="MSE") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count","Spectral"));q15



## Conti
q21 = ggplot(summary32[summary32$model =="model1",], aes(x=fraction, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Observation fraction", y="MSE") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count","Spectral"));q21


q22 = ggplot(summary32[summary32$model =="model2",], aes(x=fraction, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Observation fraction", y="MSE") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count","Spectral"));q22


q23 = ggplot(summary32[summary32$model =="model3",], aes(x=fraction, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Observation fraction", y="MSE") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count","Spectral"));q23

q24 = ggplot(summary32[summary32$model =="model4",], aes(x=fraction, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Observation fraction", y="MSE") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count","Spectral"));q24

q25 = ggplot(summary32[summary32$model =="model5",], aes(x=fraction, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Observation fraction", y="MSE") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count","Spectral"));q25






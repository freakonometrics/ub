
################################## Arthur Charpentier, @freakonometrics, July 2018

library(devtools)
 devtools::install_github("pbiecek/breakDown")
 devtools::install_github("redichh/shapleyr")
 list.of.packages <- c("rTensor","HistData","DALEX","parallel","quantreg","auditor","mlr","lime","breakDown","lpSolve","expectreg","VGAMdata","McSpatial","MASS","VIM","locfit","bsplines","rgeos","rgdal","maptools","cartography","geosphere","RDDtools","RColorBrewer","glmnet","caret","splines","factoextra","ROCR","nnet","neuralnet","quadprog","rpart","rpart.plot","gbm","networkD3")
 new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
 if(length(new.packages)) install.packages(new.packages)

library(breakDown,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(plyr,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(DALEX,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(dplyr,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(tm,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(VGAM,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(data.table,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(caret,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(glmnet,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(caret,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(recipes,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)
library(rgl,verbose=FALSE,quietly=TRUE,warn.conflicts=FALSE)

 
 #### --------- #1 INTRODUCTION
 
download.file("https://www.macalester.edu/~abeverid/data/stormofswords.csv","got.csv")
GoT=read.csv("got.csv")
library(networkD3, quietly=TRUE)
simpleNetwork(GoT[,1:2])

M=(rbind(as.matrix(GoT[,1:2]),as.matrix(GoT[,2:1])))
nodes=unique(M[,1])
friends=function(x) as.character(M[which(M[,1]==x),2])
friends_of_friends=function(y) (Vectorize(function(x) length(friends(x)))(friends(y)))
nb_friends = Vectorize(function(x) length(friends(x)))
nb_friends_of_friends = Vectorize(function(x) mean(friends_of_friends(x)))
Nb=nb_friends(nodes)
Nb2=nb_friends_of_friends(nodes)
hist(Nb,breaks=0:40,col=rgb(1,0,0,.2),border="white",probability = TRUE)
hist(Nb2,breaks=0:40,col=rgb(0,0,1,.2),border="white",probability = TRUE,add=TRUE)
lines(density(Nb),col="red",lwd=2)
lines(density(Nb2),col="blue",lwd=2)

chicago = read.table("http://freakonometrics.free.fr/chicago.txt",header=TRUE,sep=";")
tail(chicago)

myocarde = read.table("http://freakonometrics.free.fr/myocarde.csv",head=TRUE, sep=";")
tail(myocarde)
myocarde$PRONO = (myocarde$PRONO=="SURVIE")*1

x1 = c(.4,.55,.65,.9,.1,.35,.5,.15,.2,.85)
x2 = c(.85,.95,.8,.87,.5,.55,.5,.2,.1,.3)
y  = c(1,1,1,1,1,0,0,1,0,0)
df = data.frame(x1=x1,x2=x2,y=as.factor(y))
tail(df)

Davis=read.table(
  "http://socserv.socsci.mcmaster.ca/jfox/Books/Applied-Regression-2E/datasets/Davis.txt")
Davis[12,c(2,3)]=Davis[12,c(3,2)]
tail(Davis)
Y=Davis$height

library(HistData)
attach(Galton)
Galton$count <- 1
df <- aggregate(Galton, by=list(parent,child), FUN=sum)[,c(1,2,5)]
plot(df[,1:2],cex=sqrt(df[,3]/3),xlab="Height of father",ylab="Height to son")
abline(a=0,b=1,lty=2)
abline(lm(child~parent,data=Galton),col="red",lwd=2)
coefficients(lm(child~parent,data=Galton))[2]

plot(chicago$X_2,chicago$Fire)
reg <- lm(Fire~X_2, data=chicago)
vx <- seq(min(chicago$X_2)-5,max(chicago$X_2)+5,length=151)
vy <- predict(reg, newdata=data.frame(X_2=vx))
lines(vx, vy, lwd=2, col="red")

summary(reg)
confint(reg)
beta1 = summary(reg)$coefficients[2,]
vbx <- seq(qnorm(.001,beta1[1],beta1[2]),
           qnorm(.999,beta1[1],beta1[2]),length=201)
vby <- dnorm(vbx, beta1[1],beta1[2])
plot(vbx,vby,lwd=2,type="l")
abline(v=0,lty=2,col="red")

vy2 <- predict(reg, newdata=data.frame(X_2=vx), interval = "confidence")
plot(chicago$X_2,chicago$Fire)
polygon(c(vx,rev(vx)),c(vy2[,2],rev(vy2[,3])),col=rgb(1,0,0,.4),border=NA)
lines(vx, vy2[,1], lwd=2, col="red")
lines(vx, vy2[,2], lty=2, col="red")
lines(vx, vy2[,3], lty=2, col="red")

reg=lm(Fire~X_2+X_3,data=chicago)
summary(reg)
y=function(x2,x3) predict(reg,newdata=data.frame(X_2=x2,X_3=x3))
VX2=seq(0,80)
VX3=seq(5,25)
VY=outer(VX2,VX3,y)
image(VX2,VX3,VY,xlab="",ylab="")
contour(VX2,VX3,VY,add=TRUE)

persp(VX2,VX3,VY,theta=30,ticktype="detailed")

library(ellipse)
plot(ellipse(reg,which=c(2,3)),type="l",col="red",lwd=2)
abline(v=confint(reg)[2,],col="blue",lwd=2)
abline(h=confint(reg)[3,],col="blue",lwd=2)

reg=lm(Fire~.,data=chicago)
summary(reg)
step(reg)

library("DALEX")
apartments_lm_model <- lm(m2.price ~ construction.year + surface + floor + no.rooms + district, data = apartments)
summary(apartments_lm_model)
library(auditor)
explainer_lm <- DALEX::explain(apartments_lm_model, data = apartmentsTest[,2:6], y = apartmentsTest$m2.price)

predicted_mi2_lm <- predict(apartments_lm_model, apartmentsTest)
sqrt(mean((predicted_mi2_lm - apartmentsTest$m2.price)^2))

sv_lm  <- variable_response(explainer_lm, variable =  "construction.year", type = "pdp")
plot(sv_lm)
sv_lm  <- variable_response(explainer_lm, variable =  "surface", type = "pdp")
plot(sv_lm)

svd_lm  <- variable_response(explainer_lm, variable = "district", type = "factor")
plot(svd_lm)

audit_lm <- audit(explainer_lm)
plotResidual(audit_lm, variable = "construction.year")
vi_lm <- variable_importance(explainer_lm, loss_function = loss_root_mean_square)
vi_lm
plot(vi_lm)

library(mlr)
apartments_task <- makeRegrTask(data = apartments, target = "m2.price")
apartments_lm_mlr <- mlr::train("regr.lm", apartments_task)
library(lime)
explained_prediction <- apartments[836, ]
lime_explainer <- lime(apartments, model = apartments_lm_mlr)
lime_explanation <- lime::explain(apartments[836, ], explainer = lime_explainer,  n_features = 5)
plot_features(lime_explanation)
 
library(live)
set.seed(33)
new_dataset <- sample_locally2(data = apartments,
                               explained_instance = apartments[836, ],
                               explained_var = "m2.price",
                               size = 1500)
with_predictions <- add_predictions2(new_dataset, apartments_lm_model)
live_explanation <- fit_explanation2(with_predictions, "regr.lm")
live_explanation
 
plot_explanation2(live_explanation, "forest")
plot_explanation2(live_explanation, "waterfall")

breakdown_linear <- single_prediction(explainer_lm, apartments[836, ])
plot(breakdown_linear)
br <- broken(apartments_lm_model, apartments[836, ], baseline = "Intercept")
br
plot(br) 

library(breakDown)
predict.function <- function(model, new_observation) stats::predict.lm(model, newdata=new_observation)
explain_1 <- broken(apartments_lm_model, apartments[836, ], data = apartments, predict.function = predict.lm, direction = "down")
plot(explain_1)

mp_lm <- model_performance(explainer_lm)
mp_lm
plot(mp_lm, geom = "boxplot")

ggplot(mp_lm, aes(observed, diff)) + geom_point() + geom_smooth(se = FALSE) +
  xlab("Observed") + ylab("Predicted - Observed") +
  ggtitle("Diagnostic plot for the linear model") + theme_mi2()

new_apartment <- apartmentsTest[1234, ]
new_apartment
new_apartment_lm <- single_prediction(explainer_lm,
                                      observation = new_apartment)

   #### --------- #2 MONTE CARLO

runif(30)
runif(30)
set.seed(1)
runif(30)
set.seed(1)
runif(30)

h=function(u) cos(u*pi/2)
integrate(h,0,1)
mean(h(runif(1e6)))
M=rep(NA,1000)
for(i in 1:1000) M[i]=mean(h(runif(1e6)))
mean(M)
sd(M)

F = ecdf(Y)
F(180)
library(MASS)
fitdistr(Y,"normal")
mean(Y)
sd(Y)
y0 = pnorm(140:205,170,9)
y = Vectorize(F)(140:205)
plot(140:205,y,col="red",type="s")
D = rep(NA,200)
for(s in 1:200){
  X = rnorm(length(Y),170,9)
  y = Vectorize(ecdf(X))(140:205)
  lines(140:205,y,col=rgb(0,0,1,.4))
  D[s] = max(y-y0)
}
lines(140:205,y,col="red",type="s")

hist(D,probability = TRUE	)
lines(density(D),col="blue")
(demp = max(abs(Vectorize(ecdf(Y))(140:205)-y0)))
mean(D>demp)
ks.test(Y, "pnorm",170,9)

n = 20
ns = 1e6
xbar = rep(NA,ns)
for(i in 1:ns){
  x = (rchisq(n,df=1)-1)/sqrt(2)
  xbar[i] = mean(x)
}
u = seq(-.7,.8,by=.001)
v = dnorm(u,sd=1/sqrt(20))
plot(u,v,col="black")
lines(density(xbar),col="red")
set.seed(1)
x = (rchisq(n,df=1)-1)/sqrt(2)
for(i in 1:ns){
  xs = sample(x,size=n,replace=TRUE)
  xbar[i] = mean(xs)
}
lines(density(xbar),col="blue")

VS=matrix(NA,15,3)
for(s in 1:15){
  simu=function(n = 10){
    get_i = function(i){
      x = rnorm(n, sd=sqrt(6));
      S = matrix(sample(x, size=n*1000, replace=TRUE),ncol=1000)
      ThetaBoot = exp(colMeans(S))
      Bias = mean(ThetaBoot)-exp(mean(x))
      theta=exp(mean(x))/exp(.5*var(x)/n)
      c(exp(mean(x)),exp(mean(x))-Bias, theta)
    }
    res = lapply(1:1000, get_i)
    res = do.call(rbind, res)
    bias = colMeans(res-1)
    return(bias)
  }
  VS[s,]=simu(10*s)
}
VS

plot(cars)
reg=lm(dist~speed,data=cars)
abline(reg,col="red")
x=21
predict(reg,interval="confidence", level=.9,newdata=data.frame(speed=x))
Yx=rep(NA,500)
for(s in 1:500){
  indice=sample(1:n,size=n,replace=TRUE)
  base=cars[indice,]
  regb=lm(dist~speed,data=base)
  abline(regb,col="light blue")
  points(x,predict(regb,newdata=data.frame(speed=x)))
  Yx[s]=predict(reg,newdata=data.frame(speed=x))
}

predict(reg,interval="confidence", level=.9,newdata=data.frame(speed=x))
hist(Yx,proba=TRUE)
boxplot(Yx,horizontal=TRUE)
lines(density(Yx))
quantile(Yx,c(.05,.95))

plot(cars)
reg=lm(dist~speed,data=cars)
abline(reg,col="red")
x=21
predict(reg,interval="confidence", level=.9,newdata=data.frame(speed=x))
base=cars
Yx=rep(NA,500)
for(s in 1:500){
  indice=sample(1:n,size=n,replace=TRUE)
  base$dist=predict(reg)+residuals(reg)[indice]
  regb=lm(dist~speed,data=base)
  abline(reg,col="light blue")
  points(x,predict(reg,newdata=data.frame(speed=x)))
  Yx[s]=predict(reg,newdata=data.frame(speed=x))
}

predict(reg,interval="confidence", level=.9,newdata=data.frame(speed=x))
hist(Yx,proba=TRUE)
boxplot(Yx,horizontal=TRUE)
lines(density(Yx))

pvf = function(t) mean((1-pf(t,1,length(t)-2))<.05)
pvq = function(t) mean((1-pchisq(t,1)<.05))
TABLE= function(n=30){
  ns = 5000
  x = c(1.0001,rep(1,n-1))
  e = matrix(rnorm(n*ns),n)
  e2 = matrix(runif(n*ns,-3,3),n)
  e3 = matrix(rt(n*ns,2),n)
  get_i = function(i){
    r1 = lm(e[,i]~x)
    r2 = lm(e2[,i]~x)
    r3 = lm(e3[,i]~x)
    t1 = r1$coef[2]^2/vcov(r1)[2,2]
    t2 = r2$coef[2]^2/vcov(r2)[2,2]
    t3 = r3$coef[2]^2/vcov(r3)[2,2]
    c(t1,t2,t3)}
  
  # library(parallel)	
  # t = mclapply(1:ns, get_i, mc.cores=50)
  t = lapply(1:ns, get_i)
  t = simplify2array(t)
  rj1 = pvf(t[,1])
  rj2 = pvf(t[,2])
  rj3 = pvf(t[,3])
  rj12 = pvq(t[,1])
  rj22 = pvq(t[,2])
  rj32 = pvq(t[,3])
  ans = rbind(c(rj1,rj2,rj3),c(rj12,rj22,rj32))
  return(ans) }
TABLE(30)

ns=1e5
PROP=matrix(NA,ns,6)
n=30
VN=seq(10,140,by=10)
for(s in 1:ns){
  X=rnorm(n)
  E=rnorm(n)
  Y=1+X+E
  reg=lm(Y~X)
  T=(coefficients(reg)[2]-1)^2/ vcov(reg)[2,2]
  PROP[s,1]=T>qf(.95,1,n-2)
  PROP[s,2]=T>qchisq(.95,1)
  E=rt(n,df=3)
  Y=1+X+E
  reg=lm(Y~X)
  T=(coefficients(reg)[2]-1)^2/ vcov(reg)[2,2]
  PROP[s,3]=T>qf(.95,1,n-2)
  PROP[s,4]=T>qchisq(.95,1)
  E=runif(n)*4-2
  Y=1+X+E
  reg=lm(Y~X)
  T=(coefficients(reg)[2]-1)^2/ vcov(reg)[2,2]
  PROP[s,5]=T>qf(.95,1,n-2)
  PROP[s,6]=T>qchisq(.95,1)
}
apply(PROP,2,mean)

reg=lm(dist~speed,data=cars)
E=residuals(reg)
Y=predict(reg)
beta1=rep(NA,2000)
MB=MV=MQ=matrix(NA,10,2000)
for(i in 1:10){
  BIAS=VAR=QUANT=beta1
  for(b in 1:2000){
    carsS=data.frame(speed=cars$speed,
                     dist=Y+sample(E,size=50,replace=TRUE))
    beta1[b]=lm(dist~speed,data=carsS)$coefficients[2]
    BIAS[b]=mean(beta1[1:b])-reg$coefficients[2]
    VAR[b]=var(beta1[1:b])
    QUANT[b]=quantile(beta1[1:b],.95)
  }
  MB[i,]=BIAS
  MV[i,]=VAR
  MQ[i,]=QUANT
}

x1 = c(.4,.55,.65,.9,.1,.35,.5,.15,.2,.85)
x2 = c(.85,.95,.8,.87,.5,.55,.5,.2,.1,.3)
y  = c(1,1,1,1,1,0,0,1,0,0)
df = data.frame(x1=x1,x2=x2,y=as.factor(y))
L_logit = list()
n = nrow(df)
for(s in 1:100){
  df_s = df[sample(1:n,size=n,replace=TRUE),]
  L_logit[[s]] = glm(y~., df_s, family=binomial)}
p = function(x){
  nd = data.frame(x1=x[1], x2=x[2]) 
  unlist(lapply(1:1000,function(z)
    predict(L_logit[[z]],newdata=nd,type="response")))}
```

```{r}
L_logit = list()
n = nrow(df)
reg = glm(y~x1+x2, df, family=binomial)
for(s in 1:100){
  df_s = df
  df_s$y = factor(rbinom(n,size=1,prob=predict(reg,type="response")),labels=0:1)
  L_logit[[s]] = glm(y~., df_s, family=binomial)
}

#### --------- #3 LOSS FUNCTION

# ols <- lm(y~x, data=df)
# library(quantreg)
# lad <- rq(y~x, data=df, tau=.5)

n = 101 
set.seed(1)
y = rlnorm(n)
median(y)
md=Vectorize(function(m) sum(abs(y-m)))
optim(mean(y),md)$par
library(lpSolve)
# in case of problem, see https://stackoverflow.com/questions/36974206/error-maximal-number-of-dlls-reached
A1 = cbind(diag(2*n),0) 
A2 = cbind(diag(n), -diag(n), 1)
r = lp("min", c(rep(1,2*n),0),
rbind(A1, A2),c(rep("&gt;=", 2*n), rep("=", n)), c(rep(0,2*n), y))
tail(r$solution,1) 

tau = .3
quantile(x,tau)
A1 = cbind(diag(2*n),0) 
A2 = cbind(diag(n), -diag(n), 1)
r = lp("min", c(rep(tau,n),rep(1-tau,n),0),
rbind(A1, A2),c(rep("&gt;=", 2*n), rep("=", n)), c(rep(0,2*n), y))
tail(r$solution,1) 

x <- rnorm(99)
library(expectreg)
# e <- expectile(x, probs = seq(0, 1, 0.1))

library(quantreg)
fit =  rq(dist~speed, data=cars, tau=.5)
which(predict(fit)==cars$dist)

base=read.table("http://freakonometrics.free.fr/rent98_00.txt",header=TRUE)
require(lpSolve) 
tau = .3
n=nrow(base)
X = cbind( 1, base$area)
y = base$rent_euro
A1 = cbind(diag(2*n), 0,0) 
A2 = cbind(diag(n), -diag(n), X) 
r = lp("min",
       c(rep(tau,n), rep(1-tau,n),0,0), rbind(A1, A2),
       c(rep("&gt;=", 2*n), rep("=", n)), c(rep(0,2*n), y)) 
 tail(r$solution,2)
 library(quantreg)
rq(rent_euro~area, tau=tau, data=base)

plot(base$area,base$rent_euro,xlab=expression(paste("surface (",m^2,")")),
     ylab="rent (euros/month)",col=rgb(0,0,1,.4),cex=.5)
sf=0:250
yr=r$solution[2*n+1]+r$solution[2*n+2]*sf
lines(sf,yr,lwd=2,col="blue")
tau = .9
r = lp("min",c(rep(tau,n), rep(1-tau,n),0,0), rbind(A1, A2),c(rep(">=", 2*n), rep("=", n)), c(rep(0,2*n), y)) 
 yr=r$solution[2*n+1]+r$solution[2*n+2]*sf
 lines(sf,yr,lwd=2,col="blue")

tau = 0.3
n = nrow(base)
X = cbind(1, base$area, base$yearc)
y = base$rent_euro
r = lp("min",
c(rep(tau, n), rep(1 - tau, n), rep(0, 2 * 3)),
cbind(diag(n), -diag(n), X, -X),
rep("=", n),y)
beta = tail(r$solution, 6)
beta = beta[1:3] - beta[3 + 1:3]
beta
library(quantreg)
rq(rent_euro~area+yearc, tau=tau, data=base)

base = read.table("http://freakonometrics.free.fr/natality2005.txt",header=TRUE,sep=";")
u = seq(.05,.95,by=.01)
library(quantreg)
coefstd = function(u) summary(rq(WEIGHT~SEX+SMOKER+WEIGHTGAIN+BIRTHRECORD+AGE+ BLACKM+ BLACKF+COLLEGE,data=sbase,tau=u))$coefficients[,2]
coefest = function(u) summary(rq(WEIGHT~SEX+SMOKER+WEIGHTGAIN+BIRTHRECORD+AGE+ BLACKM+ BLACKF+COLLEGE,data=sbase,tau=u))$coefficients[,1]
CS = Vectorize(coefstd)(u)
CE = Vectorize(coefest)(u)

 base = read.table("http://freakonometrics.free.fr/BWeight.csv")
 base = read.table("http://freakonometrics.free.fr/rent98_00.txt")

 library(VGAMdata)
 data(xs.nz)

 library(McSpatial)
data(cookdata)
fit <- qregcpar(LNFAR~DCBD, nonpar=~LATITUDE+LONGITUDE, taumat=c(.10,.90), kern="bisq", window=.30, distance="LATLONG", data=cookdata)

# library(expectreg)
# coefstd=function(u) summary(expectreg.ls(WEIGHT~SEX+SMOKER+WEIGHTGAIN+BIRTHRECORD+AGE+ BLACKM+ BLACKF+COLLEGE,data=sbase,expectiles=u,ci = TRUE))[,2]
# coefest=function(u) summary(expectreg.ls(WEIGHT~SEX+SMOKER+WEIGHTGAIN+BIRTHRECORD+AGE+ BLACKM+ BLACKF+COLLEGE,data=sbase,expectiles=u,ci = TRUE))[,1]
# CS = Vectorize(coefstd)(u)
# CE = Vectorize(coefest)(u)

# fit <- expectreg.ls(rent_euro ~ area, data=munich, expectiles=.75)
# fit <- expectreg.ls(rent_euro ~ rb(area,"pspline"), data=munich, expectiles=.75)

#### --------- #4 NON LINEARITIES

height = Davis$height
hist(height)
plot(density(height, kernel = "rectangular"))
plot(density(height, kernel = "triangular"))
plot(density(height, kernel = "epanechnikov"))
plot(density(height, kernel = "gaussian"))

bw.nrd0(cars$speed)
bw.nrd(cars$speed)

library(VIM)
y= tao[,c("Air.Temp","Humidity")]
histMiss(y) 
y= tao[,c("Humidity","Air.Temp")]
histMiss(y) 

tao_kNN <- kNN(tao, k = 5)
vars <- c("Air.Temp","Humidity","Air.Temp_imp","Humidity_imp")
marginplot(tao_kNN[,vars], delimiter="imp", alpha=0.6)

loess(dist ~ speed, cars,span=0.75,degree=1)
predict(REG, data.frame(speed = seq(5, 25, 0.25)), se = TRUE)

library(locfit)
locfit(dist~speed,data=cars)

reg <- lm(y~x,data=db)
reg <- lm(y~poly(x,5),data=db)
reg <- lm(y~poly(x,25),data=db)

positive_part <- function(x) ifelse(x>0,x,0)
reg <- lm(Y~X+positive_part(X-s), data=db)

reg <- lm(Y~X+positive_part(X-s1)+positive_part(X-s2), data=db)
library(bsplines)

reg1 <- lm(dist~speed+positive_part(speed-15), data=cars)
reg2 <- lm(dist~bs(speed,df=2,degree=1), data=cars)
summary(reg1)
summary(reg2)

library(rgeos)
library(rgdal)
library(maptools)
library(cartography)
download.file("http://bit.ly/2G3KIUG","zonier.RData")
load("zonier.RData")
cols=rev(carto.pal(pal1="red.pal",n1=10,pal2="green.pal",n2=10))
download.file("http://bit.ly/2GSvzGW","FRA_adm0.rds")
download.file("http://bit.ly/2FUZ0Lz","FRA_adm2.rds")
FR=readRDS("FRA_adm2.rds")
donnees_carte=data.frame(FR@data)

FR0=readRDS("FRA_adm0.rds")
plot(FR0)
bk = seq(-5,4.5,length=21)
cuty = cut(simbase$Y,breaks=bk,labels=1:20)
points(simbase$long,simbase$lat,col=cols[cuty],pch=19,cex=.5)

A=aggregate(x = simbase$Y,by=list(simbase$dpt),mean)
names(A)=c("dpt","y")
d=donnees_carte$CCA_2
d[d=="2A"]="201"
d[d=="2B"]="202"
donnees_carte$dpt=as.numeric(as.character(d))
donnees_carte=merge(donnees_carte,A,all.x=TRUE)
donnees_carte=donnees_carte[order(donnees_carte$OBJECTID),]
bk=seq(-2.75,2.75,length=21)
donnees_carte$cuty=cut(donnees_carte$y, breaks=bk,labels=1:20)
plot(FR, col=cols[donnees_carte$cuty],xlim=c(-5.2,12))

bk=seq(-2.75,2.75,length=5)
donnees_carte$cuty=cut(donnees_carte$y,breaks=bk,labels=1:4)
plot(FR, col=cols[c(3,8,12,17)][donnees_carte$cuty],xlim=c(-5.2,12))

P1 = FR0@polygons[[1]]@Polygons[[355]]@coords
P2 = FR0@polygons[[1]]@Polygons[[27]]@coords
plot(FR0,border=NA)
polygon(P1)
polygon(P2)
grille<-expand.grid(seq(min(simbase$long),max(simbase$long),length=101),seq(min(simbase$lat),max(simbase$lat),length=101))
paslong=(max(simbase$long)-min(simbase$long))/100
paslat=(max(simbase$lat)-min(simbase$lat))/100

f=function(i){ (point.in.polygon (grille[i, 1]+paslong/2 , grille[i, 2]+paslat/2 , P1[,1],P1[,2])>0)+(point.in.polygon (grille[i, 1]+paslong/2 , grille[i, 2]+paslat/2 , P2[,1],P2[,2])>0) }
indic=unlist(lapply(1:nrow(grille),f))
grille=grille[which(indic==1),]
points(grille[,1]+paslong/2,grille[,2]+paslat/2)

 library(geosphere)
 knn=function(i,k=20){
   d=distHaversine(grille[i,1:2],simbase[,c("long","lat")], r=6378.137)
    r=rank(d)
    ind=which(r<=k)
    mean(simbase[ind,"Y"])
   }
 grille$y=Vectorize(knn)(1:nrow(grille))
 bk=seq(-2.75,2.75,length=21)
 grille$cuty=cut(grille$y,breaks=bk,labels=1:20)
 points(grille[,1]+paslong/2,grille[,2]+paslat/2,col=cols[grille$cuty],pch=19)

 bk=seq(-2.75,2.75,length=5)
 grille$cuty=cut(grille$y,breaks=bk,labels=1:4)
 plot(FR0,border=NA)
 polygon(P1)
 polygon(P2)
 points(grille[,1]+paslong/2,grille[,2]+paslat/2,col=cols[c(3,8,12,17)][grille$cuty],pch=19)

 for(t in 1:100){fit = rpart(yr~x,data=df)
 yp = predict(fit,newdata=df)
 df$yr = df$yr - v*yp)}
 
 Fstats(dist ~ speed,data=cars,from=7/50)
 
 Fstats(dist ~ speed,data=cars,from=2/50)
 
 cusum <- efp(dist ~ speed, type = "OLS-CUSUM",data=cars)
 plot(cusum,ylim=c(-2,2))
 plot(cusum, alpha = 0.05, alt.boundary = TRUE,ylim=c(-2,2))
 
 library(RDDtools)
 data(Lee2008) 
 
 idx1 = (Lee2008$x>0)
 reg1 = lm(y~poly(x,4),data=Lee2008[idx1,])
 idx2 = (Lee2008$x<0)
 reg2 = lm(y~poly(x,4),data=Lee2008[idx2,])
 s1=predict(reg1,newdata=data.frame(x=0))
 s2=predict(reg2,newdata=data.frame(x=0))
 abs(s1-s2)
 
 reg_para <- RDDreg_lm(RDDdata(y = Lee2008$y, x = Lee2008$x, cutpoint = 0), order = 4)
 reg_para
 
 reg1 = ksmooth(Lee2008$x[idx1], Lee2008$y[idx1], kernel = "normal", bandwidth = 0.1)
 reg2 = ksmooth(Lee2008$x[idx2], Lee2008$y[idx2], kernel = "normal", bandwidth = 0.1)
 s1 = reg1$y[1]
 s2 = reg2$y[length(reg2$y)]
 abs(s1-s2)
 
 reg_nonpara <- RDDreg_np(RDDobject = Lee2008_rdd, bw = .1)
 print(reg_nonpara)
 
 library(glmnet)
 chicago=read.table("http://freakonometrics.free.fr/chicago.txt",header=TRUE,sep=";")
 standardize <-  function(x)  {(x-mean(x))/sd(x)}
 z0 <- standardize(chicago[, 1])
 z1 <- standardize(chicago[, 3])
 z2 <- standardize(chicago[, 4])
 ridge <-glmnet(cbind(z1, z2), z0, alpha=0, intercept=FALSE, lambda=1)
 lasso <-glmnet(cbind(z1, z2), z0, alpha=1, intercept=FALSE, lambda=1)
 elastic <-glmnet(cbind(z1, z2), z0, alpha=.5, intercept=FALSE, lambda=1)
 
 fit <- lm(y ~ x, data = df)
 
 fit <- lm(y ~ poly(x, k), data = df)
 fit <- lm(y ~ bs(x), data = df)
 
 fit <- glm(y ~ x, data = df, family = poisson(link = "log"))
            
 y = myocarde$PRONO
 X = cbind(1,as.matrix(myocarde[,1:7]))
 negLogLik = function(beta){
   -sum(-y*log(1 + exp(-(X%*%beta))) - (1-y)*log(1 + exp(X%*%beta)))
 }
 beta_init = lm(PRONO~.,data=myocarde)$coefficients
 logistic_opt = optim(par = beta_init, negLogLik, hessian=TRUE, method = "BFGS", control=list(abstol=1e-9))
 
 logistic_opt$par
 simu = function(i){
   logistic_opt_i = optim(par = rnorm(8,0,3)*beta_init, 
                          negLogLik, hessian=TRUE, method = "BFGS", 
                          control=list(abstol=1e-9))
   logistic_opt_i$par[2:3]
 }
 v_beta = t(Vectorize(simu)(1:1000))
 
 plot(v_beta)
 par(mfrow=c(1,2))
 hist(v_beta[,1],xlab=names(myocarde)[1])
 hist(v_beta[,2],xlab=names(myocarde)[2])
 
 Y=myocarde$PRONO
 X=cbind(1,as.matrix(myocarde[,1:7]))
 colnames(X)=c("Inter",names(myocarde[,1:7]))
 beta=as.matrix(lm(Y~0+X)$coefficients,ncol=1)
 for(s in 1:9){
   pi=exp(X%*%beta[,s])/(1+exp(X%*%beta[,s]))
   gradient=t(X)%*%(Y-pi)
   omega=matrix(0,nrow(X),nrow(X));diag(omega)=(pi*(1-pi))
   Hessian=-t(X)%*%omega%*%X
   beta=cbind(beta,beta[,s]-solve(Hessian)%*%gradient)}
 beta[,8:10]
 
 df = myocarde
 beta_init = lm(PRONO~.,data=df)$coefficients
 X = cbind(1,as.matrix(myocarde[,1:7]))
 beta = beta_init
 for(s in 1:1000){
   p = exp(X %*% beta) / (1+exp(X %*% beta))
   omega = diag(nrow(df))
   diag(omega) = (p*(1-p))
   df$Z = X %*% beta + solve(omega) %*% (df$PRONO - p)
   beta = lm(Z~.,data=df[,-8], weights=diag(omega))$coefficients
 }  
 beta
 
 Y=myocarde$PRONO
 X=cbind(1,as.matrix(myocarde[,1:7]))
 colnames(X)=c("Inter",names(myocarde[,1:7]))
 beta=as.matrix(lm(Y~0+X)$coefficients,ncol=1)
 for(s in 1:9){
   pi=exp(X%*%beta[,s])/(1+exp(X%*%beta[,s]))
   gradient=t(X)%*%(Y-pi)
   omega=matrix(0,nrow(X),nrow(X));diag(omega)=(pi*(1-pi))
   Hessian=-t(X)%*%omega%*%X
   beta=cbind(beta,beta[,s]-solve(Hessian)%*%gradient)}
 beta[,8:10]
 
 f = myocarde
 beta_init = lm(PRONO~.,data=df)$coefficients
 X = cbind(1,as.matrix(myocarde[,1:7]))
 beta = beta_init
 for(s in 1:1000){
   p = exp(X %*% beta) / (1+exp(X %*% beta))
   omega = diag(nrow(df))
   diag(omega) = (p*(1-p))
   df$Z = X %*% beta + solve(omega) %*% (df$PRONO - p)
   beta = lm(Z~.,data=df[,-8], weights=diag(omega))$coefficients
 }  
 beta
 
 summary( lm(Z~.,data=df[,-8], weights=diag(omega)))
 summary(glm(PRONO~.,data=myocarde,family=binomial(link = "logit")))
 
 pos = function(x,s) (x-s)*(x>=s)
 
 reg = glm(PRONO~INSYS+pos(INSYS,15)+pos(INSYS,25),data=myocarde,family=binomial)
 summary(reg)
 
 library(splines)
 clr6 = c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02")
 x = seq(0,60,by=.25)
 B = bs(x,knots=c(15,25),Boundary.knots=c(5,55),degre=1)
 matplot(x,B,type="l",lty=1,lwd=2,col=clr6)
 
 reg = glm(PRONO~bs(INSYS,knots=c(15,25), Boundary.knots=c(5,55),degre=1), data=myocarde,family=binomial)
 summary(reg)
 
 pos2 = function(x,s) (x-s)^2*(x>=s)
 reg = glm(PRONO~poly(INSYS,2)+pos2(INSYS,15)+pos2(INSYS,25), data=myocarde,family=binomial)
 summary(reg)
 
 x = seq(0,60,by=.25)
 B=bs(x,knots=c(15,25),Boundary.knots=c(5,55),degre=2)
 matplot(x,B,type="l",xlab="INSYS",col=clr6)
 
 reg = glm(PRONO~bs(INSYS,knots=c(15,25), Boundary.knots=c(5,55),degre=2),data=myocarde, family=binomial)
 
 B=bs(x,knots=c(15,25),Boundary.knots=c(5,55),degre=3)
 matplot(x,B,type="l",lwd=2,col=clr6,lty=1,ylim=c(-.2,1.2))
 abline(v=c(5,15,25,55),lty=2)
 
 reg = glm(PRONO~1+bs(INSYS,degree=1,df=4),data=myocarde,family=binomial)
 attr(reg$terms, "predvars")[[3]]
 bs(INSYS, degree = 1L, knots = c(15.8, 21.4, 27.15), Boundary.knots = c(8.7, 54), intercept = FALSE)
 quantile(myocarde$INSYS,(0:4)/4)
 
 B = bs(x,degree=1,df=4)
 B = cbind(1,B)
 y = B%*%coefficients(reg)
 plot(x,y,type="l",col="red",lwd=2)
 abline(v=quantile(myocarde$INSYS,(0:4)/4),lty=2)
 
 reg = glm(y~bs(x1,degree=1,df=3)+bs(x2,degree=1,df=3),data=df,family=binomial(link = "logit"))
 u = seq(0,1,length=101)
 p = function(x,y) predict.glm(reg,newdata=data.frame(x1=x,x2=y),type="response")
 v = outer(u,u,p)
 image(u,u,v,xlab="Variable 1",ylab="Variable 2",col=clr10,breaks=(0:10)/10)
 points(df$x1,df$x2,pch=19,cex=1.5,col="white")
 points(df$x1,df$x2,pch=c(1,19)[1+(df$y=="1")],cex=1.5)
 contour(u,u,v,levels = .5,add=TRUE)
 
 mean_x = function(x,bw){
   w = dnorm((myocarde$INSYS-x)/bw, mean=0,sd=1)
   weighted.mean(myocarde$PRONO,w)}
 u = seq(5,55,length=201)
 v = Vectorize(function(x) mean_x(x,3))(u)
 plot(u,v,ylim=0:1,type="l",col="red")
 points(myocarde$INSYS,myocarde$PRONO,pch=19)
 
 v = Vectorize(function(x) mean_x(x,2))(u)
 plot(u,v,ylim=0:1,type="l",col="red")
 points(myocarde$INSYS,myocarde$PRONO,pch=19)
 
 reg = ksmooth(myocarde$INSYS,myocarde$PRONO,"normal",bandwidth = 2*exp(1))
 plot(reg$x,reg$y,ylim=0:1,type="l",col="red",lwd=2,xlab="INSYS",ylab="")
 points(myocarde$INSYS,myocarde$PRONO,pch=19)
 
 u = seq(0,1,length=101)
 p = function(x,y){
   bw1 = .2; bw2 = .2
   w = dnorm((df$x1-x)/bw1, mean=0,sd=1)*
     dnorm((df$x2-y)/bw2, mean=0,sd=1)
   weighted.mean(df$y=="1",w)
 }
 v = outer(u,u,Vectorize(p))
 image(u,u,v,col=clr10,breaks=(0:10)/10)
 points(df$x1,df$x2,pch=19,cex=1.5,col="white")
 points(df$x1,df$x2,pch=c(1,19)[1+(df$y=="1")],cex=1.5)
 contour(u,u,v,levels = .5,add=TRUE)
 
 Sigma = var(myocarde[,1:7])
 Sigma_Inv = solve(Sigma)
 d2_mahalanobis = function(x,y,Sinv){as.numeric(x-y)%*%Sinv%*%t(x-y)}
 k_closest = function(i,k){
   vect_dist = function(j) d2_mahalanobis(myocarde[i,1:7],myocarde[j,1:7],Sigma_Inv)
   vect = Vectorize(vect_dist)((1:nrow(myocarde))) 
   which((rank(vect)))}
 
 k_majority = function(k){
   Y=rep(NA,nrow(myocarde))
   for(i in 1:length(Y)) Y[i] = sort(myocarde$PRONO[k_closest(i,k)])[(k+1)/2]
   return(Y)}
 
 k_mean = function(k){
   Y=rep(NA,nrow(myocarde))
   for(i in 1:length(Y)) Y[i] = mean(myocarde$PRONO[k_closest(i,k)])
   return(Y)}
 
 Sigma_Inv = solve(var(df[,c("x1","x2")]))
 u = seq(0,1,length=51)
 p = function(x,y){
   k = 6
   vect_dist = function(j)  d2_mahalanobis(c(x,y),df[j,c("x1","x2")],Sigma_Inv)
   vect = Vectorize(vect_dist)(1:nrow(df)) 
   idx  = which(rank(vect)<=k)
   return(mean((df$y==1)[idx]))}
 
 v = outer(u,u,Vectorize(p))
 image(u,u,v,col=clr10,breaks=(0:10)/10)
 points(df$x1,df$x2,pch=19,cex=1.5,col="white")
 points(df$x1,df$x2,pch=c(1,19)[1+(df$y=="1")],cex=1.5)
 contour(u,u,v,levels = .5,add=TRUE)
 
 PennegLogLik = function(bbeta,lambda=0){
   b0   = bbeta[1]
   beta = bbeta[-1]
   -sum(-y*log(1 + exp(-(b0+X%*%beta))) - (1-y)*
          log(1 + exp(b0+X%*%beta)))+lambda*sum(beta^2)}
 
 pt_ridge = function(lambda){
   beta_init = lm(PRONO~.,data=myocarde)$coefficients
   logistic_opt = optim(par = beta_init*0, function(x) PennegLogLik(x,lambda), method = "BFGS", control=list(abstol=1e-9))
   logistic_opt$par[-1]}
 v_lambda = c(exp(seq(-2,5,length=61)))
 est_ridge = Vectorize(opt_ridge)(v_lambda)
 plot(v_lambda,est_ridge[1,])
 
 Y = myocarde$PRONO
 X = myocarde[,1:7]
 for(j in 1:7) X[,j] = (X[,j]-mean(X[,j]))/sd(X[,j])
 X = as.matrix(X)
 X = cbind(1,X)
 colnames(X) = c("Inter",names(myocarde[,1:7]))
 beta = as.matrix(lm(Y~0+X)$coefficients,ncol=1)
 for(s in 1:9){
   pi = exp(X%*%beta[,s])/(1+exp(X%*%beta[,s]))
   Delta = matrix(0,nrow(X),nrow(X));diag(Delta)=(pi*(1-pi))
   z = X%*%beta[,s] + solve(Delta)%*%(Y-pi)
   B = solve(t(X)%*%Delta%*%X+2*lambda*diag(ncol(X))) %*% (t(X)%*%Delta%*%z)
   beta = cbind(beta,B)}
 beta[,8:10]
 
 newton_ridge = function(lambda=1){
   beta = as.matrix(lm(Y~0+X)$coefficients,ncol=1)*runif(8)
   for(s in 1:20){
     pi = exp(X%*%beta[,s])/(1+exp(X%*%beta[,s]))
     Delta = matrix(0,nrow(X),nrow(X));diag(Delta)=(pi*(1-pi))
     z = X%*%beta[,s] + solve(Delta)%*%(Y-pi)
     B = solve(t(X)%*%Delta%*%X+2*lambda*diag(ncol(X))) %*% (t(X)%*%Delta%*%z)
     beta = cbind(beta,B)}
   Varz = solve(Delta)
   Varb = solve(t(X)%*%Delta%*%X+2*lambda*diag(ncol(X))) %*% t(X)%*% Delta %*% Varz %*%
     Delta %*% X %*% solve(t(X)%*%Delta%*%X+2*lambda*diag(ncol(X)))
   return(list(beta=beta[,ncol(beta)],sd=sqrt(diag(Varb))))}
 
 v_lambda=c(exp(seq(-2,5,length=61)))
 est_ridge=Vectorize(function(x) newton_ridge(x)$beta)(v_lambda)
 library("RColorBrewer")
 colrs=brewer.pal(7,"Set1")
 plot(v_lambda,est_ridge[1,],col=colrs[1],type="l")
 for(i in 2:7) lines(v_lambda,est_ridge[i,],col=colrs[i])
 
 v_lambda=c(exp(seq(-2,5,length=61)))
 est_ridge=Vectorize(function(x) newton_ridge(x)$sd)(v_lambda)
 library("RColorBrewer")
 colrs=brewer.pal(7,"Set1")
 plot(v_lambda,est_ridge[1,],col=colrs[1],type="l")
 for(i in 2:7) lines(v_lambda,est_ridge[i,],col=colrs[i],lwd=2)
 
 y = myocarde$PRONO
 X = myocarde[,1:7]
 for(j in 1:7) X[,j] = (X[,j]-mean(X[,j]))/sd(X[,j])
 X = as.matrix(X)
 library(glmnet)
 glm_ridge = glmnet(X, y, alpha=0)
 plot(glm_ridge,xvar="lambda",col=colrs,lwd=2)
 
 library(factoextra)
 pca = princomp(X)
 pca_X = get_pca_ind(pca)$coord
 
 library(glmnet)
 glm_ridge = glmnet(pca_X, y, alpha=0)
 plot(glm_ridge,xvar="lambda",col=colrs,lwd=2)
 
 plot(glm_ridge,col=colrs,lwd=2)
 
 df0 = df
 df0$y=as.numeric(df$y)-1
 plot_lambda = function(lambda){
   m = apply(df0,2,mean)
   s = apply(df0,2,sd)
   for(j in 1:2) df0[,j] = (df0[,j]-m[j])/s[j]
   reg = glmnet(cbind(df0$x1,df0$x2), df0$y==1, alpha=0,lambda=lambda)
   u = seq(0,1,length=101)
   p = function(x,y){
     xt = (x-m[1])/s[1]
     yt = (y-m[2])/s[2]
     predict(reg,newx=cbind(x1=xt,x2=yt),type='response')}
   v = outer(u,u,p)
   image(u,u,v,col=clr10,breaks=(0:10)/10)
   points(df$x1,df$x2,pch=c(1,19)[1+z],cex=1.5)
   contour(u,u,v,levels = .5,add=TRUE)
 }

   reg = glmnet(cbind(df0$x1,df0$x2), df0$y==1, alpha=0)
   par(mfrow=c(1,2))
   plot(reg,xvar="lambda",col=c("blue","red"),lwd=2)
   abline(v=log(.2))
   plot_lambda(.2)
   
   reg = glmnet(cbind(df0$x1,df0$x2), df0$y==1, alpha=0)
   par(mfrow=c(1,2))
   plot(reg,xvar="lambda",col=c("blue","red"),lwd=2)
   abline(v=log(1.2))
   plot_lambda(1.2) 
   
   y = myocarde$PRONO
   X = myocarde[,1:7]
   for(j in 1:7) X[,j] = (X[,j]-mean(X[,j]))/sd(X[,j])
   X = as.matrix(X)
   LogLik = function(bbeta){
     b0=bbeta[1]
     beta=bbeta[-1]
     sum(-y*log(1 + exp(-(b0+X%*%beta))) - 
           (1-y)*log(1 + exp(b0+X%*%beta)))}
   u = seq(-4,4,length=251)
   v = outer(u,u,function(x,y) LogLik(c(1,x,y)))
   image(u,u,v,col=rev(heat.colors(25)))
   contour(u,u,v,add=TRUE)
   polygon(c(-1,0,1,0),c(0,1,0,-1),border="blue")
   
   PennegLogLik = function(bbeta,lambda=0){
     b0=bbeta[1]
     beta=bbeta[-1]
     -sum(-y*log(1 + exp(-(b0+X%*%beta))) - 
            (1-y)*log(1 + exp(b0+X%*%beta)))+lambda*sum(abs(beta))
   }
   opt_lasso = function(lambda){
     beta_init = lm(PRONO~.,data=myocarde)$coefficients
     logistic_opt = optim(par = beta_init*0, function(x) PennegLogLik(x,lambda), 
                          hessian=TRUE, method = "BFGS", control=list(abstol=1e-9))
     logistic_opt$par[-1]
   }
   v_lambda=c(exp(seq(-4,2,length=61)))
   est_lasso=Vectorize(opt_lasso)(v_lambda)
   library("RColorBrewer")
   colrs=brewer.pal(7,"Set1")
   plot(v_lambda,est_lasso[1,],col=colrs[1],type="l")
   for(i in 2:7) lines(v_lambda,est_lasso[i,],col=colrs[i],lwd=2)
   
   library(glmnet)
   glm_lasso = glmnet(X, y, alpha=1)
   plot(glm_lasso,xvar="lambda",col=colrs,lwd=2)
   
   plot(glm_lasso,col=colrs,lwd=2)
   glmnet(X, y, alpha=1,lambda=exp(-4))$beta
   
   library(factoextra)
   pca = princomp(X)
   pca_X = get_pca_ind(pca)$coord
   glm_lasso = glmnet(pca_X, y, alpha=1)
   plot(glm_lasso,xvar="lambda",col=colrs)
   plot(glm_lasso,col=colrs)
   
   soft_thresholding = function(x,a){
     result = numeric(length(x))
     result[which(x > a)] <- x[which(x > a)] - a
     result[which(x < -a)] <- x[which(x < -a)] + a
     return(result)
   }
   
   lasso_coord_desc = function(X,y,beta,lambda,tol=1e-6,maxiter=1000){
     beta = as.matrix(beta)
     X = as.matrix(X)
     omega = rep(1/length(y),length(y))
     obj = numeric(length=(maxiter+1))
     betalist = list(length(maxiter+1))
     betalist[[1]] = beta
     beta0list = numeric(length(maxiter+1))
     beta0 = sum(y-X%*%beta)/(length(y))
     beta0list[1] = beta0
     for (j in 1:maxiter){
       for (k in 1:length(beta)){
         r = y - X[,-k]%*%beta[-k] - beta0*rep(1,length(y))
         beta[k] = (1/sum(omega*X[,k]^2))*soft_thresholding(t(omega*r)%*%X[,k],length(y)*lambda)
       }
       beta0 = sum(y-X%*%beta)/(length(y))
       beta0list[j+1] = beta0
       betalist[[j+1]] = beta
       obj[j] = (1/2)*(1/length(y))*norm(omega*(y - X%*%beta - 
                                                  beta0*rep(1,length(y))),'F')^2 + lambda*sum(abs(beta))
       if (norm(rbind(beta0list[j],betalist[[j]]) - rbind(beta0,beta),'F') < tol) { break } 
     } 
     return(list(obj=obj[1:j],beta=beta,intercept=beta0)) }
 
   lasso_coord_desc = function(X,y,beta,lambda,tol=1e-6,maxiter=1000){
     beta = as.matrix(beta)
     X = as.matrix(X)
     obj = numeric(length=(maxiter+1))
     betalist = list(length(maxiter+1))
     betalist[[1]] = beta
     beta0 = sum(y-X%*%beta)/(length(y))
     p = exp(beta0*rep(1,length(y)) + X%*%beta)/(1+exp(beta0*rep(1,length(y)) + X%*%beta))
     z = beta0*rep(1,length(y)) + X%*%beta + (y-p)/(p*(1-p))
     omega = p*(1-p)/(sum((p*(1-p))))
     beta0list = numeric(length(maxiter+1))
     beta0 = sum(y-X%*%beta)/(length(y))
     beta0list[1] = beta0
     for (j in 1:maxiter){
       for (k in 1:length(beta)){
         r = z - X[,-k]%*%beta[-k] - beta0*rep(1,length(y))
         beta[k] = (1/sum(omega*X[,k]^2))*soft_thresholding(t(omega*r)%*%X[,k],length(y)*lambda)
       }
       beta0 = sum(y-X%*%beta)/(length(y))
       beta0list[j+1] = beta0
       betalist[[j+1]] = beta
       obj[j] = (1/2)*(1/length(y))*norm(omega*(z - X%*%beta - 
                                                  beta0*rep(1,length(y))),'F')^2 + lambda*sum(abs(beta))
       p = exp(beta0*rep(1,length(y)) + X%*%beta)/(1+exp(beta0*rep(1,length(y)) + X%*%beta))
       z = beta0*rep(1,length(y)) + X%*%beta + (y-p)/(p*(1-p))
       omega = p*(1-p)/(sum((p*(1-p))))
       if (norm(rbind(beta0list[j],betalist[[j]]) - 
                rbind(beta0,beta),'F') < tol) { break } 
     } 
     return(list(obj=obj[1:j],beta=beta,intercept=beta0)) }
   
   df0 = df
   df0$y = as.numeric(df$y)-1
   plot_lambda = function(lambda){
     m = apply(df0,2,mean)
     s = apply(df0,2,sd)
     for(j in 1:2) df0[,j] &lt;- (df0[,j]-m[j])/s[j]
     reg = glmnet(cbind(df0$x1,df0$x2), df0$y==1, alpha=1,lambda=lambda)
     u = seq(0,1,length=101)
     p = function(x,y){
       xt = (x-m[1])/s[1]
       yt = (y-m[2])/s[2]
       predict(reg,newx=cbind(x1=xt,x2=yt),type="response")}
     v = outer(u,u,p)
     image(u,u,v,col=clr10,breaks=(0:10)/10)
     points(df$x1,df$x2,pch=19,cex=1.5,col="white")
     points(df$x1,df$x2,pch=c(1,19)[1+z],cex=1.5)
     contour(u,u,v,levels = .5,add=TRUE)}
   
   reg = glmnet(cbind(df0$x1,df0$x2), df0$y==1, alpha=1)
   par(mfrow=c(1,2))
   plot(reg,xvar="lambda",col=c("blue","red"),lwd=2)
   abline(v=exp(-2.8))
   plot_lambda(exp(-2.8))
   
   reg = glmnet(cbind(df0$x1,df0$x2), df0$y==1, alpha=1)
   par(mfrow=c(1,2))
   plot(reg,xvar="lambda",col=c("blue","red"),lwd=2)
   abline(v=exp(-2.1))
   plot_lambda(exp(-2.1))
   
   m0 = apply(myocarde[myocarde$PRONO=="0",1:7],2,mean)
   m1 = apply(myocarde[myocarde$PRONO=="1",1:7],2,mean)
   Sigma = var(myocarde[,1:7])
   omega = solve(Sigma)%*%(m1-m0)
   omega
   
   x = c(.4,.55,.65,.9,.1,.35,.5,.15,.2,.85)
   y = c(.85,.95,.8,.87,.5,.55,.5,.2,.1,.3)
   z = c(1,1,1,1,1,0,0,1,0,0)
   df = data.frame(x1=x,x2=y,y=as.factor(z))
   m0 = apply(df[df$y=="0",1:2],2,mean)
   m1 = apply(df[df$y=="1",1:2],2,mean)
   Sigma = var(df[,1:2])
   omega = solve(Sigma)%*%(m1-m0)
   omega
   
   library(MASS)
   fit_lda = lda(y ~x1+x2 , data=df)
   fit_lda
   
   b = (t(m1)%*%solve(Sigma)%*%m1-t(m0)%*%solve(Sigma)%*%m0)/2
   
   predlda = function(x,y) predict(fit_lda, data.frame(x1=x,x2=y))$class==1
   vv=outer(vu,vu,predlda)
   contour(vu,vu,vv,add=TRUE,lwd=2,levels = .5)
   
   fit_qda = qda(y ~x1+x2 , data=df)
   plot(df$x1,df$x2,pch=19,
        col=c("blue","red")[1+(df$y=="1")])
   predqda=function(x,y) predict(fit_qda, data.frame(x1=x,x2=y))$class==1
   vv=outer(vu,vu,predlda)
   contour(vu,vu,vv,add=TRUE,lwd=2,levels = .5)
   
   sigmoid = function(x) 1 / (1 + exp(-x))
   
   weights_0 = lm(PRONO~.,data=myocarde)$coefficients
   X = as.matrix(cbind(1,myocarde[,1:7]))
   y_5_1 = sigmoid(X %*% weights_0)
   library(ROCR)
   pred = ROCR::prediction(y_5_1,myocarde$PRONO)
   perf = ROCR::performance(pred,"tpr", "fpr")
   plot(perf,col="blue",lwd=2)
   reg = glm(PRONO~.,data=myocarde,family=binomial(link = "logit"))
   y_0 = predict(reg,type="response")
   pred0 = ROCR::prediction(y_0,myocarde$PRONO)
   perf0 = ROCR::performance(pred0,"tpr", "fpr")
   plot(perf0,add=TRUE,col="red")
   
   loss = function(weights){
     mean( (myocarde$PRONO-sigmoid(X %*% weights))^2) }
   
   weights_1 = optim(weights_0,loss)$par
   y_5_2 = sigmoid(X %*% weights_1)
   pred = ROCR::prediction(y_5_2,myocarde$PRONO)
   perf = ROCR::performance(pred,"tpr", "fpr")
   plot(perf,col="blue",lwd=2)
   plot(perf0,add=TRUE,col="red")
   
   weights_1 = lm(PRONO~1+FRCAR+INCAR+INSYS+PAPUL+PVENT,data=myocarde)$coefficients
   X1 = as.matrix(cbind(1,myocarde[,c("FRCAR","INCAR","INSYS","PAPUL","PVENT")]))
   weights_2 = lm(PRONO~1+INSYS+PRDIA,data=myocarde)$coefficients
   X2 = as.matrix(cbind(1,myocarde[,c("INSYS","PRDIA")]))
   weights_3 = lm(PRONO~1+PAPUL+PVENT+REPUL,data=myocarde)$coefficients
   X3 = as.matrix(cbind(1,myocarde[,c("PAPUL","PVENT","REPUL")]))
   X = cbind(sigmoid(X1 %*% weights_1), sigmoid(X2 %*% weights_2), sigmoid(X3 %*% weights_3))
   weights = c(1/3,1/3,1/3)
   y_5_3 = sigmoid(X %*% weights)
   pred = ROCR::prediction(y_5_3,myocarde$PRONO)
   perf = ROCR::performance(pred,"tpr", "fpr")
   plot(perf,col="blue",lwd=2)
   plot(perf0,add=TRUE,col="red")
   
   center = function(z) (z-mean(z))/sd(z)
   
   loss = function(weights){
     weights_1 = weights[0+(1:7)]
     weights_2 = weights[7+(1:7)]
     weights_3 = weights[14+(1:7)]
     weights_  = weights[21+1:4]
     X1=X2=X3=as.matrix(myocarde[,1:7])
     Z1 = center(X1 %*% weights_1)
     Z2 = center(X2 %*% weights_2)
     Z3 = center(X3 %*% weights_3)
     X = cbind(1,sigmoid(Z1), sigmoid(Z2), sigmoid(Z3))
     mean( (myocarde$PRONO-sigmoid(X %*% weights_))^2)}
   
   pca = princomp(myocarde[,1:7])
   W = get_pca_var(pca)$contrib
   weights_0 = c(W[,1],W[,2],W[,3],c(-1,rep(1,3)/3))
   weights_opt = optim(weights_0,loss)$par
   weights_1 = weights_opt[0+(1:7)]
   weights_2 = weights_opt[7+(1:7)]
   weights_3 = weights_opt[14+(1:7)]
   weights_  = weights_opt[21+1:4]
   X1=X2=X3=as.matrix(myocarde[,1:7])
   Z1 = center(X1 %*% weights_1)
   Z2 = center(X2 %*% weights_2)
   Z3 = center(X3 %*% weights_3)
   X = cbind(1,sigmoid(Z1), sigmoid(Z2), sigmoid(Z3))
   y_5_4 = sigmoid(X %*% weights_)
   pred = ROCR::prediction(y_5_4,myocarde$PRONO)
   perf = ROCR::performance(pred,"tpr", "fpr")
   plot(perf,col="blue",lwd=2)
   plot(perf,add=TRUE,col="red")
   
   library(nnet)
   myocarde_minmax = myocarde
   minmax = function(z) (z-min(z))/(max(z)-min(z))
   for(j in 1:7) myocarde_minmax[,j] = minmax(myocarde_minmax[,j])
   model_nnet = nnet(PRONO~.,data=myocarde_minmax,size=3)
   summary(model_nnet)
   library(devtools)
   source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')
   plot.nnet(model_nnet)
   library(neuralnet)
   model_nnet = neuralnet(formula(glm(PRONO~.,data=myocarde_minmax)),
                          myocarde_minmax,hidden=3, act.fct = sigmoid)
   plot(model_nnet)
   
   library(neuralnet)
   model_nnet = neuralnet(formula(glm(PRONO~.,data=myocarde_minmax)),
                          myocarde_minmax,hidden=3, act.fct = sigmoid)
   plot(model_nnet)
   
   n = length(myocarde[,"PRONO"])
   myocarde0 = myocarde
   myocarde0$PRONO = myocarde$PRONO*2-1
   C = .5
   f = function(param){
     w  = param[1:7]
     b  = param[8]
     xi = param[8+1:nrow(myocarde)]
     .5*sum(w^2) + C*sum(xi)}
   grad_f = function(param){
     w  = param[1:7]
     b  = param[8]
     xi = param[8+1:nrow(myocarde)]
     c(2*w,0,rep(C,length(xi)))}
   
   U = rbind(cbind(myocarde0[,"PRONO"]*as.matrix(myocarde[,1:7]),diag(n),myocarde0[,"PRONO"]),
             cbind(matrix(0,n,7),diag(n,n),matrix(0,n,1)))
   C = c(rep(1,n),rep(0,n))
   
   constrOptim(theta=p_init, f, grad_f, ui = U,ci = C)
   
   library(quadprog)
   eps = 5e-4
   y = myocarde[,"PRONO"]*2-1
   X = as.matrix(cbind(1,myocarde[,1:7]))
   n = length(y)
   D = diag(n+7+1)
   diag(D)[8+0:n] = 0 
   d = matrix(c(rep(0,7),0,rep(C,n)), nrow=n+7+1)
   A = Ui
   b = Ci
   
   sol = solve.QP(D+eps*diag(n+7+1), d, t(A), b, meq=1)
   qpsol = sol$solution
   (omega = qpsol[1:7])
   (b     = qpsol[n+7+1])
   
   y_pred = 2*((as.matrix(myocarde0[,1:7])%*%omega+b)>0)-1
   
   library(quadprog)
   eps = 5e-4
   y = myocarde[,"PRONO"]*2-1
   X = as.matrix(cbind(1,myocarde[,1:7]))
   n = length(y)
   Q = sapply(1:n, function(i) y[i]*t(X)[,i])
   D = t(Q)%*%Q
   d = matrix(1, nrow=n)
   A = rbind(y,diag(n),-diag(n))
   C = .5
   b = c(0,rep(0,n),rep(-C,n))
   sol = solve.QP(D+eps*diag(n), d, t(A), b, meq=1, factorized=FALSE)
   qpsol = sol$solution
   
   omega = apply(qpsol*y*X,2,sum)
   omega
   
   n.samples = nrow(X)
   n.features = ncol(X)
   K = matrix(rep(0, n.samples*n.samples), nrow=n.samples)
   for (i in 1:n.samples){
     for (j in 1:n.samples){
       K[i,j] = X[i,] %*% X[j,] }}
   Dmat = outer(y,y) * K
   Dmat = as.matrix(nearPD(Dmat)$mat) 
   dvec = rep(1, n.samples)
   Amat = rbind(y, diag(n.samples), -1*diag(n.samples))
   bvec = c(0, rep(0, n.samples), rep(-C, n.samples))
   res = solve.QP(Dmat,dvec,t(Amat),bvec=bvec, meq=1)
   a = res$solution 
   bomega = apply(a*y*X,2,sum)
   return(bomega)
   }

linear.kernel = function(x1, x2) {
  return (x1%*%x2)
}
svm.fit = function(X, y, FUN=linear.kernel, C=NULL) {
  n.samples = nrow(X)
  n.features = ncol(X)
  K = matrix(rep(0, n.samples*n.samples), nrow=n.samples)
  for (i in 1:n.samples){
    for (j in 1:n.samples){
      K[i,j] = FUN(X[i,], X[j,])
    }
  }
  Dmat = outer(y,y) * K
  Dmat = as.matrix(nearPD(Dmat)$mat) 
  dvec = rep(1, n.samples)
  Amat = rbind(y, diag(n.samples), -1*diag(n.samples))
  bvec = c(0, rep(0, n.samples), rep(-C, n.samples))
  res = solve.QP(Dmat,dvec,t(Amat),bvec=bvec, meq=1)
  a = res$solution 
  bomega = apply(a*y*X,2,sum)
  return(bomega)
}

VM2 = ksvm(y ~ x1 + x2, data = df0, C=2, kernel = "vanilladot" , prob.model=TRUE, type="C-svc")
pred_SVM2 = function(x,y){
  return(predict(SVM2,newdata=data.frame(x1=x,x2=y), type="probabilities")[,2])}
plot(df$x1,df$x2,pch=c(1,19)[1+(df$y=="1")],
     cex=1.5,xlab="",
     ylab="",xlim=c(0,1),ylim=c(0,1))
vu = seq(-.1,1.1,length=251)
vv = outer(vu,vu,function(x,y) pred_SVM2(x,y))
contour(vu,vu,vv,add=TRUE,lwd=2,levels = .5,col="red")

SVM3 = ksvm(y ~ x1 + x2, data = df0, C=2, kernel = "rbfdot" , prob.model=TRUE, type="C-svc")
pred_SVM3 = function(x,y){
  return(predict(SVM3,newdata=data.frame(x1=x,x2=y), type="probabilities")[,2])}
plot(df$x1,df$x2,pch=c(1,19)[1+(df$y=="1")],
     cex=1.5,xlab="",
     ylab="",xlim=c(0,1),ylim=c(0,1))
vu = seq(-.1,1.1,length=251)
vv = outer(vu,vu,function(x,y) pred_SVM2(x,y))
contour(vu,vu,vv,add=TRUE,lwd=2,levels = .5,col="red")

library(rpart)
cart = rpart(PRONO~.,data=myocarde)
library(rpart.plot)
prp(cart,type=2,extra=1)

gini = function(y,classe){
  T. = table(y,classe)
  nx = apply(T,2,sum)
  n. = sum(T)
  pxy = T/matrix(rep(nx,each=2),nrow=2)
  omega = matrix(rep(nx,each=2),nrow=2)/n
  g. = -sum(omega*pxy*(1-pxy))
  return(g)}

-2*mean(myocarde$PRONO)*(1-mean(myocarde$PRONO))
gini(y=myocarde$PRONO,classe=myocarde$PRONO<Inf)
gini(y=myocarde$PRONO,classe=myocarde[,1]<=100)

entropy = function(y,classe){
  T  = table(y,classe)
  nx = apply(T,2,sum)
  pxy = T/matrix(rep(nx,each=2),nrow=2)
  omega = matrix(rep(nx,each=2),nrow=2)/sum(T)
  g  = sum(omega*pxy*log(pxy))
  return(g)}

mat_gini = mat_v=matrix(NA,7,101)
for(v in 1:7){
  variable=myocarde[,v]
  v_seuil=seq(quantile(myocarde[,v],
                       6/length(myocarde[,v])),
              quantile(myocarde[,v],1-6/length(
                myocarde[,v])),length=101)
  mat_v[v,]=v_seuil
  for(i in 1:101){
    CLASSE=variable<=v_seuil[i]
    mat_gini[v,i]=
      gini(y=myocarde$PRONO,classe=CLASSE)}}
-(gini(y=myocarde$PRONO,classe=(myocarde[,3]<19))-
    gini(y=myocarde$PRONO,classe=(myocarde[,3]<Inf)))/
  gini(y=myocarde$PRONO,classe=(myocarde[,3]<Inf))

idx = which(myocarde$INSYS<19)
mat_gini = mat_v = matrix(NA,7,101)
for(v in 1:7){
  variable = myocarde[idx,v]
  v_seuil = seq(quantile(myocarde[idx,v],
                         7/length(myocarde[idx,v])),
                quantile(myocarde[idx,v],1-7/length(
                  myocarde[idx,v])), length=101)
  mat_v[v,] = v_seuil
  for(i in 1:101){
    CLASSE = variable<=v_seuil[i]
    mat_gini[v,i]=
      gini(y=myocarde$PRONO[idx],classe=CLASSE)}}
par(mfrow=c(3,2))
for(v in 2:7){
  plot(mat_v[v,],mat_gini[v,]) 
}

idx = which(myocarde$INSYS>=19)
mat_gini = mat_v = matrix(NA,7,101)
for(v in 1:7){
  variable=myocarde[idx,v]
  v_seuil=seq(quantile(myocarde[idx,v],
                       6/length(myocarde[idx,v])),
              quantile(myocarde[idx,v],1-6/length(
                myocarde[idx,v])), length=101)
  mat_v[v,]=v_seuil
  for(i in 1:101){
    CLASSE=variable<=v_seuil[i]
    mat_gini[v,i]=
      gini(y=myocarde$PRONO[idx],
           classe=CLASSE)}}
par(mfrow=c(3,2))
for(v in 2:7){
  plot(mat_v[v,],mat_gini[v,])
}

cart = rpart(PRONO~., myocarde)
(split = summary(cart)$splits)

cart$variable.importance

n_iter = 100
y = (myocarde[,"PRONO"]==1)*2-1
x = myocarde[,1:7]
error = rep(0,n_iter) 
f = rep(0,length(y)) 
w = rep(1,length(y)) #
alpha = 1
library(rpart)
for(i in 1:n_iter){
  w = exp(-alpha*y*f) *w 
  w = w/sum(w)
  rfit = rpart(y~., x, w, method="class")
  g = -1 + 2*(predict(rfit,x)[,2]>.5) 
  e = sum(w*(y*g<0))
  alpha = .5*log ( (1-e) / e )
  alpha = 0.1*alpha 
  f = f + alpha*g
  error[i] = mean(1*f*y<0)
}
plot(seq(1,n_iter),error)

set.seed(123)
id_train = sample(1:nrow(myocarde), size=45, replace=FALSE)
train_myocarde = myocarde[id_train,]
test_myocarde = myocarde[-id_train,]

for(i in 1:n_iter){
  w_train = w_train*exp(-alpha*y_train*f_train) 
  w_train = w_train/sum(w_train)
  rfit = rpart(y_train~., x_train, w_train, method="class")
  g_train = -1 + 2*(predict(rfit,x_train)[,2]>.5)
  g_test = -1 + 2*(predict(rfit,x_test)[,2]>.5)
  e_train = sum(w_train*(y_train*g_train<0))
  alpha = .5*log ( (1-e_train) / e_train )
  alpha = 0.1*alpha 
  f_train = f_train + alpha*g_train
  f_test = f_test + alpha*g_test
  train_error[i] = mean(1*f_train*y_train<0)
  test_error[i] = mean(1*f_test*y_test<0)}

library(gbm)
gbmWithCrossValidation = gbm(PRONO ~ .,distribution = "bernoulli",
                             data = myocarde,n.trees = 2000,shrinkage = .01,cv.folds = 5,n.cores = 1)
bestTreeForPrediction = gbm.perf(gbmWithCrossValidation)


   
   
   
   
   
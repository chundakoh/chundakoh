load("data_2.RData")

nts = dt1[,1:5]
kts = dt1[,6:9]
rf = dt1[,10]

nts.mean = colMeans(nts-rf)*12
nts.sd = apply(nts,2,sd) * 12^0.5
nts.sr = nts.mean/ nts.sd

q1 = cbind(nts.mean, nts.sd, nts.sr)
q1

reg_capm = NULL
nts.names = colnames(nts)
for(nts.name in nts.names){
  lm.temp = lm(unlist(nts[nts.name] - rf) ~ kts$Mkt.RF)
  sum = summary(lm.temp)
  ab = sum$coefficients[, "Estimate"]
  ab.se = sum$coefficients[, "Std. Error"]
  R = sum$r.squared
  ab.h0 = c(0,1)
  ab.tstat = (ab - ab.h0) / ab.se
  reg = c(ab, R, ab.tstat)
  reg_capm = rbind(reg_capm, reg)
}

reg_capm[1, 4:5] = NA
rownames(reg_capm) = nts.names
colnames(reg_capm) = c('alpha', 'beta', 'R-squared', 'tstat.alpha', 'tstat.beta')
reg_capm = as.data.frame(reg_capm)
reg_capm$alpha = paste0(round(reg_capm$alpha * 100, 3), '%')
reg_capm

sys = reg_capm$'R-squared' #systematic variance
idio = 1 - sys
var = data.frame(nts.mean, nts.sd, nts.sr, sys, idio)
var


reg_4fac = NULL

for(nts.name in nts.names) {
  lm.temp = lm( unlist(nts[nts.name] - rf)~ kts$Mkt.RF + kts$SMB + kts$HML + kts$MOM)
  sum = summary(lm.temp)
  coef = sum$coefficients[,"Estimate"]
  coef.se = sum$coefficients[,"Std. Error"]
  R = sum$r.squared
  coef.h0 = c(0,1,0,0,0)
  coef.tstat = (coef - coef.h0) / coef.se
  reg = c(coef, R,coef.tstat)
  reg_4fac = rbind(reg_4fac, reg)
}

reg_4fac[,2:5] = round(reg_4fac[,2:5], 3)
reg_4fac[,6:11] = round(reg_4fac[,6:11],2)
reg_4fac[1,7:11] = NA
rownames(reg_4fac) = nts.names
colnames(reg_4fac) = c('alpha', 'Mkt-Rf', 'SMB', 'HML', 'MOM', 'R-squared',
                       'tstat_alpha', 'tstat.Mkt-Rf', 'tstat.SMB', 'tstat.HML','tstat.MOM')
reg_4fac = as.data.frame(reg_4fac)
reg_4fac$alpha = paste0(round(reg_4fac$alpha *100, 3), '%')
reg_4fac  

reg_bkr= NULL
full=rep(TRUE,nrow(dt2)) # indicator of full sample
sub1=substr(rownames(dt2),1,4) %in% as.character(1980:2007)
sub2=!sub1
ind=cbind(full,sub1,sub2)
for(i in 1:ncol(ind)){
  samp=dt2[ind[,i],] #sample
  rf2=samp[,6] #risk free
  BRK=samp[,1] #Berkshire hathway A
  kts2= samp[,2:5] #known trading strategy
  msr=mean(kts2[,1]) / sd(kts2[,1]+ rf2) #market sharpe ratio
  sr=mean(BRK-rf2)/sd(BRK)
  lm.temp=lm(BRK-rf2~kts2$Mkt.RF+kts2$SMB+kts2$HML+kts2$MOM)
  sum=summary(lm.temp)
  coef=sum$coefficients[,"Estimate"] #alpha and beta
  coef.se=sum$coefficients[,"Std. Error"] #se of alpha and beta
  R=sum$r.squared
  coef.h0=0
  coef.tstat=(coef-coef.h0)/coef.se
  reg=c(msr,sr,coef,R,coef.tstat)
  reg_bkr=rbind(reg_bkr,reg)
}

reg_bkr=as.data.frame(reg_bkr)
rownames(reg_bkr)=c('Full','1980-2007','2008-2013')
colnames(reg_bkr)=c('Market Sharpe Ratio','Sharpe Ratio','alpha','rm-rf.beta',
                    'SMB.beta','HML.beta','MOM.beta','R-squared','tstat.alpha'
                    ,'tstat.Mkt-RF','tstat.SMB','tstat.HML','tstat.MOM')
reg_bkr$alpha=paste0(round(reg_bkr$alpha * 100,3),'%')
reg_bkr

library(stringr)
library(Rsolnp)

y=dt3$BRK.A
x1=dt3$Vanguard.S.P.500.Index.Inv..VFINX.
x2=dt3$Vanguard.Small.Cap.Index.Inv..NAESX.
x3=dt3$Vanguard.Value.Index.Inv..VIVAX.

ind_sub = str_sub(rownames(dt3) , -4 , -1) %in% as.character ((2008 : 2003))

#define loss function
f_loss = function(b) { #b is the vector containing parameters
  p=y #portfolio
  b = b[1] + b[2]*x1 + b[3]*x2 + b[4]*x3 # benchmark
  return (sum((b-p)^2))
}

#equality constraints function
eq = function(b){
  return(b[2]+b[3]+b[4])
}

#random start points for parameter optimization
theta= c(0.1,0.5,0.25,0.25)

f=function(dt,sub=FALSE,Shorts=FALSE) {
  if(sub) {
    dt=dt[ind_sub,]
  }
  y=dt[,1]
  x1=dt[,2]
  x2=dt[,3]
  x3=dt[,4]
  #define loss function
  f_loss=function(b) {
    p=y #portfolio
    b = b[1] + b[2]*x1 + b[3]*x2 + b[4]*x3 # benchmark
    return (sum((b-p)^2)) #tracing error
  }
  if (Shorts) {
    m=solnp(theta,f_loss,eqfun=eq,eqB=1,LB=c(-Inf,-Inf,-Inf,-Inf))
  } else {
    m=solnp(theta,f_loss,eqfun=eq,eqB=1,LB=c(-Inf,0,0,0))
  }
  w=m$pars
  y_b=as.matrix(dt) %*% w
  mar=mean(y-y_b) #mean active return
  te=sd(y-y_b) #tracking error
  return(c(w,mar,te))
}

# no shorts unit exposure
r1=f(dt3,sub=F,Shorts=F)
# Shorts Unit Exposure
r2=f(dt3,sub=F,Shorts=T)
# recent no shorts exposure
r3=f(dt3,sub=T,Shorts=F)
# recent shorts unit exposure
r4=f(dt3,sub=T,Shorts=T)

r= as.data.frame(rbind(r1,r2,r3,r4))
rownames(r)=c('1993-2013 No Shorts Unit Exposure',
              '1993-2013 Shorts Unit Exposure',
              '1993-2013 recent no shorts exposure',
              '1993-2013 recent shorts unit exposure')
colnames(r)=c('alpha','Large','Small','Value','Mean active return','Error')
r$alpha=paste0(round(r$alpha * 100,3),'%') #convert to percentage
r[,-1]=round(r[,-1],3) #modify digits
r




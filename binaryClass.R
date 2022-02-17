library(readODS)
library(ROCR)
library(aod)
library(latex2exp)

#hist = NULL
#load('/home/og/Documents/AAA-articles/fowCAM_binary/dataSetBinaryClass.RData')
col.det = 'brown'
col.phy = 'green2'

extract <- function(datafile,taxa,data.date){
  data_raw <- read_ods(datafile,sheet=taxa,skip=1)
  data <- NULL
  
  data$ID <- data_raw$Id # adim
  data$ESD0 <- data_raw$`Diameter (ESD)` # um
  data$ESD <- data_raw$`Diameter (ABD)` # um
  data$AR <- data_raw$`Aspect Ratio` # adim
  data$AREA <- data_raw$`Area (ABD)` # um2
  data$VOL <- data_raw$`Volume (ESD)` # um3
  data$TRANS <- data_raw$Transparency # adim
  data$ELON <- data_raw$Elongation # adim
  data$COMPACT <- data_raw$Compactness # adim
  data$CH1 <- data_raw$`Ch1 Peak` # um
  data$CH2 <- data_raw$`Ch2 Peak` # um
  data$CONP <- data_raw$`Convex Perimeter` # um
  data$PER <- data_raw$Perimeter # um
  data$INT <- data_raw$Intensity	
  data$SIGMAI <- data_raw$`Sigma Intensity`
  data$SUMIN <- data_raw$`Sum Intensity`
  data$ROUGHN <- data_raw$Roughness
  
  rm(data_raw)
  
  data <- sapply(data, as.numeric)
  data <- as.data.frame(data)
  print(taxa)
  str(data)
  
  data$TAXA <- rep(as.character(taxa),times=length(data$ID))
  data$DATE <- rep(as.character(data.date),times=length(data$ID))
  
  data <- data[!is.na(data$ID),]
  return(data)  
}

construct<- function(data.date){
  dataFile <- paste0('./data/orig_data_',data.date,'.ods')
  
  diat <- extract(dataFile,'diat',data.date)
  dinos <- extract(dataFile,'dinos',data.date)
  coccos <- extract(dataFile,'coccos',data.date)
  flag <- extract(dataFile,'flag',data.date)
  #cil <- extract(dataFile,'cil',data.date)
  undef <- extract(dataFile,'undefined',data.date)
  detritus <- extract(dataFile,'detritus',data.date)
  
  total <- diat
  total <- rbind(total,dinos)
  total <- rbind(total,coccos)
  total <- rbind(total,flag)
  total <- rbind(total,undef)
  total <- rbind(total,detritus)
  
  return(as.data.frame(total))
}

construct<- function(data.date){
  dataFile <- paste0('./data/orig_data_',data.date,'.ods')
  taxa.list=list_ods_sheets(dataFile)
  total = NULL
  for(i in taxa.list) total = rbind(total,extract(dataFile,i,data.date))
  return(as.data.frame(total))
}

if(!exists('dd.raw')){
  setwd("/home/og/Documents/mixotrophy/flowcam")
  dataList <- read_ods('./dateList.ods')
  dd.raw = NULL
  for(t in dataList$date){
    este<-construct(t)
    dd.raw = rbind(dd.raw,este)
  }
  if(T)for(t in dataList$noused){
    este<-construct(t)
    dd.raw = rbind(dd.raw,este)
  }
  dd=dd.raw
}

####
####
## Test zone
#construct(data.date = '20070802')


####
####

library(rpart)

eval = function(actual,predicted){
  result = NULL
  kkk = 0*1:length(predicted[,1])
  for(i in 1:length(kkk)) kkk[i]=names(which.max(predicted[i,]))
  for(j in levels(as.factor(actual))){
    P = sum(actual==j)
    N = sum(actual!=j)
    TP = sum((actual==j)&(actual==kkk))
    TN = sum((actual!=j)&(actual==kkk))
    FP = sum((actual!=j)&(actual!=kkk))
    FN = sum((actual==j)&(actual!=kkk))
    sensitivity = TP/P
    specificity = TN/N
    precision = TP/(TP+FP)
    accuracy = (TP+TN)/(P+N)
    FP.rate = FP/(P+N)
    FN.rate = FN/(P+N)
    balanced.acc = (sensitivity+specificity)/2
    result[j]=list(c(sensitivity, specificity,precision,FP.rate,FN.rate,accuracy,balanced.acc))
  }
  return(result)
}

####
####

dd = as.data.frame(dd.raw)
#dd = dd[dd$TAXA!='undefined',]

if(T){
  dd = dd[dd$ROUGHN!=0,]
  dd = dd[dd$COMPACT!=0,]
  dd = dd[dd$SIGMAI!=0,]
  dd = dd[dd$CONP!=0,]
  dd = dd[dd$INT!=0,]
  dd = dd[dd$SUMIN!=0,]
}

dd$PR = (dd$CH1+1e-6)/(dd$CH2+1e-6)
dd$DET = dd$TAXA%in%c('detritus')
#dd$DET = dd$TAXA%in%c('coccos')

names(dd)
levels(as.factor(dd$TAXA))
dd0 = dd

###Panel plot
if(F){
  par(mfrow=c(4,4),las=1,mai=c(.4,.2,.2,.2),oma=c(4,4,0,0),bty='n',font.main=1,adj=0)
  
  plot(density(log(dd$VOL)),lwd=2,ylim=c(0,.8),xlab='',ylab='',main=TeX('log(volume/$\\mu$m$^3$)'),type='n',yaxt='n')
  lines(density(log(dd$VOL[dd$DET])),col=col.det,lwd=3)
  lines(density(log(dd$VOL[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,0.8))
  
  plot(density(log(dd$AREA)),lwd=2,ylim=c(0,.8),xlab='',ylab='',main=TeX('log(area/$\\mu$m$^2$)'),type='n',yaxt='n')
  lines(density(log(dd$AREA[dd$DET])),col=col.det,lwd=3)
  lines(density(log(dd$AREA[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,0.8))
  
  plot(density(log(dd$ESD0)),lwd=2,ylim=c(0,1.7),xlab='',ylab='',main=TeX('log(diameter/$\\mu$m)'),type='n',yaxt='n')
  lines(density(log(dd$ESD0[dd$DET])),col=col.det,lwd=3)
  lines(density(log(dd$ESD0[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,1.7))
  
  plot(density(log(dd$PER)),ylim=c(0,1),xlab='',ylab='',main=TeX('log(perimeter/$\\mu$m)'),type='n',yaxt='n')
  lines(density(log(dd$PER[dd$DET])),col=col.det,lwd=3)
  lines(density(log(dd$PER[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,1))
  
  
  plot(density(log(dd$ELON)),ylim=c(0,.5),xlab='',ylab='',main='log(elongation)',type='n',yaxt='n')
  lines(density(log(dd$ELON[dd$DET])),col=col.det,lwd=3)
  lines(density(log(dd$ELON[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,.5))
  
  plot(density(log(dd$COMPACT)),lwd=2,ylim=c(0,1),xlab='',ylab='',main='log(compactness)',type='n',yaxt='n')
  lines(density(log(dd$COMPACT[dd$DET])),col=col.det,lwd=3)
  lines(density(log(dd$COMPACT[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,1.))
  
  plot(density(log(dd$CONP)),ylim=c(0,1.8),xlab='',ylab='',main='log(convex perimeter)',type='n',yaxt='n')
  lines(density(log(dd$CONP[dd$DET])),col=col.det,lwd=3)
  lines(density(log(dd$CONP[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,1.8))
  
  plot(density((dd$AR)),lwd=2,ylim=c(0,2.5),xlab='',ylab='',main='aspect ratio',type='n',yaxt='n')
  lines(density((dd$AR[dd$DET])),col=col.det,lwd=3)
  lines(density((dd$AR[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,2.5))
  
  
  plot(density(dd$INT),ylim=c(0,0.05),xlab='',ylab='',main='intensity',type='n',yaxt='n')
  lines(density(dd$INT[dd$DET]),col=col.det,lwd=3)
  lines(density(dd$INT[!dd$DET]),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,0.05))
  
  
  plot(density(dd$SIGMAI),ylim=c(0,0.09),xlab='',ylab='',main='sigma intensity',type='n',yaxt='n')
  lines(density(dd$SIGMAI[dd$DET]),col=col.det,lwd=3)
  lines(density(dd$SIGMAI[!dd$DET]),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,0.09))
  
  
  plot(density(log(dd$SUMIN)),ylim=c(0,0.6),xlab='',ylab='',main='sum intensity',type='n',yaxt='n')
  lines(density(log(dd$SUMIN[dd$DET])),col=col.det,lwd=3)
  lines(density(log(dd$SUMIN[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,0.6))
  
  
  plot(density(log(dd$ROUGHN)),lwd=2,ylim=c(0,3),xlab='',ylab='',main='log(roughness)',type='n',yaxt='n')
  lines(density(log(dd$ROUGHN[dd$DET])),col=col.det,lwd=3)
  lines(density(log(dd$ROUGHN[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,3))
  
  
  plot(density(dd$TRANS),ylim=c(0,4.5),xlab='',ylab='',main='transparency',type='n',yaxt='n')
  lines(density(dd$TRANS[dd$DET]),col=col.det,lwd=3)
  lines(density(dd$TRANS[!dd$DET]),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,4.5))
  
  
  plot(density(dd$CH1*1e-3),ylim=c(0,0.8),xlab='',ylab='',main='fluorescence ch. 1',type='n',yaxt='n')
  lines(density(dd$CH1[dd$DET]*1e-3),col=col.det,lwd=3)
  lines(density(dd$CH1[!dd$DET]*1e-3),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,0.8))
  
  
  plot(density(dd$CH2*1e-3),ylim=c(0,1.5),xlab='',ylab='',main='fluorescence ch. 2',type='n',yaxt='n')
  lines(density(dd$CH2[dd$DET]*1e-3),col=col.det,lwd=3)
  lines(density(dd$CH2[!dd$DET]*1e-3),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,1.5))
  
  plot(density(log(dd$PR)),ylim=c(0,1.2),xlab='',ylab='',main='fluorescence peak ratio',type='n',yaxt='n')
  lines(density(log(dd$PR[dd$DET])),col=col.det,lwd=3)
  lines(density(log(dd$PR[!dd$DET])),col=col.phy,lwd=3)
  axis(side=2,labels=F,at=c(0,1.2))
  
  
  mtext('density',side=2,outer=T,las=0,line=1,cex=1.)
  mtext('variable value',side=1,outer=T,las=0,line=1,cex=1.)
  
  barplot(c(sum(dd$DET),sum(!dd$DET))/length(dd$DET),yaxt='n',
          main='fraction of observations',
          xlim=c(-0.25,2.5),
          ylim=c(0,0.6),
          col=c(col.det,col.phy),
          names.arg = c('detritus','phytoplankton'))
  
  axis(side=2,labels=T,at=c(0,0.3,0.6))
  
  
  length(dd$ID)}
####
####
####

par(mfrow=c(2,1),las=1,mai=c(.4,.2,.2,.2),oma=c(4,4,0,0),bty='n',font.main=1,adj=0)

length(dd$ID)
levels(as.factor(dd$TAXA))
#how.many = as.integer(0.01*length(dd$ID))
how.many = as.integer(runif(1,0.01,0.3)*length(dd$ID))
#how.many = as.integer(sample(100:10000,1))

the.filter = 1:length(dd$ID)%in%sample(1:length(dd$ID),how.many)
train.data = as.data.frame(dd)[the.filter,]
test.data = as.data.frame(dd)[!the.filter,]

### TO use rpart models
if(F){
  ### TO identify taxas
  fit=rpart(TAXA~ESD0+AR+AREA+TRANS+ELON+COMPACT+CH1+CH2+CONP+INT+SIGMAI+SUMIN+ROUGHN,data=train.data,
  control = rpart.control(cp=0.01,maxdepth=30))
  fit=rpart(TAXA~ESD0+AR+AREA+TRANS+ELON+COMPACT+CH1+CH2,data=train.data,
            control = rpart.control(cp=0.001,maxdepth=30))
  
  plot(fit)
  text(fit, use.n=T,cex=0.75)
  #print(fit)
  #summary(fit)
  
  ss=predict(fit,test.data)
  tt=predict(fit,train.data)
  rr=predict(fit,dd)
  
  #sensitivity, specificity, precision, FP.rate, FN.rate, accuracy, balanced.acc  
  eval(test.data$TAXA,ss)
  #eval(train.data$TAXA,tt)
  #eval(dd$TAXA,rr)
  
  plot(performance(pred, "tpr", "fpr"),type='n')
  abline(0, 1, lty = 2)
  palette(hcl.colors(8,'Zis'))
  for(i in 1:8){
    which.one = levels(as.factor(dd$TAXA))[i]
    pred <- ROCR::prediction(predict(fit,test.data,type='prob')[,i], test.data$TAXA == which.one)
    plot(performance(pred, "tpr", "fpr"),add=T,col=i,lwd=3,type='l')
  }
  
  legend('bottomright',legend=levels(as.factor(dd$TAXA)),lwd=3,col=1:8,bty='n')
  
  ### TO identify detritus
  
  fit=rpart(as.character(DET)~ESD0+AR+AREA+TRANS+ELON+COMPACT+CH1+CH2+CONP+PER+INT+SIGMAI+SUMIN+ROUGHN+PR,data=train.data,
            control = rpart.control(cp=0.001,maxdepth=30))
  #fit=rpart(as.character(DET)~AR+TRANS+ELON+COMPACT+CH1+CH2+CONP+INT+SIGMAI+SUMIN+ROUGHN,data=train.data,
  #          control = rpart.control(cp=0.001,maxdepth=30))
  fit.rt = fit
  
  plot(fit)
  text(fit, use.n=T, cex=0.75)
  #print(fit)
  summary(fit)
  
  ss=predict(fit,test.data)
  tt=predict(fit,train.data)
  rr=predict(fit,dd)
  
  #sensitivity, specificity, precision, accuracy, balanced.acc  
  eval(test.data$DET,ss)
  #eval(train.data$TAXA,tt)
  #eval(dd$TAXA,rr)
  
  pred <- ROCR::prediction(predict(fit,test.data,type='prob')[,2], test.data$DET)
  plot(performance(pred, "tpr", "fpr"))
  abline(0, 1, lty = 2)
  
  plot(ss[,2],ylim=c(0,1))
  points(test.data$DET,col='red')
}
### logit model

fit = glm(DET~ESD0+AR+AREA+TRANS+ELON+COMPACT+CH1+CH2+CONP+PER+INT+SIGMAI+SUMIN+ROUGHN+PR,data=train.data,family='binomial')
summary(fit)
#fasdfsa

### individual test
if(F){
  ss0=predict(fit,test.data,type = "response")
  plot(ss0,ylim=c(0,1))
  points(test.data$DET,col='red')
  
  ss=NULL
  ss$`FALSE` = 1-ss0
  ss$`TRUE` = ss0
  ss=as.data.frame(ss)
  
  eval(paste0(as.character(test.data$DET),'.'),ss)
  
  pred <- ROCR::prediction(predict(fit,test.data,type='response'), test.data$DET)
  plot(performance(pred, "tpr", "fpr"))
  abline(0, 1, lty = 2)
}

####
####
####
dd = test.data

#the.filter = 1:length(dd0$ID)%in%sample(1:length(dd0$ID),10000)
#dd = as.data.frame(dd0)[the.filter,]
#dd = dd0

the.lim = 0.5
the.bw = 0.1

if(T){
  ss0=predict(fit,dd,type = "response")
  ss=NULL
  ss$`FALSE` = 1-ss0
  ss$`TRUE` = ss0
  ss=as.data.frame(ss)  
  
  kkk=eval(paste0(as.character(dd$DET),'.'),ss)
  hist = rbind(hist,c(how.many,kkk$`TRUE`))

}else{
  ss=predict(fit.rt,dd)
  eval(dd$DET,ss)
  ss=as.data.frame(ss)
  ss$`TRUE.`=ss$`TRUE`
  ss$`FALSE.`=ss$`FALSE`
}

#fasd

ss$TRUE. = ss$TRUE.**1
ss$FALSE. = (1.-ss$TRUE.)

#plot(cumsum(as.numeric(ss$TRUE.>the.lim)),col='red',type='l')
#lines(cumsum(as.numeric(dd$DET)),col='black')

####
####
####

#plot
if(F){
  dd.TP=dd[(ss$TRUE.>the.lim) & dd$DET,]
  dd.FP=dd[(ss$TRUE.>the.lim) & !dd$DET,]
  dd.TN=dd[(ss$FALSE.>the.lim) & !dd$DET,]
  dd.FN=dd[(ss$FALSE.>the.lim) & dd$DET,]
  
  par(mfrow=c(2,2),las=1,mai=c(.4,.2,.2,.2),oma=c(4,4,0,0),bty='n',font.main=1,adj=0)
  
  plot(w<-density(log(dd$ESD),weights=(1+0*dd$VOL)/length(dd$VOL),bw=the.bw),lwd=6,col='gray',
       log='y',
       xlim=log(c(1,300)),
       ylim=c(2e-5,2),xlab='',ylab='',main=TeX('diameter ($\\mu$m) by particle number'),type='l',yaxt='n',xaxt='n')
  lines(density(log(dd.TP$ESD),weights=(1+0*dd.TP$VOL)/length(dd$VOL),bw=w$bw),col=col.det,lwd=3)
  lines(density(log(dd.TN$ESD),weights=(1+0*dd.TN$VOL)/length(dd$VOL),bw=w$bw),col=col.phy,lwd=3)
  #lines(density(log(dd.FP$ESD),weights=(1+0*dd.FP$VOL)/length(dd$VOL),bw=w$bw),col='red',lwd=1,lty=1)
  #lines(density(log(dd.FN$ESD),weights=(1+0*dd.TN$VOL)/length(dd$VOL),bw=w$bw),col='blue',lwd=1,lty=1)
  axis(side=2,labels=c(0,TeX('10$^{-4}$'),TeX('10$^{-2}$'),TeX('10$^{0}$')),at=c(0,0.0002,0.02,2))
  axis(side=2,labels=F,at=2*c(0,1e-4,1e-3,1e-2,1e-1))
  axis(side=1,labels=F,at=log(c(1,3,10,30,100,300)))
  
  plot(density(log(dd$ESD),weights=(1+0*dd$VOL)/length(dd$VOL),bw=the.bw),lwd=6,col='gray',
       log='y',
       xlim=log(c(1,300)),
       ylim=c(2e-5,2),xlab='',ylab='',main='',type='l',yaxt='n',xaxt='n')
  #lines(density(log(dd.TP$ESD),weights=(1+0*dd.TP$VOL)/length(dd$VOL),bw=w$bw),col=col.det,lwd=3)
  #lines(density(log(dd.TN$ESD),weights=(1+0*dd.TN$VOL)/length(dd$VOL),bw=w$bw),col=col.phy,lwd=3)
  lines(density(log(dd.FP$ESD),weights=(1+0*dd.FP$VOL)/length(dd$VOL),bw=w$bw),col='red',lwd=3,lty=1)
  lines(density(log(dd.FN$ESD),weights=(1+0*dd.FN$VOL)/length(dd$VOL),bw=w$bw),col='blue',lwd=3,lty=1)
  axis(side=2,labels=F,at=c(0,2))
  axis(side=2,labels=F,at=2*c(0,1e-4,1e-3,1e-2,1e-1))
  axis(side=1,labels=F,at=log(c(1,3,10,30,100,300)))
  
  plot(w<-density(log(dd$ESD),weights=dd$VOL/sum(dd$VOL),bw=the.bw),lwd=6,col='gray',
       log='y',
       xlim=log(c(1,300)),
       ylim=c(2e-5,2.),xlab='',ylab='',main=TeX('diameter ($\\mu$m) by biomass'),type='l',yaxt='n',xaxt='n')
  lines(density(log(dd.TP$ESD),weights=dd.TP$VOL/sum(dd$VOL),bw=w$bw),col=col.det,lwd=3)
  lines(density(log(dd.TN$ESD),weights=dd.TN$VOL/sum(dd$VOL),bw=w$bw),col=col.phy,lwd=3)
  #lines(density(log(dd.FP$ESD),weights=dd.FP$VOL/sum(dd.FP$VOL),bw=w$bw),col='red',lwd=1,lty=1)
  #lines(density(log(dd.FN$ESD),,weights=dd.FN$VOL/sum(dd.FN$VOL),bw=w$bw),col='blue',lwd=1,lty=1)
  axis(side=2,labels=c(0,TeX('10$^{-4}$'),TeX('10$^{-2}$'),TeX('10$^{0}$')),at=c(0,0.0002,0.02,2))
  axis(side=2,labels=F,at=2*c(0,1e-4,1e-3,1e-2,1e-1))
  axis(side=1,labels=c(1,3,10,30,100,300),at=log(c(1,3,10,30,100,300)))
  
  plot(w<-density(log(dd$ESD),weights=dd$VOL/sum(dd$VOL),bw=the.bw),lwd=6,col='gray',
       log='y',
       xlim=log(c(1,300)),
       ylim=c(2e-5,2.),xlab='',ylab='',main='',type='l',yaxt='n',xaxt='n')
  #lines(density(log(dd.TP$ESD),weights=dd.TP$VOL/sum(dd.TP$VOL),bw=w$bw),col=col.det,lwd=3)
  #lines(density(log(dd.TN$ESD),weights=dd.TN$VOL/sum(dd.TN$VOL),bw=w$bw),col=col.phy,lwd=3)
  lines(density(log(dd.FP$ESD),weights=dd.FP$VOL/sum(dd$VOL),bw=w$bw),col='red',lwd=3,lty=1)
  lines(density(log(dd.FN$ESD),weights=dd.FN$VOL/sum(dd$VOL),bw=w$bw),col='blue',lwd=3,lty=1)
  axis(side=2,labels=F,at=c(0,2))
  axis(side=2,labels=F,at=2*c(0,1e-4,1e-3,1e-2,1e-1))
  axis(side=1,labels=c(1,3,10,30,100,300),at=log(c(1,3,10,30,100,300)))
  
  mtext('density (log-scale)',side=2,outer=T,las=0,line=1,cex=1.)
  mtext(TeX('equivalent spherical diameter, $\\mu$m'),side=1,outer=T,las=0,line=1,cex=1.)
  
  sum(dd.TP$VOL)+sum(dd.FP$VOL)+sum(dd.TN$VOL)+sum(dd.FN$VOL) - sum(dd$VOL)
  length(dd.TP$VOL)+length(dd.FP$VOL)+length(dd.TN$VOL)+length(dd.FN$VOL) - length(dd$VOL)
  length(dd.FP$VOL)/length(dd$VOL)
  length(dd.FN$VOL)/length(dd$VOL)
  
  #pp<-princomp(scale(train.data[c(2:16)]))
  #biplot(pp)
  #plot(pp)
  #summary(pp)
  
  par(mfrow=c(1,2),las=1,mai=c(.4,.2,.2,.2),oma=c(4,4,0,0),bty='n',font.main=1,adj=0)
  
  plot(w<-density(log(dd$ESD),weights=dd$VOL/sum(dd$VOL),bw=the.bw),lwd=6,col='gray',
       log='y',
       xlim=log(c(1,300)),
       ylim=c(1e-5,1.),xlab='',ylab='',main='',type='n',yaxt='n',xaxt='n')
  lines(density(log(dd$ESD[ss$TRUE.>the.lim]),weights=dd$VOL[ss$TRUE.>the.lim]/sum(dd$VOL),bw=w$bw),col=col.det,lwd=3)
  lines(density(log(dd$ESD[dd$DET]),weights=dd$VOL[dd$DET]/sum(dd$VOL),bw=w$bw),col='darkorange',lwd=1,lty=1)
  axis(side=2,labels=c(0,TeX('10$^{-4}$'),TeX('10$^{-2}$'),TeX('10$^{0}$')),at=c(0,0.0001,0.01,1))
  axis(side=2,labels=F,at=1*c(0,1e-4,1e-3,1e-2,1e-1))
  axis(side=1,labels=c(1,3,10,30,100,300),at=log(c(1,3,10,30,100,300)))
  
  plot(w<-density(log(dd$ESD),weights=dd$VOL/sum(dd$VOL),bw=the.bw),lwd=6,col='gray',
       log='y',
       xlim=log(c(1,300)),
       ylim=c(1e-5,1.),xlab='',ylab='',main='',type='n',yaxt='n',xaxt='n')
  lines(density(log(dd$ESD[ss$FALSE.>the.lim]),weights=dd$VOL[ss$FALSE.>the.lim]/sum(dd$VOL),bw=w$bw),col=col.phy,lwd=3)
  lines(density(log(dd$ESD[!dd$DET]),weights=dd$VOL[!dd$DET]/sum(dd$VOL),bw=w$bw),col='DARKGREEN',lwd=1,lty=1)
  axis(side=2,labels=F,at=1*c(0,1e-4,1e-3,1e-2,1e-1,1))
  axis(side=1,labels=c(1,3,10,30,100,300),at=log(c(1,3,10,30,100,300)))
  
  
  mtext('density (log-scale)',side=2,outer=T,las=0,line=1,cex=1.)
  mtext(TeX('equivalent spherical diameter, $\\mu$m'),side=1,outer=T,las=0,line=1,cex=1.)
  
  ####
  ####
  ####
  
  par(mfrow=c(1,2),las=1,mai=c(.4,.2,.2,.2),oma=c(4,4,0,0),bty='n',font.main=1,adj=0)
  
  plot(w<-density(log(dd$ESD),bw=the.bw),lwd=6,col='gray',
       log='y',
       xlim=log(c(1,300)),
       ylim=c(1e-5,1.),xlab='',ylab='',main='',type='n',yaxt='n',xaxt='n')
  lines(density(log(dd$ESD[ss$TRUE.>the.lim]),bw=w$bw),col=col.det,lwd=4)
  lines(density(log(dd$ESD[dd$DET]),bw=w$bw),col='darkorange',lwd=1.5,lty=1)
  axis(side=2,labels=c(0,TeX('10$^{-4}$'),TeX('10$^{-2}$'),TeX('10$^{0}$')),at=c(0,0.0001,0.01,1))
  axis(side=2,labels=F,at=1*c(0,1e-4,1e-3,1e-2,1e-1))
  axis(side=1,labels=c(1,3,10,30,100,300),at=log(c(1,3,10,30,100,300)))
  
  plot(w<-density(log(dd$ESD),bw=the.bw),lwd=6,col='gray',
       log='y',
       xlim=log(c(1,300)),
       ylim=c(1e-5,1.),xlab='',ylab='',main='',type='n',yaxt='n',xaxt='n')
  lines(density(log(dd$ESD[ss$FALSE.>the.lim]),bw=w$bw),col=col.phy,lwd=4)
  lines(density(log(dd$ESD[!dd$DET]),bw=w$bw),col='DARKGREEN',lwd=1.5,lty=1)
  axis(side=2,labels=F,at=1*c(0,1e-4,1e-3,1e-2,1e-1,1))
  axis(side=1,labels=c(1,3,10,30,100,300),at=log(c(1,3,10,30,100,300)))
  
  mtext('density (log-scale)',side=2,outer=T,las=0,line=1,cex=1.)
  mtext(TeX('equivalent spherical diameter, $\\mu$m'),side=1,outer=T,las=0,line=1,cex=1.)
}
####
####
####

#sensitivity, specificity, precision, FP.rate, FN.rate, accuracy, balanced.acc
par(mfrow=c(2,2),las=1,mai=c(.4,.2,.2,.2),oma=c(4,4,2,0),bty='n',font.main=1,adj=0)
palette(hcl.colors(8,'Spec'))

plot(hist[,1],hist[,2],log='',ylim=c(0.8,0.9),xlim=c(2,8000),type='n',xaxt='n')
axis(side=1,at=c(00,2000,4000,6000,8000),labels = T)
axis(side=3,at=c(1,5,10,20,30)/100*24578,labels=c('1','5',10,20,30))
mtext('percentage of observations, %',side=3,outer=F,las=0,line=2,cex=1.)

#for(i in c(2:6,8)) points(hist[,1],hist[,i],bg=i,pch=21,cex=1.,lwd=0.5)
for(i in c(2:6,8)) lines(supsmu(hist[,1],hist[,i]),col='black',lwd=5)
for(i in c(2:6,8)) lines(supsmu(hist[,1],hist[,i]),col=i,lwd=3)
for(i in c(2:6,8)) print(c(mean(hist[,i]),sd(hist[,i])))
points(hist[,1],hist[,8],bg=8,pch=21,cex=1.,lwd=0.5)

mtext('observations using for fitting',side=1,outer=F,las=0,line=3,cex=1.)
mtext('fraction of observations',side=2,outer=F,las=0,line=3,cex=1.)



#plot.new()
#plot(hist[,1],hist[,2],log='x',ylim=c(0.05,0.1),xlim=c(200,10000),type='n',xaxt='n')
#axis(side=1,at=c(100,300,1000,3000,10000))
##for(i in 2:8) points(hist[,1],hist[,i],bg=i,pch=21,cex=1.5)
#for(i in 2:8) lines(supsmu(hist[,1],hist[,i]),col=i,lwd=3)
kkk<-hist(ss$TRUE.,col=col.det)
hist(ss$FALSE.,add=T,col=col.phy)
plot(kkk$breaks,c(cumsum(kkk$counts)/sum(kkk$counts),1),type='s')
points(0.5,0.5)

### Wald test
if(F){
  k=order(dd$ID)
  plot(ss[,2][k],ylim=c(0,1))
  points(dd$DET[k],col='red')
  
  summary(fit)
  for(i in 1:16)print(aod::wald.test(b=coef(fit),Sigma=vcov(fit),Terms=i,verbose=F))
  aod::wald.test(b=coef(fit),Sigma=vcov(fit),Terms=1:16,verbose=T)
}

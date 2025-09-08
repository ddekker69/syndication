library(statnet)
library(relevent)
library(tcltk2)
library(dplyr)
library(lubridate)
getwd()
dn<-gsub("/","\\\\",tclvalue(tkchooseDirectory() ))
setwd(dn)

funding_investor_el = read.csv('C:/Users/david/Dropbox/EBS-HW/Research/venture capital/syndication/funding_investor_edgelist.csv')
names(funding_investor_el)
el<-funding_investor_el[,2:4]
nodes<-unique(c(el$Sender, el$Receiver))
el['Sender']<-match(el$Sender, nodes)
el['Receiver']<-match(el$Receiver, nodes)
names(el)<-c("time", "src", "dest")
dte<-as.Date(el$time)
dttm <- seconds_to_period(seq(from=3600,by=24, length.out=length(el$time)))
dttm<-sprintf("%02i:%02i:%02i", hour(dttm), minute(dttm), second(dttm))
el$time <- as.numeric(as.POSIXct(paste(dte,dttm), format="%Y-%m-%d %H:%M:%S", tz='CET'))
el<-el[sort.list(el$time),]
attr(el, "n")<-length(nodes)
fiel<-as.sociomatrix.eventlist(el, length(nodes))

##Generate some simple sample data based on fixed effects 
roweff<-rnorm(10)#Buildratematrix 
roweff<-roweff-roweff[1]#Adjustforlaterconvenience 
coleff<-rnorm(10) 
coleff<-coleff-coleff[1]
lambda<-exp(outer(roweff,coleff,"+")) 
diag(lambda)<-0 
ratesum<-sum(lambda) 
esnd<-as.vector(row(lambda))#Listofsenders/receivers 
erec<-as.vector(col(lambda)) 
time<-0 
edgelist<-vector() 
while(time<15){#Observethesystemfor15timeunits 
  drawsr<-sample(1:100,1,prob=as.vector(lambda))#Drawfrommodel 
  time<-time+rexp(1,ratesum) 
  if(time<=15)#Censorat15 
    edgelist<-rbind(edgelist,c(time,esnd[drawsr],erec[drawsr])) 
  else edgelist<-rbind(edgelist,c(15,NA,NA)) 
  }

#Fit the model, ordinal BPM 
effects<-c("NTDegSnd","NTDegRec", "OSPSnd", "ISPSnd", "OTPSnd", "ITPSnd", "RRecSnd", "RSndSnd")
#covar<-list((funding_investor_el$Is.Deadpooled=="Yes")*1)

fit.ord<-rem.dyad(el,length(nodes), effects=effects ,hessian=TRUE) 
summary(fit.ord) 
par(mfrow=c(1,2))#Checkthecoefficients 
plot(roweff[-1],fit.ord$coef,asp=1) 
abline(0,1) 
plot(coleff[-1],fit.ord$coef[10:18],asp=1) 
abline(0,1)

#Now,findthetemporalBPM 
fit.time<-rem.dyad(el,length(nodes),effects=effects, ordinal=FALSE,hessian=TRUE) 
summary(fit.time) 
plot(fit.ord$coef,fit.time$coef,asp=1)#Similarresults 
abline(0,1) 

#Finally,trytheBSIRmethod(note:amuchlargerexpansionfactor #isrecommendedinpractice) 
fit.bsir<-rem.dyad(edgelist,length(nodes),effects=effects,fit.method="BSIR", sir.draws=100,sir.expand=5) 
summary(fit.bsir) 
par(mfrow=c(3,3))#Examinetheapproximateposteriormarginals 
for(i in 1:9){ 
  hist(fit.bsir$post[,i],main=names(fit.bsir$coef)[i],prob=TRUE) 
  abline(v=roweff[i+1],col=2,lwd=3) 
  } 
for(i in 10:18){ 
  hist(fit.bsir$post[,i],main=names(fit.bsir$coef)[i],prob=TRUE) 
  abline(v=coleff[i-8],col=2,lwd=3) 
  } 

#Simulateaneventsequencefromthetemporalmodel 
sim<-simulate(fit.time,nsim=50000)#Simulate50000events 
head(sim)#Showtheeventlist 
par(mfrow=c(1,2))#Checkthebehavior 
esnd<-exp(c(0,fit.time$coef[1:8])) 
esnd<-esnd/sum(esnd)*5e4#Expectedsendingcount
erec<-exp(c(0,fit.time$coef[10:18])) 
erec<-erec/sum(erec)*5e4#Expectedsendingcount 
plot(esnd,tabulate(sim[,2]),xlab="ExpectedOut-events",ylab="Out-events") 
abline(0,1,col=2) 
plot(erec,tabulate(sim[,3]),xlab="ExpectedIn-events",ylab="In-events") 
abline(0,1,col=2) 

#Keepthefirst10eventsofthesimulatedsequence,andproduce10more 
sim.pre<-sim[1:10,] 
sim2<-simulate(fit.time,nsim=20,edgelist=sim.pre) 
sim.pre#Seethefirst10events 
sim2#First10eventspreserved 
all(sim2[1:10,]==sim.pre)#AllTRUE 

#Repeat,butredrawingpartoftheinputsequence 
sim2.t<-simulate(fit.time,nsim=20,edgelist=sim.pre,redraw.timing=TRUE) 
sim2.e<-simulate(fit.time,nsim=20,edgelist=sim.pre,redraw.events=TRUE) 
sim2.t#Eventskept,timingsnot 
sim2.t[1:10,]==sim.pre#SecondtwocolumnsTRUE 
sim2.e#Timingkept,eventsnot 
sim2.e[1:10,]==sim.pre#(Note:someeventsmayrepeatbychance!)

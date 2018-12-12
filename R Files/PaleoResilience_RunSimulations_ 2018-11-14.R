#######################################
# Paleoecological Resilience Indicator functions
# Stegner et al. 
# updated 08 July 2018
#######################################

source("Grass-Wood model 30June2018.R")
source("PaleoResilience functions 07July2018.R")

# Set Grass-Wood model parameters________________________
h = 0.5
r=0.25
delta_t = 1
gens = 10000
K_Start = 1
K_Pulse_amt = -0.4
pulse_time = 1000
sigma_sd = 0.005
V0 = 1
beta_ps<-estBetaParams(mu=0.15, var=0.015)
phi = 0.05

# Set taphonomic parameters______________________________
exRStime=6000
nreps=5
nsamp=200
AC.buff=0.1
samp.freq2=0.4
steps<-c(1,1,1,1)
a0=2
a1=0.025
sd.pct=0.05
AC.samp=0.4

# Run time series iterations________________________________
# linear (5yr/cm) 
TAtop=5
TAbottom=5
agemodel="linearTA"
windows<-c(2500,600,50,50)

# run simulations for gradually-forced critical transitions
system.time(Xct1<-rep.ews(TStype="TSct",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0,a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

# run simulations for gradually-forced non-critical transitions
system.time(Xdc1<-rep.ews(TStype="TSdc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=1,sd.pct=sd.pct,AC.samp=AC.samp))

# run simulations for abruptly-forced critical transitions
system.time(Xrs1<-rep.ews(TStype="TSrs",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

# run simulations for no change scenario
system.time(Xnc1<-rep.ews(TStype="TSnc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

# linear (20yr/cm) #
TAtop=20
TAbottom=20
agemodel="linearTA"
windows<-c(2500,150,50,50)

system.time(Xct2<-rep.ews(TStype="TSct",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0,a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xdc2<-rep.ews(TStype="TSdc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=1,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xrs2<-rep.ews(TStype="TSrs",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xnc2<-rep.ews(TStype="TSnc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

# broken stick 1 #
TAtop=5
TAbottom=20
breakpoint=2500
agemodel="brokenstick"
windows<-c(2500,400,50,50)

system.time(Xct3<-rep.ews(TStype="TSct",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0,a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xdc3<-rep.ews(TStype="TSdc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=1,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xrs3<-rep.ews(TStype="TSrs",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xnc3<-rep.ews(TStype="TSnc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

# broken stick 2 #
TAtop=5
TAbottom=20
breakpoint=4000
agemodel="brokenstick"
windows<-c(2500,300,50,50)

system.time(Xct4<-rep.ews(TStype="TSct",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0,a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xdc4<-rep.ews(TStype="TSdc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=1,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xrs4<-rep.ews(TStype="TSrs",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xnc4<-rep.ews(TStype="TSnc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

# Plot Figure 4 (SD for each age model)_____________________________________
# settings for Figures 4 and 5
colorCT<-rgb(0.1,0.3,0.4,1)
colorDC<-rgb(0.1,0.3,0.4,0.7)
colorRS<-rgb(0.1,0.3,0.4,0.5)
colorNC<-rgb(0.1,0.3,0.4,0.3)
mains2<-c("","","","")

ind="sd"
dev.new(width=9,height=5)
par(oma=c(6,4,3,1),mar=c(0.75,0.5,0,0.5))
nf<-layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),nrow=4,ncol=5,byrow=F))

letters.list<-c("a)","b)","c)","d)")
plot.taph.hists(Xct1,Xdc1,Xrs1,Xnc1,indicator=ind,yaxis=T,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Untransformed",type.label=F,taph.ind=1)
 
letters.list<-c("e)","f)","g)","h)")
plot.taph.hists(Xct1,Xdc1,Xrs1,Xnc1,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Linear 5 yrs/cm",type.label=F,taph.ind=2)

letters.list<-c("i)","j)","k)","l)")
plot.taph.hists(Xct2,Xdc2,Xrs2,Xnc2,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Linear 20 yrs/cm",type.label=F,taph.ind=2)

letters.list<-c("m)","n)","o)","p)")
plot.taph.hists(Xct3,Xdc3,Xrs3,Xnc3,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Broken Stick (2500)",type.label=F,taph.ind=2)

letters.list<-c("q)","r)","s)","t)")
plot.taph.hists(Xct4,Xdc4,Xrs4,Xnc4,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Broken Stick (4000)",type.label=T,taph.ind=2)

mtext("Frequency",2,line=2,outer=T,cex=1.2) 
mtext("Kendall's tau",1,line=2.5,outer=T,cex=1.2) 

# Plot Figure 5 (AC for each age model)_____________________________________
# settings for Figures 4 and 5
colorCT<-rgb(0.1,0.3,0.4,1)
colorDC<-rgb(0.1,0.3,0.4,0.7)
colorRS<-rgb(0.1,0.3,0.4,0.5)
colorNC<-rgb(0.1,0.3,0.4,0.3)
mains2<-c("","","","")

ind<-"ac"
dev.new(width=9,height=5)
par(oma=c(6,4,3,1),mar=c(0.75,0.5,0,0.5))
nf<-layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),nrow=4,ncol=5,byrow=F))

# Plot col 1: untransformed time series
letters.list<-c("a)","b)","c)","d)")
plot.taph.hists(Xct1,Xdc1,Xrs1,Xnc1,indicator=ind,yaxis=T,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Untransformed",type.label=F,taph.ind=1)

letters.list<-c("e)","f)","g)","h)")
plot.taph.hists(Xct1,Xdc1,Xrs1,Xnc1,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Linear 5 yrs/cm",type.label=F,taph.ind=2)

letters.list<-c("i)","j)","k)","l)")
plot.taph.hists(Xct2,Xdc2,Xrs2,Xnc2,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Linear 20 yrs/cm",type.label=F,taph.ind=2)

letters.list<-c("m)","n)","o)","p)")
plot.taph.hists(Xct3,Xdc3,Xrs3,Xnc3,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Broken Stick (2500)",type.label=F,taph.ind=2)

letters.list<-c("q)","r)","s)","t)")
plot.taph.hists(Xct4,Xdc4,Xrs4,Xnc4,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Broken Stick (4000)",type.label=T,taph.ind=2)

mtext("Frequency",2,line=2,outer=T,cex=1.2) 
mtext("Kendall's tau",1,line=2.5,outer=T,cex=1.2) 

# Plot Figure 6 (SD for age models and subsampling)____________________________________________
ind<-"sd"
mains2<-c("","","","")
ymax<-1.05

dev.new(width=10,height=6)
par(oma=c(6,4,4,1),mar=c(0.75,0.5,0,0.5))
nf<-layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),nrow=4,ncol=6,byrow=F))

colorCT<-rgb(0,0,1,1)
colorDC<-rgb(0,0,1,0.7)
colorRS<-rgb(0,0,1,0.4)
colorNC<-rgb(0,0,1,0.2)

title2<-"Age-Depth"
letters.list<-c("a)","b)","c)","d)")
plot.taph.hists(Xct2,Xdc2,Xrs2,Xnc2,indicator=ind,yaxis=T,mains=mains2,ymax=ymax,labs2=NULL,letters=letters.list,title=title2,type.label=F,taph.ind=2)

letters.list<-c("e)","f)","g)","h)")
title2<-"AD+Even Samp"
plot.taph.hists(Xct2,Xdc2,Xrs2,Xnc2,indicator=ind,yaxis=F,mains=mains2,ymax=ymax,labs2=NULL,letters=letters.list,title=title2,type.label=F,taph.ind=3)

letters.list<-c("i)","j)","k)","l)")
title2<-"AD+Targeted Samp"
plot.taph.hists(Xct2,Xdc2,Xrs2,Xnc2,indicator=ind,yaxis=F,mains=mains2,ymax=ymax,labs2=NULL,letters=letters.list,title=title2,type.label=F,taph.ind=4)

mtext("Linear 20 yrs/cm",3,line=2.25,outer=T,cex=1.2,adj=0.2)
mtext("Frequency",2,line=2,outer=T,cex=1.2)  
mtext("Kendall's tau",1,line=2.25,outer=T,cex=1.2)  
  
colorCT<-rgb(1,0,0,1)
colorDC<-rgb(1,0,0,0.7)
colorRS<-rgb(1,0,0,0.4)
colorNC<-rgb(1,0,0,0.2)

title2<-"Age-Depth"
letters.list<-c("m)","n)","o)","p)")
plot.taph.hists(Xct4,Xdc4,Xrs4,Xnc4,indicator=ind,yaxis=F,mains=mains2,ymax=ymax,labs2=NULL,letters=letters.list,title=title2,type.label=F,taph.ind=2) 

title2<-"AD+Even Samp"
letters.list<-c("q)","r)","s)","t)")
plot.taph.hists(Xct4,Xdc4,Xrs4,Xnc4,indicator=ind,yaxis=F,mains=mains2,ymax=ymax,labs2=NULL,letters=letters.list,title=title2,type.label=F,taph.ind=3)

title2<-"AD+Targeted Samp"
letters.list<-c("u)","v)","w)","x)")
plot.taph.hists(Xct4,Xdc4,Xrs4,Xnc4,indicator=ind,yaxis=F,mains=mains2,ymax=ymax,labs2=NULL,letters=letters.list,title=title2,type.label=T,taph.ind=4)
mtext("Broken Stick (4000)",3,line=2,outer=T,cex=1.1,adj=0.8)


# Generate Supplemental_____________________________________________
# create table
Supp1<-rbind(Kt.summary.stats(Xct1,Xdc1,Xrs1,Xnc1,3,"sd"),
	Kt.summary.stats(Xct2,Xdc2,Xrs2,Xnc2,3,"sd"),
	Kt.summary.stats(Xct3,Xdc3,Xrs3,Xnc3,3,"sd"),
	Kt.summary.stats(Xct4,Xdc4,Xrs4,Xnc4,3,"sd"))

Supp2<-rbind(Kt.summary.stats(Xct1,Xdc1,Xrs1,Xnc1,3,"ac"),
	Kt.summary.stats(Xct2,Xdc2,Xrs2,Xnc2,3,"ac"),
	Kt.summary.stats(Xct3,Xdc3,Xrs3,Xnc3,3,"ac"),
	Kt.summary.stats(Xct4,Xdc4,Xrs4,Xnc4,3,"ac"))

supp1s<-Supp1[,c(1,2,3,5,6,8,9,10)]
supp2s<-Supp2[,c(1,2,3,5,6,8,9,10)]

# Plot Supplemental
plot.supp(Xct1,Xdc1,Xrs1,Xnc1,"sd","Linear 5 yr/cm, Standard Deviation")
plot.supp(Xct2,Xdc2,Xrs2,Xnc2,"sd","Linear 20 yr/cm, Standard Deviation")
plot.supp(Xct3,Xdc3,Xrs3,Xnc3,"sd","Broken Stick (2500), Standard Deviation")
plot.supp(Xct4,Xdc4,Xrs4,Xnc4,"sd","Broken Stick (4000), Standard Deviation")

plot.supp(Xct1,Xdc1,Xrs1,Xnc1,"ac","Linear 5 yr/cm, Autocorrelation Time")
plot.supp(Xct2,Xdc2,Xrs2,Xnc2,"ac","Linear 20 yr/cm, Autocorrelation Time")
plot.supp(Xct3,Xdc3,Xrs3,Xnc3,"ac","Broken Stick (2500), Autocorrelation Time")
plot.supp(Xct4,Xdc4,Xrs4,Xnc4,"ac","Broken Stick (4000), Autocorrelation Time")


# Set up for plotting figures 2, 3___________________________________
# single runs of GW model, all 4 time series types
single_test = single_run(r=r, gens=gens, delta_t=delta_t, K_Start=K_Start, K_Pulse_amt=K_Pulse_amt, V0=V0, pulse_time=pulse_time,driver_press_topo="gradual",q=5)
TS<-single_test[,3]
TSct<-cbind(c(1:length(TS)),TS)
Kct<-single_test[,2]
bp.out<-CE.Normal.Mean(as.data.frame(TS),Nmax=1)
timeCT<-bp.out$BP.Loc	

single_test = single_run(r=r, gens=gens, delta_t=delta_t, K_Start=K_Start, K_Pulse_amt=K_Pulse_amt, V0=V0, pulse_time=pulse_time,driver_press_topo="gradual",q=1)
TS<-single_test[,3]
TSdc<-cbind(c(1:length(TS)),TS)
Kdc<-single_test[,2]			

single_test = single_run(r=r, gens=gens, delta_t=delta_t, K_Start=K_Start, K_Pulse_amt=K_Pulse_amt, V0=V0, pulse_time=exRStime,driver_press_topo="abrupt",q=5)
TS<-single_test[,3]
TSrs<-cbind(c(1:length(TS)),TS)
Krs<-single_test[,2]
			
single_test = single_run(r=r, gens=gens, delta_t=delta_t, K_Start=K_Start, K_Pulse_amt=0, V0=V0, pulse_time=50,driver_press_topo="gradual",q=5)
TS<-single_test[,3]
TSnc<-cbind(c(1:length(TS)),TS)
Knc<-single_test[,2]

#Age models
# single runs of GW model, all 4 time series types
single_test = single_run(r=0.25, gens=300, delta_t=1, K_Start=1, K_Pulse_amt=-0.4, V0=V0, pulse_time=3,driver_press_topo="gradual",q=5)
TS<-single_test[,3]
TSctexp<-cbind(c(1:length(TS)),TS)
CTexp<-trimtoRS2(TSctexp,cutoff=150,cutoff2=50,trim.type="to.RS",start=60)
oTSexp<-CTexp$trimTS
tailexp<-CTexp$tail.length
am1<-Carpenter_timeavgTS(oTSexp,a0,a1,tailexp)

XX<-trimtoRS2(TSct,cutoff=5000,cutoff2=3000,trim.type="to.RS",start=1000)
oTS<-XX$trimTS
tail<-XX$tail.length
am2<-timeavgTS(oTS,TAbottom=5,TAtop=5,sd.pct=sd.pct)
am3<-timeavgTS(oTS,TAbottom=20,TAtop=20,sd.pct=sd.pct)
am4<-brokenstick.timeavgTS(oTS,TAbottom=20,TAtop=5,breakpoint=2500,sd.pct=sd.pct)
am5<-brokenstick.timeavgTS(oTS,TAbottom=20,TAtop=5,breakpoint=4000,sd.pct=sd.pct)

adTS<-am3$adTS
regTS<-sampleTS(adTS,sample.by="distribute.samples",samp.freq1=NULL,nsamp,timeCT)
acTS<-sampleTSatAC(adTS,AC.buffer=0.10,AC.samp=0.4,timeCT,XX$tail.length)

# Plot Figure 2__________________________________________________
dev.new(width=10,height=6)
X<-c(1,1,1,2,2,2,3,3,3,4,4,4,
1,1,1,2,2,2,3,3,3,4,4,4,
5,5,5,5,5,6,6,7,7,7,7,7,
8,8,8,8,8,9,9,10,10,10,10,10,11)
matx<-matrix(X,nrow=12,ncol=4)
nf<-layout(matx)
layout.show(nf)
par(mar=c(0,5,0.25,0.5),oma=c(6,2,3,2),mgp=c(2.4,1,0))

# plot 4 time series types with CPT results
plot(TSnc,type="l",ylab="",ylim=c(-0.05,1.1),xaxt="n",xlab="",las=1,col="gray40",lwd=1.5,xlim=c(-300,10000),yaxt="n")
text(-350,1,"a)",cex=1.75)
axis(2,at=seq(0,1,0.5),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
tempK<-cbind(1:gens,Knc)
int<-seq(1,gens,275)
points(tempK[int,],pch="-",col="red",cex=2)
#lines(1:gens,Knc,col="red",lty=6,lwd=2)
legend("bottomright",c("Tree cover","K parameter"),lty=c(1,2),col=c("gray40","red"),bty="n")

plot(TSrs,type="l",ylab="",ylim=c(-0.05,1.1),xaxt="n",xlab="",las=1,col="gray40",lwd=1.5,xlim=c(-300,10000),yaxt="n")
text(-350,1,"b)",cex=1.75)
axis(2,at=seq(0,1,0.5),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
tempK<-cbind(1:gens,Krs)
int<-seq(1,gens,275)
points(tempK[int,],pch="-",col="red",cex=2)
#lines(1:gens,Krs,col="red",lty=3,lwd=2)

plot(TSdc,type="l",ylab="",ylim=c(-0.05,1.1),xaxt="n",xlab="",las=1,col="gray40",lwd=1.5,xlim=c(-300,10000),yaxt="n")
text(-350,1,"c)",cex=1.75)
axis(2,at=seq(0,1,0.5),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
tempK<-cbind(1:gens,Kdc)
int<-seq(1,gens,275)
points(tempK[int,],pch="-",col="red",cex=2)
#lines(1:gens,Kdc,col="red",lty=3,lwd=2)

plot(TSct,type="l",ylab="",xlab="time steps",ylim=c(-0.05,1.1),las=1,col="gray40",lwd=1.5,xlim=c(-300,10000),yaxt="n",xaxt="n")
text(-350,1,"d)",cex=1.75)
axis(2,at=seq(0,1,0.5),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
axis(1,hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
mtext("Time Steps (Years)",1,line=2,cex=1)
tempK<-cbind(1:gens,Kct)
int<-seq(1,gens,275)
points(tempK[int,],pch="-",col="red",cex=2)
#lines(1:gens,Kct,col="red",lty=3,lwd=2)

mtext("Tree cover (proportional) (gray)   ",2,line=-1,outer=T,adj=0.2)
mtext("K parameter (red)",2,line=-1,outer=T,col="red",adj=0.85)

# Plot age models
plot(am2$adTS[,1],am2$TAbins[1:(length(am2$TAbins)-1)],type="l",ylim=c(0,26),ylab="",xlab="Time steps",lwd=1.5,col="black",lty=1,xlim=c(0,7000),xaxt="n",las=1,yaxt="n")
axis(2,at=seq(0,35,10),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
axis(1,hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
lines(am3$adTS[,1],am3$TAbins[1:(length(am3$TAbins)-1)],lwd=1.5,lty=1,col="gray40")
text(300,24.5,"e)",cex=1.75)
mtext("Time averaging (yrs/cm)",2,line=2.25,cex=1,outer=F)
mtext("Time Steps (Years)",1,line=2.5,cex=1,outer=F)

# plot spacer
plot(c(0,1),c(0,1),type="n",bty="n",xaxt="n",yaxt="n",ylab="",xlab="")

plot(am4$adTS[,1],am4$TAbins[1:(length(am4$TAbins)-1)],type="l",ylim=c(0,26),ylab="",xlab="Time steps",lwd=1.5,col="black",lty=1,xlim=c(0,7000),xaxt="n",las=1,yaxt="n")
axis(2,at=seq(0,35,10),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
axis(1,hadj=0.6,las=1,tcl=-0.25,,cex.axis=1.5)
text(300,24.5,"f)",cex=1.75)
mtext("Time averaging (yrs/cm)",2,line=2.25,cex=1,outer=F)
mtext("Time Steps (Years)",1,line=2.5,cex=1,outer=F)

plot(am5$adTS[,1],am5$TAbins[1:(length(am5$TAbins)-1)],type="l",ylim=c(0,26),ylab="",xlab="Time steps",lwd=1.5,col="black",lty=1,xlim=c(0,7000),xaxt="n",las=1,yaxt="n")
axis(2,at=seq(0,35,10),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
axis(1,at=seq(1000,7000,2000),labels=seq(2000,8000,2000),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
text(300,24.5,"g)",cex=1.75)
mtext("Time Steps (Years)",1,line=2.5,cex=1,outer=F)
mtext("Time averaging (yrs/cm)",2,line=2.25,cex=1,outer=F)

# plot spacer
plot(c(0,1),c(0,1),type="n",bty="n",xaxt="n",yaxt="n",ylab="",xlab="")

plot(am1$adTS[,1],am1$TAvect,type="l",xlim=c(330,75),las=1,ylab="",xlab="",lwd=1.5,xaxt="n",yaxt="n")
axis(2,at=seq(0,8,2),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
axis(1,at=seq(100,250,50),labels=seq(100,250,50),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
text(340,6.6,"h)",cex=1.75,pos=4)
mtext("Time averaging (yrs/cm)",2,line=2.25,cex=1,outer=F)
mtext("Time Steps (Years)",1,line=2.5,cex=1,outer=F)

# Plot Figure 3__________________________________________________
XX<-trimtoRS2(TSct,cutoff=5000,cutoff2=3000,trim.type="to.RS",start=1000)
oTS<-XX$trimTS
tail<-XX$tail.length
adTS<-am5$adTS
regTS<-sampleTS(adTS,sample.by="distribute.samples",samp.freq1=NULL,nsamp,timeCT)
acTS<-sampleTSatAC(adTS,AC.buffer=0.10,AC.samp=0.4,timeCT,XX$tail.length)

variant.cols<-c("gray20","gray35","gray60","gray75")

dev.new(width=10,height=5.1)
par(mar=c(0,3.5,0.45,1),oma=c(4,3,3,2),mgp=c(2.4,1,0))
plot.matx<-c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,5,5,6,6,7,7,8,8,9,9,9,9,10,10,10,10,9,9,9,9,10,10,10,10)
nf<-layout(matrix(plot.matx, nrow = 8, ncol = 6, byrow = FALSE))

# Plot examples of paleo transformation
plot(oTS[,1],oTS[,2],pch=16,cex=0.7,col=variant.cols[1],ylab="",las=1,yaxt="n",xlim=c(500,9500),xaxt="n",ylim=c(0,1.2))
axis(2,at=seq(0.2,1,0.4),las=1,cex.axis=1.5)
text(700,1.25,"a)",cex=1.75,pos=1)

Xs<-XX$tail.length+regTS[,1]
Ys<-rep(1.16,length(Xs))
points(Xs,Ys,pch="|",cex=0.8)
text(max(Xs)-75,Ys[1],"E",pos=4,cex=1)

Xs<-XX$tail.length+acTS[,1]
Ys<-rep(1.02,length(Xs))
points(Xs,Ys,pch="|",col="black",cex=0.8)
text(max(Xs)-75,Ys[1],"T",pos=4,cex=1)

plot(adTS[,1]+oTS[1,1],adTS[,2],pch=16,cex=0.7,col=variant.cols[2],ylab="",las=1,yaxt="n",xlim=c(500,9500),xaxt="n",ylim=c(0,1.2))
axis(2,at=seq(0.2,1,0.4),las=1,cex.axis=1.5)
text(700,1.25,"b)",cex=1.75,pos=1)
plot(regTS[,1]+oTS[1,1],regTS[,2],pch=16,cex=0.7,col=variant.cols[3],ylab="",las=1,yaxt="n",xlim=c(500,9500),xaxt="n",ylim=c(0,1.2))
axis(2,at=seq(0.2,1,0.4),las=1,cex.axis=1.5)
text(700,1.25,"c)",cex=1.75,pos=1)
plot(acTS[,1]+oTS[1,1],acTS[,2],pch=16,cex=0.7,col=variant.cols[4],ylab="",las=1,yaxt="n",xlim=c(500,9500),xaxt="n",ylim=c(0,1.2))
axis(2,at=seq(0.2,1,0.4),las=1,cex.axis=1.5)
text(700,1.25,"d)",cex=1.75,pos=1)
axis(1,tcl=-0.25,padj=-0.5,cex.axis=1.5)

mtext("Time Steps (Years)",1,line=2.5,cex=1.2,outer=F)
mtext("Tree Cover (proportion)",2,line=-0.25,cex=1.2,out=T)

#plot detrended time series
oTS<-detrendTS(TSct,method="gaussian")
plot(oTS[,1],oTS[,2],col=variant.cols[1],xaxt="n",pch=16,cex=0.7,las=1,ylab="",yaxt="n",ylim=c(-0.25,0.25))
axis(2,seq(-0.2,0.2,0.2),las=1,tcl=-0.5,hadj=0.75,cex.axis=1.5)
text(min(oTS[,1])*100,0.19,"e)",cex=1.75)

adTS<-detrendTS(adTS,method="gaussian")
plot(adTS[,1]+TSct[1,1],adTS[,2],col=variant.cols[2],xaxt="n",pch=16,cex=0.7,las=1,ylab="",yaxt="n",ylim=c(-0.25,0.25))
axis(2,seq(-0.2,0.2,0.2),las=1,tcl=-0.5,hadj=0.75,cex.axis=1.5)
text(min(oTS[,1])*100,0.19,"f)",cex=1.75)
mtext("Detrended Proportional Tree Cover",2,line=2.75,cex=1.2,adj=0.75)

regTS<-detrendTS(regTS,method="gaussian")
plot(regTS[,1]+TSct[1,1],regTS[,2],col=variant.cols[3],xaxt="n",pch=16,cex=0.7,las=1,ylab="",yaxt="n",ylim=c(-0.25,0.25))
axis(2,seq(-0.2,0.2,0.2),las=1,tcl=-0.5,hadj=0.75,cex.axis=1.5)
text(min(oTS[,1])*100,0.19,"g)",cex=1.75)

acTS<-detrendTS(acTS,method="gaussian")
plot(acTS[,1]+TSct[1,1],acTS[,2],col=variant.cols[4],xaxt="n",pch=16,cex=0.7,las=1,ylab="",yaxt="n",ylim=c(-0.25,0.25))
axis(2,seq(-0.2,0.2,0.2),las=1,tcl=-0.5,hadj=0.75,cex.axis=1.5)
axis(1,tcl=-0.25,padj=-0.5,cex.axis=1.5)
text(min(oTS[,1])*100,0.19,"h)",cex=1.75)
mtext("Time Steps (Years)",1,line=2.5,cex=1.2,outer=F)

#plot SD EWS
Ktau<-matrix(NA,nrow=5,ncol=2)
Ktau[1,2]<-"tau"
Xsd<-std.sd(oTS,2500,1)
temp.matx<-cbind(Xsd$std.SD,Xsd$windows2[,1])
temp.matx2<-temp.matx[which(temp.matx[,2]<timeCT),]
Ktau[2,2]<-round(cor(temp.matx2[,1],temp.matx2[,2],method="kendall"),digits=2)

plot(Xsd$midpoints,Xsd$std.SD,xaxt="n",type="l",las=1,ylab="",ylim=c(0,5),col=variant.cols[1],yaxt="n",xlim=c(1000,9000),lwd=3)
text(1200,4.8,"i)",cex=1.75)
axis(2,seq(0,5,1),las=1,tcl=-0.25,hadj=0.5,cex.axis=1.5)
mtext("Standardized SD",2,line=2,cex=1.2,outer=F)
abline(v=timeCT-TSct[1,1])

Xsd<-std.sd(adTS,150,1)
temp.matx<-cbind(Xsd$std.SD,Xsd$windows2[,1])
temp.matx2<-temp.matx[which(temp.matx[,2]<timeCT),]
Ktau[3,2]<-round(cor(temp.matx2[,1],temp.matx2[,2],method="kendall"),digits=2)
lines(Xsd$midpoints,Xsd$std.SD,col=variant.cols[2],lwd=3)

Xsd<-std.sd(regTS,100,1)
temp.matx<-cbind(Xsd$std.SD,Xsd$windows2[,1])
temp.matx2<-temp.matx[which(temp.matx[,2]<timeCT),]
Ktau[4,2]<-round(cor(temp.matx2[,1],temp.matx2[,2],method="kendall"),digits=2)
lines(Xsd$midpoints,Xsd$std.SD,col=variant.cols[3],lwd=3)

Xsd<-std.sd(acTS,100,1)
temp.matx<-cbind(Xsd$std.SD,Xsd$windows2[,1])
temp.matx2<-temp.matx[which(temp.matx[,2]<timeCT),]
Ktau[5,2]<-round(cor(temp.matx2[,1],temp.matx2[,2],method="kendall"),digits=2)
lines(Xsd$midpoints,Xsd$std.SD,col=variant.cols[4],lwd=3)

legend("topright",c(Ktau[,2]),col=c(rgb(0,0,0,0),variant.cols),pch=16,bty="n")

#plot AC EWS
Ktau<-matrix(NA,nrow=5,ncol=2)
Ktau[1,2]<-"tau"
Xac<-ACtime(oTS,2500,1)
temp.matx<-cbind(Xac$ACtime,Xac$windows2[,1])
temp.matx2<-temp.matx[which(temp.matx[,2]<timeCT),]
Ktau[2,2]<-round(cor(temp.matx2[,1],temp.matx2[,2],method="kendall"),digits=2)

plot(Xac$midpoints,Xac$ACtime,xaxt="n",type="l",las=1,ylab="",yaxt="n",ylim=c(0,40),col=variant.cols[1],xlim=c(1000,9000),lwd=3)
text(1200,39,"j)",cex=1.75)
abline(v=timeCT-TSct[1,1])
axis(1,tcl=-0.25,padj=-0.5,cex.axis=1.5)
axis(2,seq(0,40,5),las=1,tcl=-0.25,hadj=0.5,cex.axis=1.5)
mtext("Autocorrelation Time",2,line=2,cex=1.2,outer=F)

Xac<-ACtime(adTS,150,1)
temp.matx<-cbind(Xac$ACtime,Xac$windows2[,1])
temp.matx2<-temp.matx[which(temp.matx[,2]<timeCT),]
Ktau[3,2]<-round(cor(temp.matx2[,1],temp.matx2[,2],method="kendall"),digits=2)
lines(Xac$midpoints,Xac$ACtime,col=variant.cols[2],lwd=3)

Xac<-ACtime(regTS,100,1)
temp.matx<-cbind(Xac$ACtime,Xac$windows2[,1])
temp.matx2<-temp.matx[which(temp.matx[,2]<timeCT),]
Ktau[4,2]<-round(cor(temp.matx2[,1],temp.matx2[,2],method="kendall"),digits=2)
lines(Xac$midpoints,Xac$ACtime,col=variant.cols[3],lwd=3)

Xac<-ACtime(acTS,100,1)
temp.matx<-cbind(Xac$ACtime,Xac$windows2[,1])
temp.matx2<-temp.matx[which(temp.matx[,2]<timeCT),]
Ktau[5,2]<-round(cor(temp.matx2[,1],temp.matx2[,2],method="kendall"),digits=2)
lines(Xac$midpoints,Xac$ACtime,col=variant.cols[4],lwd=3)
mtext("Time Steps (Years)",1,line=2.5,cex=1.2,outer=F)

legend("topright",c(Ktau[,2]),col=c(rgb(0,0,0,0),variant.cols),pch=16,bty="n")

# SIMULATIONS USING THE EXPONENTIAL AGE MODEL#######
# set model paramters_____________________________________________
 h = 0.5
 r=0.25
# c=1
delta_t = 1 
gens = 300
K_Start = 1
K_Pulse_amt = -0.4 
pulse_time = 3 
sigma_sd = 0.005
V0 = 1 
FRI = 1
beta_ps<-estBetaParams(mu=0.15, var=0.015)

# set taphonomic parameters____________________________________
agemodel="Carpenter"
a0=2
a1=0.025
nsamp=25
title<-"Exponential"
windows<-c(75,20,10,10)
exRStime=210
nsamp=200
AC.buff=0.1
samp.freq2=0.4
steps<-c(1,1,1,1)
a0=2
a1=0.025
cutoff=150
cutoff2=50
start=60

 

# Run simulations____________________________________
# occasionally this x1 simulation fails. Just re-run
system.time(x1<-rep.ews(TStype="TSct",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0,a1,cutoff,cutoff2,trim.type="to.RS",start,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(x2<-rep.ews(TStype="TSdc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff,cutoff2,trim.type="to.set.bounds",start,q=1,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(x3<-rep.ews(TStype="TSrs",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff,cutoff2,trim.type="to.RS",start,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(x4<-rep.ews(TStype="TSnc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff,cutoff2,trim.type="to.set.bounds",start,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

Xct<-x1
Xdc<-x2
Xrs<-x3
Xnc<-x4

# plot Figure 7_____________________________________________________
dev.new(width=7,height=5.5)
par(oma=c(4,4,4,2),mar=c(0.75,0.5,0,0.5))
nf<-layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,ncol=3,byrow=F))

colorCT<-rgb(0.1,0.3,0.4,1)
colorDC<-rgb(0.1,0.3,0.4,0.7)
colorRS<-rgb(0.1,0.3,0.4,0.5)
colorNC<-rgb(0.1,0.3,0.4,0.3)
mains2<-c("","","","")
ind<-"sd"
main<-""

plot.taph.hists(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=T,mains=mains2,ymax=1.05,labs2=NULL,letters=c("a)","b)","c)","d)"),title="Untransformed",type.label=F,taph.ind=1)
plot.taph.hists(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=c("e)","f)","g)","h)"),title="Sedimentation",type.label=F,taph.ind=2)
plot.taph.hists(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=c("i)","j)","k)","l)"),title="Sed.+Even Sampling",type.label=T,taph.ind=3)

mtext("Frequency",2,line=2.25,outer=T,cex=1.1)
mtext("Kendall's tau",1,line=2.5,outer=T,cex=1.1) 

# Exponential Supplemental___________________________________
plot.supp(Xct,Xdc,Xrs,Xnc,"sd","Exponential, Standard Deviation")
plot.supp(Xct,Xdc,Xrs,Xnc,"ac","Exponential, Autocorrelation Time")

# Exponential Supplemental table
Supp1<-Kt.summary.stats(Xct,Xdc,Xrs,Xnc,3,"sd")
Supp2<-Kt.summary.stats(Xct,Xdc,Xrs,Xnc,3,"ac")
	
supp1<-Supp1[,c(1,2,3,5,6,8,9,10)]
supp2<-Supp2[,c(1,2,3,5,6,8,9,10)]


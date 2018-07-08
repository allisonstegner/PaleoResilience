# Load packages_______________________________

library(breakpoint)
library(moments)

source('~/Grass-Wood model_final 30June2018.R')

# FUNCTIONS#####################################################################
# std.sd___________________________________________________________________________
# calculates standard deviation in a rolling window, standardized by the length of time represented by each window

std.sd<-function(tstemp,winlen,step.size){
	
	# INPUTS:
	# tstemp = a time series, as a 2 col matrix where [,1] is time and [,2] is an ecological variable
	# winlen = length of rolling window
	# step.size = increment by which to move the rolling window
	
	# RETURN: 
	# midpoints = vector of rolling window midpoints
	# std.SD = vector of standard deviation values for each rolling window
	# windows2 = a 2 colum matrix; start and end points for each window
	
	samp.points<-seq(1,length(tstemp[,2]),step.size)
	temp1<-cbind(samp.points,(samp.points+winlen)-1)

	temp2<-temp1[which(temp1[,2]<(nrow(tstemp))+1),]
	std.SD<-c()
	midpoints<-c()
	windows2<-matrix(NA,nrow=nrow(temp2),ncol=2)
	for (j in 1:nrow(temp2)){
		rowj<-temp2[j,]
		windows2[j,]<-c(tstemp[rowj[1],1],tstemp[rowj[2],1])
		midpoints.ind<-((rowj[[2]]-rowj[[1]])/2)+rowj[[1]]
		midpoints[j]<-tstemp[midpoints.ind,1]
		wait.time<-tstemp[rowj[[2]],1]-tstemp[rowj[[1]],1]
		win.pol.vals<-tstemp[c(rowj[[1]]:rowj[[2]]),2]				
		std.SD[j]<-sd(win.pol.vals*sqrt(abs(wait.time)))	
	}
	out<-list(midpoints=midpoints, std.SD=std.SD,windows2=windows2)
	return(out)
}


# ACtime___________________________________________________________________	
# calculates autocorrelation time with missing values in a rolling window.

ACtime<-function(tstemp,winlen,step.size){
	
	# INPUTS:
	# tstemp = a time series, as a 2 col matrix where [,1] is time and [,2] is an ecological variable
	# winlen = length of rolling window
	# step.size = increment by which to move the rolling window
	
	# RETURN: 
	# midpoints = vector of rolling window midpoints
	# ACtime = vector of autocrrelation time values for each rolling window
	# windows2 = a 2 colum matrix; start and end points for each window
	
	samp.points<-seq(1,length(tstemp[,2]),step.size)
	temp1<-cbind(samp.points,(samp.points+winlen)-1)
	temp2<-temp1[which(temp1[,2]<(nrow(tstemp))+1),]
	
	#fill gaps with NAs
	pol_vals<-tstemp[,2]
	steps<-1:(max(floor(tstemp[,1])))
	xx<-match(steps,floor(tstemp[,1]),nomatch=9999)
	xxx<-c()
	for (k in 1:length(xx)) {
		if (xx[k]==9999) { xxx[k]<-NA
		} else { xxx[k]<-pol_vals[xx[k]] }
	}

	ACt<-c()
	midpoints<-c()
	windows2<-matrix(NA,nrow=nrow(temp2),ncol=2)
	for (j in 1:nrow(temp2)){
		rowj<-temp2[j,]
		windows2[j,]<-c(tstemp[rowj[1],1],tstemp[rowj[2],1])
		midpoints.ind<-((rowj[[2]]-rowj[[1]])/2)+rowj[[1]]
		midpoints[j]<-tstemp[midpoints.ind,1]
		wait.time<-tstemp[rowj[[2]],1]-tstemp[rowj[[1]],1]
		ac.test<-acf(tstemp[c(rowj[[2]]:rowj[[1]]),2],lag=1,plot=FALSE,na.action=na.pass)
		ACt[j] = -1/log(abs(ac.test$acf[2]))
	}

	out<-list(midpoints=midpoints, ACtime=ACt,windows2=windows2)
	return(out)
}

# ews.summary______________________________________________________________________
# summarizes standard deviation and autocorrelation time rolling window values, and kendall's tau for each

ews.summary<-function(tstemp,winlen,step.size,timeCT){
	
	# INPUTS:
	# tstemp = 2 col matrix; [,1] is age, [,2] is ecological variable 
	# winlen = length of rolling window
	# step.size = increment by which to move the rolling window
	# timeCT = age (in time steps, not as an index) of the regime shift. usually passed from another function
 
	# RETURN:
	# sd.kendall = Kendall's tau values for correlation between standard deviation and time
	# sd.vals = rolling window standard deviation values
	# sd.k.p = p value for standard deviation kendall's tau
	# ac.kendall = Kendall's tau values for correlation between autocorrelation time and time
	# ac.vals = rolling window autocorrelation time values
	# ac.k.p = p value for autocorrelation time kendall's tau

	#calculate SD and Kendall's tau
	tstemp.sd<-std.sd(tstemp,winlen,step.size)
	temp3<-cbind(tstemp.sd$windows2[,2],tstemp.sd$std.SD)
	temp4<-temp3[which(temp3[,1]<timeCT),]
	sd.ken.test<-cor.test(temp4[,1],temp4[,2],method="kendall")
	sd.ken.p<-sd.ken.test$p.value
	sd.ken<-sd.ken.test$estimate[[1]]
	sd.vals<-temp3
	
	#calculateACtime and Kendall's tau
	actemp<-ACtime(tstemp,winlen,step.size)
	temp3<-cbind(actemp$windows2[,2],actemp$ACtime)
	temp4<-temp3[which(temp3[,1]<timeCT),]
	ac.ken.test<-cor.test(temp4[,1],temp4[,2],method="kendall")
	ac.ken.p<-ac.ken.test$p.value
	ac.ken<-ac.ken.test$estimate[[1]]
	ac.vals<-temp3
	
	out<-list(sd.kendall=sd.ken,sd.vals=sd.vals,sd.k.p=sd.ken.p,ac.kendall=ac.ken,ac.vals=ac.vals,ac.k.p=ac.ken.p)
	return(out)
}

# detrendTS__________________________________________________________________________
# detrend time series using linear, gaussian, lowess, or no detrending. Initial function form inspried by the earlywarnings::generic_ews function

detrendTS<-function(oTS,method,bandwidth=NULL,span=NULL, degree=NULL){
	# INPUTS:
	# oTS = generally, a 2 col matrix, with [,1] as time and [,2] as ecological variable. An untransformed single run of the GW model
	# method = should detrending be "linear" (using lm), "gaussian" (using ksmooth function), "lowess" (using loess function), or "none"?
	# bandwidth = bandwidth to use for gaussian detrending. Default is NULL, in which case bandwidth is optimized using bw.nrd0 (following the earlywarnings package)
	# span = span to use for loess smoothing. default is 25/100 (following the earlywarnings package)
	# degree = degree to use for loess smoothign. default is 2 (following the earlywarnings package)
	
	# RETURN:
	# a 2 col matrix with [,1] as time, and [,2] as detrended ecological variable
	
	TS<-oTS[,2]
	timeindex<-oTS[,1]
	if (method=="lowess") {
		if (is.null(span)) {  
			span <- 25/100 
		} else { 
			span <- span/100  
		}
        if (is.null(degree)) {  
        	degree <- 2 
        } else { 
        	degree <- degree   
        }
        smYY <- loess(TS ~ timeindex, span = span, degree = degree, normalize = FALSE, family = "gaussian")
        smY <- predict(smYY, data.frame(x = timeindex), se = FALSE)
        smoothTS <- TS - smY
        out<-cbind(timeindex,smoothTS)
	} else if (method=="gaussian"){
		if (is.null(bandwidth)) { 
			bw <- round(bw.nrd0(timeindex)) 
		} else { 
			bw <- round(length(TS) * bandwidth/100) 
		}
		smYY <- ksmooth(timeindex, TS, kernel = "normal", bandwidth = bw, range.x = range(timeindex), x.points = timeindex)
		smoothTS <- TS - smYY$y
		smY <- smYY$y
		 out<-cbind(timeindex,smoothTS)
    } else if (method == "linear") {
        smoothTS <- resid(lm(TS ~ timeindex))
        smY <- fitted(lm(TS ~ timeindex))
        out<-cbind(timeindex,smoothTS)
	} else if (method=="none") {
		smY<-TS
		smoothTS<-TS
		out<-cbind(timeindex,smoothTS)
	}
	return(out)
}

# multi.ews.summary.emp________________________________________________________________
# summarize rolling window metrics for empirical data

multi.ews.summary.emp<-function(oTS,TAbins,sample.at,sample.at2,timeCT,detrend_method,yrs.preCT){
	
	# INPUTS:
	# oTS = a 2 col matrix, with [,1] as time and [,2] as ecological variable
	# TAbins = vector of time averaging values in yrs/cm
	# sample.at = time points (in years) to sample oTS
	# sample.at2 = second set of time points (in years) to sample oTS
	# timeCT = time (in years) of regime shift 
	# detrend_method = detrending method: must be "linear" (using lm), "gaussian" (using ksmooth function), "lowess" (using loess function), or "none" (passed to detrendTS)
	# yrs.preCT = number of years to inlcude prior to the regime shift (e.g., if yrs.preCT=100, trim to the 100 years preceeding the regime shift)
	
	# RETURNS:
	# sdKs = standard deviaiton kendall's taus
	# acKs = autocorrelation time kendall's taus
	# sd.vals.multi = standard deviation rolling window values
	# ac.vals.multi = autocorrelation time rolling window values
	# sdKp = standard deviaiton kendall's p values
	# acKp = autocorrelation time kendall's p values
	# oTS = trimmed, detrended original time series
	# adTS = trimmed, detrended time-averaged time series
	# regTS = trimmed, detrended time-averaged and sampled time series
	# regTS2 = trimmed, detrended time-averaged and sampled time series
	
	OTS<-oTS[round(timeCT-yrs.preCT):(timeCT),]
	
	cTAbins<-c(1,cumsum((TAbins)))
	cTAbins<-floor(cTAbins[which(cTAbins<max(OTS[,1]))])
	cTAbins<-floor(cTAbins[which(cTAbins<yrs.preCT)])
	TAts<-c()
	for (j in 1:(length(cTAbins)-1)){
		TAts[j]<-mean(OTS[cTAbins[j]:cTAbins[(j+1)],2])
	}
	adTS<-cbind((cTAbins[1:length(TAts)]),TAts)
	#adTS<-adTS[intersect(which(adTS[,1]<max(OTS[,1])),which(adTS[,1]>min(OTS[,1]))),]
		
	sample.at1<-round(sample.at[sample.at<nrow(adTS)])
	regTS<-adTS[sample.at1,]
	
	sample.at2_2<-round(sample.at2[sample.at2<nrow(adTS)])
	regTS2<-adTS[sample.at2_2,]
	
	#detrend_method="gaussian"
	oTS<-detrendTS(OTS,method=detrend_method)
	adTS<-detrendTS(adTS,method=detrend_method)
	regTS<-detrendTS(regTS,method=detrend_method)
	regTS2<-detrendTS(regTS2,method=detrend_method)
	
	win1<-floor(nrow(oTS)/2)
	win2<-floor(nrow(adTS)/2)
	win3<-floor(nrow(regTS)/2)
	win4<-win3
		
	ews.orig<-ews.summary(oTS,win1,1,timeCT)
	ews.ad<-ews.summary(adTS,win2,1,timeCT)
	ews.reg<-ews.summary(regTS,win3,1,timeCT)
	ews.reg2<-ews.summary(regTS2,win4,1,timeCT)
	
	sdKs<-c(ews.orig$sd.kendall,ews.ad$sd.kendall,ews.reg$sd.kendall,ews.reg2$sd.kendall)
	sdKp<-c(ews.orig$sd.k.p,ews.ad$sd.k.p,ews.reg$sd.k.p,ews.reg2$sd.k.p)
	sd.vals.multi<-list(ews.orig$sd.vals,ews.ad$sd.vals,ews.reg$sd.vals,ews.reg2$sd.vals)
	names(sdKs)<-c("orig","ad","reg","reg2")
	names(sd.vals.multi)<-c("orig","ad","reg","reg2")
	names(sdKp)<-c("orig","ad","reg","reg2")
		
	acKs<-c(ews.orig$ac.kendall,ews.ad$ac.kendall,ews.reg$ac.kendall,ews.reg2$ac.kendall)
	acKp<-c(ews.orig$ac.k.p,ews.ad$ac.k.p,ews.reg$ac.k.p,ews.reg2$ac.k.p)
	ac.vals.multi<-list(ews.orig$ac.vals,ews.ad$ac.vals,ews.reg$ac.vals,ews.reg2$ac.vals)
	names(acKs)<-c("orig","ad","reg","reg2")
	names(ac.vals.multi)<-c("orig","ad","reg","reg2")
	names(acKp)<-c("orig","ad","reg","reg2")
	
	out<-list(sdKs=sdKs,acKs=acKs,sd.vals.multi=sd.vals.multi,ac.vals.multi=ac.vals.multi,sdKp=sdKp,acKp=acKp,oTS=OTS,adTS=adTS,regTS=regTS,regTS2=regTS2)
	return(out)
}

# rep.ews.emp________________________________________________________________________	
# summarize resilience indicators for multiple simulations, time averaged and sampled like an empirical dataset

rep.ews.emp<-function(TStype,nreps,TAbins,sample.at,sample.at2,detrend_method,yrs.preCT,pulse_time){
	
	# INPUTS:
	# TStype = corresponding to "TSct" (gradually-forced critical transition), "TSrs" (abruptly-forced critical transition), "TSdc" (gradually forced non-critical transition), or "TSnc" (no change)
	# nreps = number of simulations
	# TAbins = vector of time averaging in yrs/cm
	# sample.at = vector of sample locations, in years before present
	# sample.at2 = second vector of sample locations, in years before present
	# detrend_method = detrending method: must be "linear" (using lm), "gaussian" (using ksmooth function), "lowess" (using loess function), or "none" (passed to multi.ews.summary, "detrend_method")
	# yrs.preCT = amount of time in the empirical record prior to the regime shift
	# pulse_time = time (in time steps) to initiate driver pulse in the simulations
	
	#RETURNS:
	# sd.kendalls = a 4 x nreps column matrix, standard deviation kendall's tau values 
	# ac.kendalls = a 4 x nreps column matrix, autocorrelation time kendall's tau values
	# rep.sd.vals = a list, standard deviation rolling window values
	# rep.ac.vals = a list, autocorrelation time rolling window values
	# timeCTs = time of regime shift for each simulation
	# rep.sd.ps = a 4 x nreps column matrix, standard deviation kendall's p values
	# rep.ac.ps = a 4 x nreps column matrix, autocorrelation time kendall's p values	
		
	rep.sdKs<-matrix(NA,nrow=nreps,ncol=4)
	rep.sd.vals<-list()
	rep.sd.ps<-matrix(NA,nrow=nreps,ncol=4)
	
	rep.acKs<-matrix(NA,nrow=nreps,ncol=4)
	rep.ac.vals<-list()
	rep.ac.ps<-matrix(NA,nrow=nreps,ncol=4)

	timeCTs<-c()
	
	if ((TStype %in% c("TSct","TSrs","TSdc","TSnc"))==FALSE) {
		print("TStype must be 'TSct', 'TSrs', 'TSnc', or 'TSdc'")
	}
	
	if (TStype=="TSrs"){
		driver_topo<-"abrupt"
	} else {
		driver_topo<-"gradual"
	}
	
	for (i in 1:nreps){
			print(i)
			
			single_test = single_run(r=r, gens=gens, delta_t=delta_t, K_Start=K_Start, K_Pulse_amt=K_Pulse_amt, V0=V0, pulse_time=pulse_time,driver_press_topo=driver_topo,q=q)
			
			TS<-single_test[,3]
			origTS<-cbind(c(1:length(TS)),TS)
			
		if (TStype=="TSct" || TStype=="TSrs"){
			bp.out<-CE.Normal.Mean(as.data.frame(TS),Nmax=1)
			timeCT<-bp.out$BP.Loc
		} else {
			timeCT<-exRStime
		}
			
		timeCTs[i]<-timeCT
		
		summary.temp<-multi.ews.summary.emp(origTS,TAbins,sample.at,sample.at2,timeCT,detrend_method="gaussian",yrs.preCT)
							
		rep.sdKs[i,]<-summary.temp$sdKs
		rep.sd.vals[[i]]<-summary.temp$sd.vals.multi
		rep.sd.ps[i,]<-summary.temp$sdKp
			
		rep.acKs[i,]<-summary.temp$acKs
		rep.ac.vals[[i]]<-summary.temp$ac.vals.multi
		rep.ac.ps[i,]<-summary.temp$acKp
	}
	
	out<-list(sd.kendalls=rep.sdKs,ac.kendalls=rep.acKs,rep.sd.vals=rep.sd.vals,rep.ac.vals=rep.ac.vals,timeCTs=timeCTs,rep.sd.ps=rep.sd.ps,rep.ac.ps=rep.ac.ps)
	
	return(out)
}


# plot.taph.hists.sub____________________________________________________________
# plot overlapping histograms and ROC data

plot.taph.hists.sub.overlay<-function(Xct,Xdc,Xrs,Xnc,indicator,yaxis,mains,ymax,labs2,letters,title,type.label,taph.ind,empiricalK,draw.xaxt,plot.empK){
	# INPUTS:
	# INPUTS:
	# Xct, Xdc, Xrs, Xnc = simulation outputs from rep.ews function. Not necessarily, but generally: Xct = return for gradually-forced critical transitions, Xdc = = return for gradually-forced non-critical transitions, Xrs = return for abruptly-forced critical transitions, Xnc = return for no change
	# indicator = "sd" to choose standard deviation, or "ac" to choose autocorrelation time
	# yaxis = T/F, should the y axis be plotted? 
	# mains = list of headings
	# ymax = maximum y axis value
	# labs2 = labels for each tile
	# letters = letters for each tile
	# title = title for entire set of plots
	# type.label = T/F, should the time series type name be plotted on the right-hand axis?
	# taph.ind = interger between 1 and 4. 1 = untransformed, 2 = age-depth transformed, 3 = age-depth and sampling transformed, 4 = age-depth and targeted sampling transformed
	# empiricalK = supply a Kendall's tau value for the empirical dataset
	# draw.xaxt = T/F, should the x axis be drawn?
	# plot.empK = T/F, should the empirical Kendall's tau be drawn?
	
	bins<-seq(-1,1,0.1)
	prob.table<-list()
	for (i in taph.ind){
		ts.list<-c("orig","ad","reg","ac")
		if (indicator=="sd"){
			kendalls<-"sd.kendalls"
		} else if (indicator=="ar") {
			kendalls<-"ar.kendalls"
		}else if (indicator=="sk"){
			kendalls<-"sk.kendalls"
		} else if (indicator=="ac"){
			kendalls<-"ac.kendalls"
		} else { print("indicator must be 'sd', 'ar', 'sk', or 'ac'")	}

		nreps<-length(Xnc[[kendalls]][,i])
		
		A<-hist(Xct[[kendalls]][,i],breaks=bins,plot=F)
		Aprops<-A$counts/nreps
		
		C<-hist(Xrs[[kendalls]][,i],plot=FALSE,breaks=bins)
		Cprops<-C$counts/nreps
		counts<-cbind(A$counts,C$counts)
			
		barplot(Aprops,ylim=c(0,ymax),main="",col=colorCT,cex.axis=1.5,cex.main=1,las=1,yaxt="n")
				
		if (plot.empK==T){
			RS.interp<-approx(c(-1,1),c(0,24),xout=empiricalK)
		abline(v=RS.interp$y,col="black",lwd=1.5)
		}
		
		IN<-Xct[[kendalls]][,i]
		OUT<-Xrs[[kendalls]][,i]
		if (min(IN)>max(OUT)) {
			text(-1,ymax*0.6,"No overlap",pos=4,cex=0.8)
		} else {
			roc<-calcROC(IN, OUT, max.len = 10000)
			if (roc$AUC<0.5){ 
				AUC<-1-round(roc$AUC,2)
			} else { AUC<-round(roc$AUC,2) }
			opt.interp<-approx(c(-1,1),c(0,24),xout=roc$optimal)
			#auc<-paste(AUC,"+-",round(roc$se.fit,2))
			auc<-paste(AUC)
			abline(v=opt.interp$y,lty=3,lwd=2,col="black")
			if (round(roc$optimal,3)<0){ xloc=-1
			} else { xloc=-1 }	
		}
		
		text(xloc,ymax*0.65,paste("AUC =", auc,sep=" "),pos=4,cex=1.1)
		text(xloc,ymax*0.45,paste("opt. =", round(roc$optimal,2),sep=" "),pos=4,cex=1.1)
						
		barplot(Cprops,ylim=c(0,ymax),main="",col=colorRS,add=T,yaxt="n")
	
		if (draw.xaxt==T){
			temp.table<-round(t(apply(counts,1,function(x) x/sum(x))),digits=4)	
			mids<-barplot(t(temp.table),plot=FALSE)
			mstep<-(mids[2]-mids[1])/2
			locs<-mids-(mids[2]-mids[1])/2
			step<-locs[2]-locs[1]
			at.vect<-c(locs,locs[length(locs)]+step)
			locs3<-at.vect[seq(1,length(at.vect),2)]
			axis(1,at=locs3,labels=seq(-1,1,0.2),tick=T,las=2,cex=0.6,line=0.25,mgp=c(3,0.5,0.5),tcl=-0.25,cex.axis=1.5)	
		}
		
	}
}

# calcROC_______________________________________________________________
# function for calculating ROC for two distributions of Kendall's tau values, modified from analogue::roc 

calcROC <- function(IN, OUT, max.len = 10000) {
		# INPUTS:
		# IN = values for the the first distribution
		# OUT = values for the the second distribution
		# max.len = numeric; length of analolgue and non-analogue vectors. Used as limit to thin points on ROC curve to (from analogue::roc)
		
		# RETURNS:
		# TPF = the true positive fraction
		# FPE = the false positive error
		# optimal  = optimal dissimilarity value where slope of the ROC curve is maximal
		# AUC = area under the ROC curve
		# se.fit = standard error of the AUC estimate
		# n.in = numeric, number of samples within the current group
		# n.out = numeric, number of samples not in the current group
		# p.value = p value of a wilcoxon rank sum test on the two sets of dissimilarities. Also known as a Mann-Whitney test
		# roc.points =  the unique dissimilarities at which the roc curve was evaluated
		# max.roc = position along the ROC curve at which the slope of the ROC curve is maximal. this is the index of this point on the curve
		# prior = numeric, length 2. Vector of observed probabilities of true analogue and true non-analogues in the group
		# analogue = a list with componenents yes and no containing the dissimilarities for true analgoe and true non-analogues in the group
		
    	n.IN <- length(IN)
    	n.OUT <- length(OUT)
    	g <- rep(c(TRUE, FALSE), times = c(n.IN, n.OUT))
    	tab <- table(c(IN, OUT), g)
    	TPF <- cumsum(tab[, 2])/sum(tab[, 2])
    	FPE <- cumsum(tab[, 1])/sum(tab[, 1])
    	roc.values <- TPF - FPE
    	roc.points <- rev(sort(as.numeric(dimnames(tab)[[1]])))
    	optimal <- as.numeric(names(max.roc <- which.max(abs(roc.values))))
    	names(FPE) <- names(TPF) <- names(roc.values) <- NULL
    	wilcox <- wilcox.test(IN, OUT, conf.int = FALSE)
    	AUC <- 1 - (wilcox$statistic/(n.IN * n.OUT))
    	 AUC2 <- AUC^2
    	q1 <- AUC/(2 - AUC)
    	q2 <- (2 * AUC2)/(1 + AUC)
    	se.fit <- AUC * (1 - AUC) + ((n.IN - 1) * (q1 - AUC2)) + ((n.OUT - 1) * (q2 - AUC2))
    	se.fit <- sqrt(se.fit/(n.IN * n.OUT))
    	p.value <- wilcox$p.value
    	prior <- c(n.IN, n.OUT)/sum(n.IN, n.OUT)
    	retval <- list(TPF = TPF, FPE = FPE, optimal = optimal, AUC = AUC, se.fit = se.fit, n.in = n.IN, n.out = n.OUT, p.value = p.value, roc.points = roc.points, max.roc = max.roc, prior = prior, analogue = list(yes = IN, no = OUT))
    	retval
    }

# Kt.summary.stats_________________________________________
# generate summary stats for simulation Kendall's tau values

Kt.summary.stats<-function(Xct,Xdc,Xrs,Xnc,digits,indicator){
	t1<-multi.ks.test(Xct,Xdc,Xrs,Xnc,ts="orig",indicator,p.correct=F,alpha=0.05)
	t2<-multi.ks.test(Xct,Xdc,Xrs,Xnc,ts="ad",indicator,p.correct=F,alpha=0.05)
	t3<-multi.ks.test(Xct,Xdc,Xrs,Xnc,ts="reg",indicator,p.correct=F,alpha=0.05)
	t4<-multi.ks.test(Xct,Xdc,Xrs,Xnc,ts="ac",indicator,p.correct=F,alpha=0.05)
	col1<-c(rep("Untransformed",4),rep("AD",4),rep("AD+SAMP",4),rep("AD+TSAMP",4))
	col2<-rep(c("CT","LD","RS","NC"),4)
	
	temp.stats<-rbind(t1$ind.stats,t2$ind.stats,t3$ind.stats,t4$ind.stats)
	temp.stats<-round(temp.stats,digits)
	out<-cbind(col1,col2,temp.stats)
	return(out)
}

# multi.ks.test_________________________________________
multi.ks.test<-function(Xct,Xdc,Xrs,Xnc,ts,indicator, p.correct,alpha){
		if (ts=="orig"){ ts.ind<-1 } 
		else if (ts=="ad"){ ts.ind<-2 } 
		else if (ts=="reg"){ ts.ind<-3 }
		else if (ts=="ac"){ ts.ind<-4 }
		else { print("ts must be 'orig','ad','reg',or'ac'") }
		
		if (indicator=="sd"){
			kens<-"sd.kendalls"
			ken.ps<-"rep.sd.ps"
		} else if (indicator=="sk"){
			kens<-"sk.kendalls"
			ken.ps<-"rep.sk.ps"
		} else if (indicator=="ar"){
			kens<-"ar.kendalls"
			ken.ps<-"rep.ar.ps"
		} else if (indicator=="ac"){
			kens<-"ac.kendalls"
			ken.ps<-"rep.ac.ps"
		} else { print("indicator must be 'sd', 'sk', 'ar', or 'ac'") }
		
		CTvDC<-ks.test(Xct[[kens]][,ts.ind],Xdc[[kens]][,ts.ind])
		CTvRS<-ks.test(Xct[[kens]][,ts.ind],Xrs[[kens]][,ts.ind])
		CTvNC<-ks.test(Xct[[kens]][,ts.ind],Xnc[[kens]][,ts.ind])
		DCvRS<-ks.test(Xdc[[kens]][,ts.ind],Xrs[[kens]][,ts.ind])
		DCvNC<-ks.test(Xdc[[kens]][,ts.ind],Xnc[[kens]][,ts.ind])
		RSvNC<-ks.test(Xrs[[kens]][,ts.ind],Xnc[[kens]][,ts.ind])

		ps<-c(CTvDC$p.value,CTvRS$p.value,CTvNC$p.value,DCvRS$p.value,DCvNC$p.value,RSvNC$p.value)
				
		if (p.correct==TRUE){
			ps.out<-p.adjust(ps, method="bonferroni")
			ps.out<-t(ps.out)
		} else { ps.out<-ps }
		
		ind.stats<-matrix(NA,nrow=4,ncol=9)
		ind.stats[1,]<-c(summary(Xct[[kens]][,ts.ind]),quantile(Xct[[kens]][,ts.ind],0.025),quantile(Xct[[kens]][,ts.ind],0.975),sum(Xct[[ken.ps]][,ts.ind]<alpha)/length(Xct[[ken.ps]][,ts.ind])*100)
		ind.stats[2,]<-c(summary(Xdc[[kens]][,ts.ind]),quantile(Xdc[[kens]][,ts.ind],0.025),quantile(Xdc[[kens]][,ts.ind],0.975),sum(Xdc[[ken.ps]][,ts.ind]<alpha)/length(Xdc[[ken.ps]][,ts.ind])*100)
		ind.stats[3,]<-c(summary(Xrs[[kens]][,ts.ind]),quantile(Xrs[[kens]][,ts.ind],0.025),quantile(Xrs[[kens]][,ts.ind],0.975),sum(Xrs[[ken.ps]][,ts.ind]<alpha)/length(Xrs[[ken.ps]][,ts.ind])*100)
		ind.stats[4,]<-c(summary(Xnc[[kens]][,ts.ind]),quantile(Xnc[[kens]][,ts.ind],0.025),quantile(Xnc[[kens]][,ts.ind],0.975),sum(Xnc[[ken.ps]][,ts.ind]<alpha)/length(Xnc[[ken.ps]][,ts.ind])*100)
		colnames(ind.stats)<-c("Min.","1st Qu","Median","Mean","3rd Qu.","Max.","2.5% CI","97.5% CI","% Significant")
		rownames(ind.stats)<-c("TSct","TSdc","TSrs","TSnc")
		
		indp<-matrix(NA,nrow=3,ncol=3)
		indp[,1]<-ps.out[c(1:3)]
		indp[,2]<-c(NA,ps.out[c(4,5)])
		indp[,3]<-c(NA,NA,ps.out[6])
		colnames(indp)<-c("CT","DC","RS")
		rownames(indp)<-c("DC","RS","NC")
		
		out<-list(ps.out=ps.out,ind.stats=ind.stats,p.table=indp)
		return(out)
	}
	
	
##########
# import Woody Cover data
steel<-read.csv("~/Steel Lake.csv")
tree<-steel[,4]/100
chron<-steel[,1]

#import quartz+feldspar data (only necessary to draw figure)
SLqf<-read.csv("~SLqf.csv")

# Calculate depostion rates
depth<-steel[,6]
depth<-depth-min(depth)+1
depth.difs<-depth[-1]-depth[-length(depth)]
age.difs<-chron[-1]-chron[-length(chron)]
dep.rates<-age.difs/depth.difs
dep.tab<-cbind(depth[-length(depth)],dep.rates)

# create a data table for LWO
data.tab<-cbind(rev(chron[-length(chron)]),rev(depth[-length(depth)]),rev(dep.rates),rev(tree[-length(tree)]))
data.tab<-data.tab[data.tab[,3]>0,]
colnames(data.tab)<-c("ages","depth","dep.rates","AP")

# find time of regime shift
bp.out<-CE.Normal.Mean(as.data.frame(tree),Nmax=4)
timeCT<-chron[bp.out$BP.Loc][3]
start<-chron[bp.out$BP.Loc][4]

trunc.data<-data.tab[(which(data.tab[,1]==start)+1):which(data.tab[,1]==timeCT),]
plot(chron,tree,col="black",pch=16)
points(trunc.data[,c(1,4)],col="red",pch=16)
abline(v=c(start,timeCT))

# find length of time prior to regime shift at LWO
yrs.preCT<-round(start-timeCT)
yrs.preCT2<-round(9400-timeCT)

# run rolling sd for full high woody cover data
win<-floor(length(trunc.data[,1])/2)
det.pol<-detrendTS(trunc.data[,c(1,4)],method="gaussian")
tree.sd<-std.sd(det.pol,winlen=win,step.size=1)
time.vect<-max(tree.sd$windows2[,2])-tree.sd$windows2[,2]
empiricalK<-cor(tree.sd$std.SD,time.vect,method="kendall")

# run rolling sd for full period of aridification woody cover data
trunc.data2<-data.tab[min(which(data.tab[,1]<9400)):which(data.tab[,1]==timeCT),]
det.pol2<-detrendTS(trunc.data2[,c(1,4)],method="gaussian")
win2<-floor(length(det.pol2[,1])/2)
tree.sd2<-std.sd(det.pol2,winlen=win2,step.size=1)
time.vect2<-max(tree.sd2$windows2[,2])-tree.sd2$windows2[,2]
empiricalK2<-cor(tree.sd2$std.SD,time.vect2,method="kendall")

#plot sd and pollen data
par(mfrow=c(3,1),mar=c(1,1,0,1),oma=c(3,3,1,1))
plot(data.tab[,c(1,4)],pch=1,xlim=c(0,12500),xaxt="n")
lines(data.tab[,c(1,4)])
points(data.tab[,c(1,4)],pch=16,cex=0.8,col="white")
points(trunc.data[,c(1,4)],pch=16)
points(trunc.data2[,c(1,4)],pch=1,col="blue")
abline(v=timeCT)
mtext("%WC",2,line=2)

plot(det.pol,xlim=c(0,12500),pch=16,xaxt="n")
lines(det.pol)
points(det.pol2,col="blue",pch=16)
lines(det.pol2,col="blue")
abline(v=timeCT)
mtext("detrended %WC",2,line=2)

plot(tree.sd$windows[,2],tree.sd$std.SD,xlim=c(0,12500),cex=1,pch=16,ylim=c(0,2.5))
points(tree.sd2$windows[,2],tree.sd2$std.SD,col="blue",pch=1,cex=1.2)
abline(v=timeCT)

legend("topright",legend=c(round(empiricalK,digits=2),round(empiricalK2,digits=2)),text.col=c("black","blue"))
mtext("SD",2,line=2)

# Simulations_______________________________________________________________________
# set number of gens, and amt. of time before the rgime shift based on observed at LWO
gens<-ceiling(max(data.tab[,1]))

# set up vector of time averaging
interpTA<-approx(data.tab[,2],data.tab[,3],xout=seq(min(data.tab[,2]),max(data.tab[,2]),1))
TAbins<-interpTA$y
cTAbins<-c(1,cumsum(TAbins))
cTAbins<-cTAbins[which(cTAbins<gens)]

# Set model parameters__________________________
h = 0.5 
r=0.25 
q=5
c=1  
delta_t = 1 
K_Start = 1
K_Pulse_amt = -0.3 
pulse_time = 2000 
sigma_sd = 0.005
phi = 0.05
V0 = 1 
FRI = 1
beta_ps<-estBetaParams(mu=0.15, var=0.015) 

# Set taphonomic parameters__________________________
exRStime=timeCT
nreps=100
steps<-c(1,1,1,1)

# run simulations for each time series type__________________________
# set up vector of time averaging
interpTA<-approx(trunc.data[,2],trunc.data[,3],xout=seq(min(trunc.data[,2]),max(trunc.data[,2]),1))
TAbins<-interpTA$y
cTAbins<-c(1,cumsum(TAbins))
cTAbins<-cTAbins[which(cTAbins<yrs.preCT)]

# set up sampling
sample.at<-max(trunc.data[,2])-trunc.data[,2]+1 #observed sampling at LWO
sample.at2<-round(seq(min(sample.at),max(sample.at),length.out=length(sample.at)*2))

Xct<-rep.ews.emp(TStype="TSct",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT,pulse_time=2000)
Xdc<-rep.ews.emp(TStype="TSdc",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT,pulse_time=2000)
Xrs<-rep.ews.emp(TStype="TSrs",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT,pulse_time = yrs.preCT+100)
Xnc<-rep.ews.emp(TStype="TSnc",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT,pulse_time=2000)

# run simulations for each time series type, for the shorter Steel Lake time series
# set up vector of time averaging
interpTA<-approx(trunc.data2[,2],trunc.data2[,3],xout=seq(min(trunc.data2[,2]),max(trunc.data2[,2]),1))
TAbins<-interpTA$y
cTAbins<-c(1,cumsum(TAbins))
cTAbins<-cTAbins[which(cTAbins<yrs.preCT2)]

# set up sampling
sample.at<-max(trunc.data2[,2])-trunc.data2[,2]+1 #observed sampling at LWO
sample.at2<-round(seq(min(sample.at),max(sample.at),length.out=length(sample.at)*2))

Xct2<-rep.ews.emp(TStype="TSct",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT2,pulse_time=2000)
Xdc2<-rep.ews.emp(TStype="TSdc",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT2,pulse_time=2000)
Xrs2<-rep.ews.emp(TStype="TSrs",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT2,pulse_time = yrs.preCT2+100)
Xnc2<-rep.ews.emp(TStype="TSnc",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT2,pulse_time=2000)


# plot results________________________________________
tree<-tree*100
ind="sd"
colorCT<-rgb(1,0.5,0.5)
colorRS<-rgb(0,0,1,0.5)
mains2<-c("","","","")

dev.new(width=6.5,height=7)

matx<-c(
1,1,1,2,2,2,3,3,4,4,5,6,6,
1,1,1,2,2,2,3,3,4,4,5,6,6,
1,1,1,2,2,2,3,3,4,4,5,6,6,
1,1,1,2,2,2,3,3,4,4,5,6,6,
1,1,1,2,2,2,3,3,4,4,5,6,6,
1,1,1,2,2,2,3,3,4,4,5,6,6,
1,1,1,2,2,2,3,3,7,7,5,8,8,
1,1,1,2,2,2,3,3,7,7,5,8,8,
1,1,1,2,2,2,3,3,7,7,5,8,8,
1,1,1,2,2,2,3,3,7,7,5,8,8,
1,1,1,2,2,2,3,3,7,7,5,8,8,
1,1,1,2,2,2,3,3,7,7,5,8,8,
1,1,1,2,2,2,3,3,9,9,9,9,9)
nf<-layout(matrix(matx,nrow=13,ncol=13,byrow=F))
layout.show(nf)
par(oma=c(4,5,1,5),mar=c(1,1,0.5,0),xpd=F)

empiricalK<-round(empiricalK,digits=2) # Kendall's tau for full high woody cover period
empiricalK2<-round(empiricalK2,digits=2) # Kendall's tau period of aridification

plot(SLqf,yaxt="n",xaxt="n",xlim=c(12000,0),type="l",ylim=c(0,80),col="blue",lwd="2")
polygon(x=c(9400,9400,8000,8000),c(-20,100,100,-20),col=rgb(1,0.8,0,0.3),border=NA)
lines(SLqf,col="blue",lwd="2")
axis(4,at=seq(0,80,20),las=1,col="blue",col.ticks="blue",col.axis="blue",cex.axis=1.5)
lab<-"% Q+F"
mtext(lab,4,line=3.25,col="blue")

par(new=T)

preRS<-data.tab[data.tab[,1]>timeCT,]
plot(data.tab[,1],data.tab[,4]*100,pch=1,xlim=c(12000,0),ylim=c(0,80),xaxt="n",las=1,cex.axis=1.5,yaxt="n")
lines(data.tab[,1],data.tab[,4]*100)
points(data.tab[,1],data.tab[,4]*100,pch=16,col="white",cex=0.8)
points(data.tab[,1],data.tab[,4]*100,pch=1)
points(trunc.data[,1],trunc.data[,4]*100,pch=16)
points(trunc.data2[,1],trunc.data2[,4]*100,pch=16,col="gray")
points(trunc.data2[,1],trunc.data2[,4]*100,pch=1)
axis(2,at=seq(0,100,20),cex.axis=1.5,las=1)
mtext("% WC",2,line=3.25,cex=1.1)
mtext("a)",3,line=-1.6,cex=1,adj=0.015)
abline(v=timeCT,lty=3)

preRSsd<-tree.sd$std.SD[tree.sd$windows2[,2]>timeCT]
plot(tree.sd$windows2[tree.sd$windows2[,2]>timeCT,2],preRSsd,xlim=c(12000,0),pch=16,las=1,yaxt="n",ylim=c(0,2.25),cex.axis=1.5)
axis(2,at=seq(0,2,1),cex.axis=1.5,las=1)
lines(tree.sd$windows2[tree.sd$windows2[,2]>timeCT,2],preRSsd)
preRSsd2<-tree.sd2$std.SD[tree.sd2$windows2[,2]>timeCT]
lines(tree.sd2$windows2[tree.sd2$windows2[,2]>timeCT,2],preRSsd2)
points(tree.sd2$windows2[tree.sd2$windows2[,2]>timeCT,2],preRSsd2,col="gray",pch=16)
points(tree.sd2$windows2[tree.sd2$windows2[,2]>timeCT,2],preRSsd2,pch=1)
mtext("SD",2,line=3.25,cex=1.1)
mtext("Calendar YBP",1,line=2.75,cex=1.1)
mtext("b)",3,line=-1.6,cex=1,adj=0.015)
abline(v=timeCT,lty=3)

legend("topright",c("11,100-8,200 YBP, Kt = -0.27","9,400-8,200 YBP, Kt = 0.82"),pt.bg=c("black","gray"),col="black",pch=21,bty="n",cex=1.1,pt.cex=1.5)

plot(c(0,1),c(0,1),type="n",bty="n",xaxt="n",yaxt="n",ylab="",xlab="") #plot spacer

ind="sd"
colorCT<-rgb(1,0.5,0.5)
colorRS<-rgb(0,0,1,0.5)
mains2<-c("","","","")

plot.taph.hists.sub.overlay(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=T,mains=mains2,ymax=0.4,labs2=NULL,letters=letters.list,title="",type.label=F,taph.ind=3,empiricalK,draw.xaxt=F,plot.empK=T)
axis(2,at=seq(0,0.6,0.2),las=1,cex.axis=1.5)
mtext("Frequency",2,line=3,outer=F,cex=1)
mtext("Observed Sampling",3,cex=1.2,line=2)
mtext("c) 11,100-8,200 YBP",3,line=0,adj=0.05)

plot(c(0,1),c(0,1),type="n",bty="n",xaxt="n",yaxt="n",ylab="",xlab="") #plot spacer

plot.taph.hists.sub.overlay(Xct2,Xdc2,Xrs2,Xnc2,indicator=ind,yaxis=T,mains=mains2,ymax=0.4,labs2=NULL,letters=letters.list,title="",type.label=F,taph.ind=3,empiricalK2,draw.xaxt=T,plot.empK=T)
axis(2,las=1,at=seq(0,0.6,0.2),cex.axis=1.5)
mtext("Frequency",2,line=3,outer=F,cex=1)
mtext("d) 9,400-8,200 YBP",3,line=0,adj=0.05)
mtext("Kendall's tau",1,line=3.5,cex=1.1)

plot.taph.hists.sub.overlay(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=0.4,labs2=NULL,letters=letters.list,title="",type.label=F,taph.ind=4,empiricalK,draw.xaxt=F,plot.empK=F)
axis(2,at=seq(0,0.6,0.2),labels=F)
mtext("e) 11,100-8,200 YBP",3,line=0,adj=0.05)
mtext("2x Observed Sampling",3,cex=1.2,line=2)

letters.list<-c("i) Gradual CT","","j) Abrupt CT","")
plot.taph.hists.sub.overlay(Xct2,Xdc2,Xrs2,Xnc2,indicator=ind,yaxis=F,mains=mains2,ymax=0.4,labs2=NULL,letters=letters.list,title="",type.label=F,taph.ind=4,empiricalK2,draw.xaxt=T,plot.empK=F)
axis(2,at=seq(0,0.6,0.2),labels=F)
mtext("f) 9,400-8,200 YBP",3,line=0,adj=0.05)
mtext("Kendall's tau",1,line=3.5,cex=1.1)

par(xpd=NA)
plot(c(0,1),c(0,1),type="n",bty="n",xaxt="n",yaxt="n",ylab="",xlab="")

legend(-0.8,0.9,c("Grad. CT","Abrupt. CT"),fill=c(colorCT,colorRS),bty="n",cex=1.1,pt.cex=1.5)
legend(-0.8,1.1,c("SL Kt","optimum"),lty=c(1,3),bty="n",cex=1.1,pt.cex=1.5,lwd=2)
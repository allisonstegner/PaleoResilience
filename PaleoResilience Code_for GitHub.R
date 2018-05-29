#######################################
# Paleoecological Resilience Indicator functions
# Stegner et al. 
# updated 29 May 2018
#######################################

# PACKAGES  #####################################################################
library(breakpoint)
library(analogue)
library(neotoma)
require(stats)
require(moments)

# used to generate the Lake West Okoboji age model
#library(Bchron) 

####FUNCTIONS###################################################################
# estBetaParam_____________________________________________________________
# given a certain mean fire mortaility and standard deviation, this function estmates alpha and beta parameters

estBetaParams <- function(mu, var) {
	# INPUTS:
	# mu = (approximate) mean of the beta distribution
	# var = standard devaition of the beta distribution
	
	# RETURN: 
	# alpha, beta parameters for a beta distribution

	alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
	beta <- alpha * (1 / mu - 1)
	return(params = list(alpha = alpha, beta = beta))
}

# growth_func_______________________________________________
#this function is outputs the change in vegetation biomass over a time, dt. 

growth_func = function(r, K, X, q, h, delta_t, fire_sev) {
	
	# INPUTS:
	# r = maximum tree growth rate
	# K = carrying capacity (driver)
	# X = tree cover
	# q = hill slope parameter for fire mortality
	# h = half saturation constant
	# delta_t = dt time step used for numerical simulations
	# fire_sev (c) = maximum proportional tree mortality
	# note that "X" is "T" and "fire_sev" is the parameter "c", both in reference to equation 2 of the main text. 
	
	# RETURN: vector of tree cover values
	dXdt = delta_t *(r*X*(K-X))-fire_sev *X*(h^q / (X^q+h^q)) 
	result  = dXdt + X
	return(result)
}

# red.noise_______________________________________________
# this function generates a vector of red (autocorrealted) noise

red.noise<-function(N,phi,sigma_sd)  {
	# INPUTS:
	# N = length of time series
	# phi = lag-1 autocorrelation
	# sigma = standard deviation of noise
	
	# RETURN: a vector of red (autocorrealted) noise
	
	white = rnorm(N,0,1)
	x.1 = rep(0,N)
	x.1[1] = white[1]
	for(i in 2:N) { x.1[i] = phi*x.1[i-1] + sqrt( (1-phi^2) )*white[i] }
	x.sigma = sigma_sd*x.1
	return(x.sigma)
}

#single_run________________________________________________________________________
# this function performs a single simulation that varies the K parameter over time

single_run = function (r, gens, delta_t, K_Start, K_Pulse_amt, V0, pulse_time, driver_press_topo, q){
	
	# INPUTS:
	# r = maximum tree growth rate
	# gens = number of time steps
	# delta_t = dt time step used for numerical simulations
	# K_Start = initial value for the K parameter
	# K_Pulse_amt = amount that K changes during the driver pulse applied to the system
	# V0 = initial tree cover
	# pulse_time = time step at which to begin driver pulse
	# driver_press_topo = whether to apply the driver abruptly (default) or gradually ("gradual")
	# q = hill slope parameter for fire mortality
	
	# RETURN: a three column matrix
	# [,1] = time step (t)
	# [,2] = K at time t
	# [,3] = tree cover at time t
	
	K_Pulse_time=gens-pulse_time
	#this matrix holds the results of the simulation
	sim_l <- ceiling(gens/delta_t)
	results = matrix(1:sim_l, nrow = sim_l, ncol = 3)
	colnames(results) = c("Time", "K_param", "Tree_cover[proportional]")
	results[1,3] = V0

	# generating a vector of values of "K" for the time series 
	K_mat = matrix(K_Start, nrow = sim_l, ncol =1)
	K_mat[(pulse_time/delta_t):(pulse_time/delta_t+K_Pulse_time/delta_t)] = K_Start+K_Pulse_amt 
	if (driver_press_topo == "gradual") K_mat[(pulse_time/delta_t):(pulse_time/delta_t+K_Pulse_time/delta_t)] =  seq(from = K_Start, to =K_Pulse_amt+K_Start, length.out = K_Pulse_time/delta_t+1)
	results[,2] = K_mat

	# generating a vector of red noise
	dt.noise = sqrt(delta_t)
	Rnoise = red.noise(sim_l,phi,sigma_sd)
	
	# generating a vector of values of "c" based on a beta-distribution
	fire_mort_vect<-rbeta(n=sim_l, shape1=beta_ps$alpha,shape2=beta_ps$beta)
	treecover<-rep(V0, sim_l)
	
	# running the simulation
	for (i in 2:(sim_l)) {
		temp <- growth_func(r=r, K=K_mat[i], X=treecover[i-1], q=q, h=h, delta_t=delta_t, fire_sev=fire_mort_vect[i]) + Rnoise[i]
		treecover[i] <- max(0, temp)
		results[i,3] <- treecover[i]
	}
	
	# this is a sub-setting procedure for simulations run with a time-step different than 1 
 	spl <- seq(1,sim_l,by = 1/delta_t)

	final_results = results[spl,]
	final_results[,1] = final_results[,1]*delta_t

	return(final_results)
}


# trimtoRS2___________________________________________________________________________
# time a time series to the portion prior to regime shift, using breakpoint analysis
# returns a two column vector, with data from prior to the regime shift

trimtoRS2<-function(TS,cutoff,cutoff2,trim.type,start){
	# INPUTS:
	# TS = 2 col matrix, with [,1] as time and [,2] as ecological variable
	# cutoff = amount of time to retain prior to regime shift 
	# cutoff2 = amount of time to retain after regime shift
	# trim.type = should the time of regime shift be determined analytically ("to.RS"), arbitarily ("to.set.bounds"), or should no trimming take place ("none")?
	# start = if trim.type = "to.set.bounds", start is the first time step to include
	
	# RETURN:
	# trimTS = trimmed time series 2 col matrix, with [,1] as time and [,2] as ecological variable
	# tail.length = amount of time trimmed from the beginning of the time series, used for lining up trimmed and untrimmed versions
	
	if (trim.type=="to.RS"){
		bp.out<-CE.Normal.Mean(as.data.frame(TS[,2]),Nmax=1)
		timeCT<-TS[bp.out$BP.Loc,1]
		trimTSt<-TS[which(TS[,1]<timeCT),]
		trimTS<-TS[c((timeCT-cutoff):(timeCT+cutoff2)),]
		tail.length<-nrow(trimTSt)-cutoff
	} else if (trim.type=="to.set.bounds"){
		trimTS<-TS[c(start:(start+cutoff+cutoff2)),]
		tail.length<-start
	} else if (trim.type=="none"){
		trimTS<-TS
		tail.length<-nrow(TS)
	} else {
		print("trim.type must be 'to.RS','set.bounds', or 'none'")
	}
	out<-list(trimTS=trimTS,tail.length=tail.length)
	return(out)
}

# brokenstick.timeavgTS_______________________________________________________
# This function applies time averaging: two constant rates of time averaging, with a single breakpoint (expressed as a time step) at which time averaging changes from one rate to the other. Also applies a small amount of noise around the # of years per cm, and calculates of centimeters in the core

brokenstick.timeavgTS<-function(oTS,TAbottom,TAtop,breakpoint,sd.pct){
	# INPUTS:
	# oTS = generally, a 2 col matrix, with [,1] as time and [,2] as ecological variable. An untransformed single run of the GW model
	# TAbottom = time averaging in number of years per cm at the bottom (older) of the core
	# TAtop = time averaging in number of years per cm at the top (younger) of the core
	# breakpoint = time step (as a irow number index) of transition from TAbottom to TAtop 
	
	# RETURN:
	# adTS = a two colum matrix, with time as [,1] (oldest age of the time aberaged bin) and time-averaging ecological variable as [,2]
	# TAbins = vector of time averaging durations generated by rnorm sampling of TAvect
	# TAvect = vector of mean time averaging used to generate TAbins
	
	TS<-oTS[,2]
	
	# divide the time series in half, according to the breakpoint 
	TSold<-TS[1:breakpoint]
	TSyoung<-TS[(breakpoint+1):length(TS)]
	
	# create a vector of mean time averaging values that is the same length as the time series
	TAvec<-c(rep(TAbottom,length(TSold)),rep(TAtop,length(TSyoung)))
	
	# add noise to the time averaging vector
	rTA<-c()
	for (i in 1:length(TAvec)){
		rTA[i]<-round(rnorm(1,TAvec[i],sd.pct*TAvec[i]))
	}
	rTA[which(rTA==0)]<-1

	
	TAbins<-rep(0,length(TS))
	TAbins[1]<-rTA[1]
	for (i in 2:length(TS)) {
		TAbins[i]<-rTA[sum(TAbins)+1]
		if (sum(TAbins,na.rm=T)>length(TS)) break
	}
	
	TAbins<-TAbins[which(TAbins>0)]
	cTAbins<-c(1,cumsum(TAbins))
	cTAbins<-cTAbins[which(cTAbins<length(TS))]

	TAts<-c()
	for (j in 1:(length(cTAbins)-1)){
		TAts[j]<-mean(TS[cTAbins[j]:cTAbins[(j+1)]])
	}
	
	adTS<-cbind((cTAbins[1:length(TAts)]),TAts)
	out<-list(adTS=adTS,TAbins=TAbins,TAvect=rev(TAvec))
	return(out)
}

# timeavgTS___________________________________________________________________________
# applies time averaging as a slow linear change from one rate to another, adds a small amount of noise around the # of years per cm, and calculates # of centimeters in the core

timeavgTS<-function(oTS,TAbottom,TAtop,sd.pct){
	# INPUTS:
	# oTS = generally, a 2 col matrix, with [,1] as time and [,2] as ecological variable. An untransformed single run of the GW model
	# TAbottom = time averaging in number of years per cm at the bottom (older) of the core
	# TAtop = time averaging in number of years per cm at the top (younger) of the core
	
	# RETURN:
	# adTS = a two colum matrix, with time as [,1] (oldest age of the time aberaged bin) and time-averaging ecological variable as [,2]
	# TAbins = vector of time averaging durations generated by rnorm sampling of TAvect
	# TAvect = vector of mean time averaging used to generate TAbins
	
	TS<-oTS[,2]

	m<-(TAbottom-TAtop)/(length(TS)-1)
	b<-TAtop-m

	TA<-c()
	rTA<-c()
	for (i in 1:length(TS)){
		TA[i]<-m*i+b
		rTA[i]<-round(rnorm(1,TA[i],sd.pct*TA[i]))
	}
	rTA[which(rTA<=0)]<-1

	TAbins<-c()
	TAbins[1]<-rTA[1]
	for (i in 2:length(TS)) {
		TAbins[i]<-rTA[sum(TAbins)+1]
		if (sum(TAbins,na.rm=T)>length(TS)) break
	}
	
	TAbins<-TAbins[which(TAbins>0)]
	cTAbins<-c(1,cumsum(TAbins))
	cTAbins<-cTAbins[which(cTAbins<length(TS))]

	TAts<-c()
	for (j in 1:(length(cTAbins)-1)){
		TAts[j]<-mean(TS[cTAbins[j]:cTAbins[(j+1)]])
	}
	
	adTS<-cbind((cTAbins[1:length(TAts)]),TAts)
	out<-list(adTS=adTS,TAbins=TAbins,TAvect=TA)
	return(out)
}

# Carpenter_timeavgTS____________________________________________________________________
# compresses a time series using exponential sedimentation, according to Taranu et al. (in review) method

Carpenter_timeavgTS<-function(oTS,a0,a1,tail){
	# INPUTS:
	# oTS = generally, a 2 col matrix, with [,1] as time and [,2] as ecological variable. An untransformed single run of the GW model
	# a0 = compression intercepts parameter
	# a1 = compression slope parameter
	# tail = tail.length from trimtoRS2; can be set to 0 if time series has not been trimmed
	
	# RETURN:
	# adTS = a two colum matrix, with time as [,1] (oldest age of the time aberaged bin) and time-averaging ecological variable as [,2]
	# TAbins = vector of time averaging durations generated by rnorm sampling of TAvect
	# TAvect = vector of mean time averaging used to generate TAbins
	
	TS<-oTS[,2]
	
	# reverse the time series to start at top of core
	Wmix.rev = rev(TS) 
	gens<-length(Wmix.rev)
		
	# Compute core length in cm
	L.core = (1/a1)*( log(a1*gens + a0) - log(a0))
	
	# Discard non-integer fraction and end of core
	L.core = floor(L.core)
	
	# Average tracer concentrations for each 1 cm slice of core
	core = rep(0,L.core)
	t0 = 0 # start at the top of the core
	t1 = 0
	# save matrix of t0 and t1 values to check results
	t01mat = matrix(0,nr=L.core,nc=2)
	
	# Do the first slice separately to accommodate starting at exactly 0
	# Compute new t1 from forumula in 'Pseudocode_Mixing+Compression_2016-08-25.doc'
	t1 = (1/a1)*( (a1*t0 + a0)*exp(a1) - a0)
	# Save t0 and t1
	t01mat[1,]=c(t0,t1)
	# average the tracer concentrations between t0 and t1
	y.slice = ceiling(t1)-floor(t0)
	wts = rep(1,y.slice)
	x = Wmix.rev[1:ceiling(t1)]
	wts[1] = ceiling(t0)-t0 # fraction of the first year in the slice
	wts[y.slice] = t1 - floor(t1) # fraction of the final year in the slice
	core[1] = x%*%wts/sum(wts)
	
	# Do the remaining slices
	for(i in 2:L.core) {  # loop over core slices
		t0 = t1 # previous t1 is now t0
		# Compute new t1 from forumula in 'Pseudocode_Mixing+Compression_2016-08-25.doc'
		t1 = (1/a1)*( (a1*t0 + a0)*exp(a1) - a0)
		# Save t0 and t1
		t01mat[i,]=c(t0,t1)
		# average the tracer concentrations between t0 and t1
		y.slice = ceiling(t1)-floor(t0)+1
		wts = rep(1,y.slice)
		x = Wmix.rev[floor(t0):ceiling(t1)]
		wts[1] = ceiling(t0)-t0 # fraction of the first year in the slice
		wts[y.slice] = t1 - floor(t1) # fraction of the final year in the slice
		core[i] = x%*%wts/sum(wts)
	}
		
		# Compute times at the center of each core slice
		t01mat[1,1]=1.e-3 # replace the zero to prevent problem with log
		T.core = exp( 0.5*(log(t01mat[,1])+log(t01mat[,2])) )
		# reverse the sequence of times to correspond with the simulated series
		T.core = (gens - floor(T.core))+tail
		
		adTS<-cbind(rev(T.core),rev(core))
		TAbins<-cbind(rev(t01mat[,1]),rev(t01mat[,2]))
		TAvect=TAbins[,2]-TAbins[,1]
		out<-list(adTS=adTS,TAbins=TAbins,TAvect=rev(TAvect))
		return(out)
}

# sampleTS___________________________________________________________________________
# this function subsamples a time series at regular intervals

sampleTS<-function(TS2,sample.by,samp.freq1,nsamp,timeCT){
	# INPUTS:
	# TS2 = a matrix with 2 cols, [,1] is ages, [,2] is ecological value, e.g. tree cover
	# sample.by = should sampling occur at predefined depth increments ("depth") across entire time sereis, or should a set number of samples be distrbuted at regular increments ("distribute.samples") prior to timeCT?
	# sample.freq1 = an integer, increment at which to sample TS2 if sample.by="depth"
	# nsamp = an integer, number of samples of TS2 to take is sample.by="distribute.samples"
	# timeCT = time of the ecological transition/regime shift. Often passed from another function. 
	
	# RETURN: 
	# sampTS = a time series. A 2 col matrix with [,1]=ages and [,2]=ecological value. Produced by subsampling TS2
	
	
	if (sample.by=="depth"){
		sTS<-TS2[seq(1,nrow(TS2),samp.freq1),]
	} else if (sample.by=="distribute.samples") {
		TStemp<-TS2
		sTS<-TStemp[round(seq(1,nrow(TStemp),length.out=nsamp)),]
	} else {
		print("sample.by must be 'depth' or 'distribute.samples'")
	}
	return(sampTS=sTS)
}

# sampleTSatAC___________________________________________________________________________
# this function takes samples from a time series at regular intervals, but adds extra samples immediately prior to a specified time point

sampleTSatAC<-function(TS2,AC.buffer,AC.samp,timeCT,start.at){		
	
	# INPUTS:
	# TS2 = a matrix with 2 cols, [,1] is ages, [,2] is ecological value, e.g. tree cover
	# sample.by = should sampling occur at predefined depth increments ("depth") across entire time sereis, or should a set number of samples be distrbuted at regular increments ("distribute.samples") prior to timeCT?
	# sample.freq1 = an integer, increment at which to sample TS2 if sample.by="depth"
	# sample.freq2: an integer, increment at which to sample prior to ecological transition
	# AC.buffer: between 0 and 1; proportion of the time series to sample intensively. E.g., if AC.buffer=0.25, function with intensively sample the 25% of the TS immediately prior to the ecological transition
	# nsamp = an integer, number of samples of TS2 to take is sample.by="distribute.samples"
	# timeCT = time of the ecological transition/regime shift. Often passed from another function. 
	
	# RETURN: 
	# acTS = a time series. A 2 col matrix with [,1]=ages and [,2]=ecological value. Produced by subsampling TS2
	
	ACbuff<-(max(TS2[,1])-min(TS2[,1]))*AC.buffer	
	rng<-c(((timeCT-ACbuff):timeCT-start.at))
	ACzone<-TS2[intersect(which(TS2[,1] > min(rng)), which(TS2[,1]<max(rng))),]
	
	pass1<-TS2[c(1:(which(TS2[,1]==min(ACzone[,1])))),]
	pass3<-TS2[c((which(TS2[,1]==max(ACzone[,1]))):nrow(TS2)),]
		
	nsamp1<-(nsamp-nsamp*AC.samp)
	nsamp2<-nsamp-nsamp1
	
	seg1<-pass1[seq(1,nrow(pass1),length.out=(nrow(pass1)/nrow(TS2))*nsamp1),]
	seg2<-ACzone[seq(1,nrow(ACzone),length.out=nsamp2),]
	seg3<-pass3[seq(1,nrow(pass3),length.out=(nrow(pass3)/nrow(TS2))*nsamp1),]
			
	AC.sampled<-rbind(seg1,seg2,seg3)
	return(acTS=AC.sampled)
}

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

# multi.ews.summary________________________________________________________________
multi.ews.summary<-function(origTS,TAbottom,TAtop,sample.by,nsamp,samp.freq1,samp.freq2,AC.buff,windows,steps,timeCT,TStype,agemodel,breakpoint,detrend_method,a0,a1,cutoff,cutoff2,trim.type,start,sd.pct,AC.samp){
	
	XX<-trimtoRS2(origTS,cutoff,cutoff2,trim.type,start)
	oTS<-XX$trimTS
	
	if (agemodel=="brokenstick"){
		tats1<-brokenstick.timeavgTS(XX$trimTS,TAbottom,TAtop,breakpoint,sd.pct)
	} else if (agemodel=="linearTA") {
		tats1<-timeavgTS(XX$trimTS,TAbottom,TAtop,sd.pct)
	} else if (agemodel=="Carpenter"){
		tats1<-Carpenter_timeavgTS(XX$trimTS,a0,a1,XX$tail.length)
	} else {
		print("agemodel must be 'brokenstick', 'linearTA', or 'Carpenter'")
	}
	
	ta.time<-tats1$adTS[,1]
	adTS<-cbind(ta.time,tats1$adTS[,2])
	regTS<-sampleTS(adTS,sample.by,samp.freq1=NULL,nsamp,timeCT)
	
	if (TStype=="TSnc" || TStype=="TSdc"){
		acTS<-sampleTSatAC(adTS,AC.buffer=AC.buff,AC.samp=AC.samp,exRStime,start.at=XX$tail.length)
	} else {
		acTS<-sampleTSatAC(adTS,AC.buffer=AC.buff,AC.samp=AC.samp,timeCT,start.at=XX$tail.length)
	}
	
	oTS<-detrendTS(oTS,method=detrend_method)
	adTS<-detrendTS(adTS,method=detrend_method)
	regTS<-detrendTS(regTS,method=detrend_method)
	acTS<-detrendTS(acTS,method=detrend_method)
	
	ews.orig<-ews.summary(oTS,windows[1],steps[1],timeCT)
	ews.ad<-ews.summary(adTS,windows[2],steps[2],timeCT)
	ews.reg<-ews.summary(regTS,windows[3],steps[3],timeCT)
	ews.ac<-ews.summary(acTS,windows[4],steps[4],timeCT)
	
	sdKs<-c(ews.orig$sd.kendall,ews.ad$sd.kendall,ews.reg$sd.kendall,ews.ac$sd.kendall)
	sdKp<-c(ews.orig$sd.k.p,ews.ad$sd.k.p,ews.reg$sd.k.p,ews.ac$sd.k.p)
	sd.vals.multi<-list(ews.orig$sd.vals,ews.ad$sd.vals,ews.reg$sd.vals,ews.ac$sd.vals)
	names(sdKs)<-c("orig","ad","reg","ac")
	names(sd.vals.multi)<-c("orig","ad","reg","ac")
	names(sdKp)<-c("orig","ad","reg","ac")

	acKs<-c(ews.orig$ac.kendall,ews.ad$ac.kendall,ews.reg$ac.kendall,ews.ac$ac.kendall)
	acKp<-c(ews.orig$ac.k.p,ews.ad$ac.k.p,ews.reg$ac.k.p,ews.ac$ac.k.p)
	ac.vals.multi<-list(ews.orig$ac.vals,ews.ad$ac.vals,ews.reg$ac.vals,ews.ac$ac.vals)
	names(acKs)<-c("orig","ad","reg","ac")
	names(ac.vals.multi)<-c("orig","ad","reg","ac")
	names(acKp)<-c("orig","ad","reg","ac")
	
	out<-list(sdKs=sdKs,acKs=acKs,sd.vals.multi=sd.vals.multi,ac.vals.multi=ac.vals.multi,sdKp=sdKp,acKp=acKp)
	return(out)
}

# multi.ews.summary.emp________________________________________________________________
multi.ews.summary.emp<-function(oTS,TAbins,sample.at,timeCT,detrend_method,yrs.preCT){
	OTS<-oTS[round(timeCT-yrs.preCT):(timeCT),]
	
	cTAbins<-c(1,cumsum(TAbins))
	cTAbins<-round(cTAbins[which(cTAbins<max(oTS[,1]))])
	TAts<-c()
	for (j in 1:(length(cTAbins)-1)){
		TAts[j]<-mean(oTS[cTAbins[j]:cTAbins[(j+1)],2])
	}
	adTS<-cbind((cTAbins[1:length(TAts)]),TAts)
	adTS<-adTS[intersect(which(adTS[,1]<max(OTS[,1])),which(adTS[,1]>min(OTS[,1]))),]
			
	sample.at<-sample.at[sample.at<nrow(adTS)]
	regTS<-adTS[sample.at,]
	
	sample.at2<-sample.at2[sample.at2<nrow(adTS)]
	regTS2<-adTS[sample.at2,]
	
	#detrend_method="gaussian"
	oTS<-detrendTS(OTS,method=detrend_method)
	adTS<-detrendTS(adTS,method=detrend_method)
	regTS<-detrendTS(regTS,method=detrend_method)
	regTS2<-detrendTS(regTS2,method=detrend_method)
	
	win1<-floor(nrow(oTS)/2)
	win2<-floor(nrow(adTS)/2)
	#win3<-floor(nrow(regTS)/2)
	win3<-11
	win4<-11
		
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
	
	out<-list(sdKs=sdKs,acKs=acKs,sd.vals.multi=sd.vals.multi,ac.vals.multi=ac.vals.multi,sdKp=sdKp,acKp=acKp,oTS=OTS,adTS=adTS,regTS=regTS,reTS2=regTS2)
	return(out)
}


# rep.ews________________________________________________________________________	
rep.ews<-function(TStype,nreps, TAbottom,TAtop,sample.by,nsamp,samp.freq1,samp.freq2,AC.buff,windows,steps,agemodel,breakpoint,det_method,a0,a1,cutoff,cutoff2,trim.type,start,q,sd.pct,AC.samp){
		
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
		pulse_time=exRStime
	} else {
		driver_topo<-"gradual"
		pulse_time=pulse_time
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
			
		summary.temp<-multi.ews.summary(origTS,TAbottom,TAtop,sample.by,nsamp,samp.freq1,samp.freq2,AC.buff,windows,steps,timeCT,TStype=TStype,agemodel,breakpoint,det_method,a0,a1,cutoff,cutoff2,trim.type,start,sd.pct,AC.samp)
			
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


# rep.ews.emp________________________________________________________________________	
rep.ews.emp<-function(TStype,nreps,TAbins,sample.at,sample.at2,detrend_method,yrs.preCT){
		
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
		pulse_time=exRStime
	} else {
		driver_topo<-"gradual"
		pulse_time=pulse_time
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
			
		summary.temp<-multi.ews.summary.emp(origTS,TAbins,sample.at,timeCT,detrend_method,yrs.preCT)
			
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




# calcROC_______________________________________________________________
# function for calculating ROC 
# modified from analogue::roc
# in order to calculate roc for 
# two distributions of Kendall's tau values

calcROC <- function(IN, OUT, max.len = 10000) {
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
    #AUC <- (wilcox$statistic/(n.IN * n.OUT))
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

# plot.taph.hists_______________________________________________________________________
plot.taph.hists<-function(Xct,Xdc,Xrs,Xnc,indicator,yaxis,mains,ymax,labs2,letters,title,type.label,taph.ind){
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
		B<-hist(Xdc[[kendalls]][,i],plot=FALSE,breaks=bins)
		Bprops<-B$counts/nreps
		C<-hist(Xrs[[kendalls]][,i],plot=FALSE,breaks=bins)
		Cprops<-C$counts/nreps
		D<-hist(Xnc[[kendalls]][,i],plot=FALSE,breaks=bins)
		Dprops<-D$counts/nreps
		counts<-cbind(A$counts,B$counts,C$counts,D$counts)
		
		type.labels<-c("Grad. CT","Grad. nonCT","Abrupt CT","No Change")
		prop.list<-list(Aprops,Bprops,Cprops,Dprops)
		colors<-c(colorCT,colorDC,colorRS,colorNC)
		OUTS<-list(c(Xrs[[kendalls]][,i],Xdc[[kendalls]][,i],Xnc[[kendalls]][,i]),Xdc[[kendalls]][,i],Xrs[[kendalls]][,i],Xnc[[kendalls]][,i])
		
		for (j in 1:4){
			
			barplot(prop.list[[j]],ylim=c(0,ymax),main="",col=colors[j],cex.axis=1.2,cex.main=1.5,las=1,yaxt="n")
						
			IN<-Xct[[kendalls]][,i]
			OUT<-OUTS[[j]]
			
			if (j==1) {
				mtext(title,3,line=0.5,cex=1)
				LTY=1
			} else { LTY=3 }
			
			
			if (min(IN)>max(OUT)) {
				text(-1,ymax*0.6,"No overlap",pos=4,cex=1)
			} else {
				roc<-calcROC(IN, OUT, max.len = 10000)
				if (roc$AUC<0.5){ 
					AUC<-1-round(roc$AUC,3)
				} else { AUC<-round(roc$AUC,3) }
				
				opt.interp<-approx(c(-1,1),c(0,24),xout=roc$optimal)
				auc<-paste(AUC,"+-",round(roc$se.fit,2))
				abline(v=opt.interp$y,lty=LTY,lwd=1.5,col="gray30")
				
				if (round(roc$optimal,3)<0){ xloc=-1
				} else { xloc=-1 }
				
				text(xloc,ymax*0.6,paste("AUC =", auc),pos=4,cex=1.1)
				text(xloc,ymax*0.45,paste("optimum =", round(roc$optimal,3)),pos=4,cex=1.1)
			}
			
			if (yaxis==T) { axis(2,las=1,hadj=0.75,tcl=-0.25,cex.axis=1.2) 
			} else { axis(2,las=1,tcl=-0.25,labels=F,cex.axis=1.2) }
			if (type.label==T){ mtext(type.labels[j],4,outer=F,cex=1) }
				
			text(-0.9,ymax*0.9,letters[j],pos=4,cex=1.5)
			#mtext(mains[1],3,line=-2,cex=0.7)
		}

		temp.table<-round(t(apply(counts,1,function(x) x/sum(x))),digits=4)	
		mids<-barplot(t(temp.table),plot=FALSE)
		mstep<-(mids[2]-mids[1])/2
		locs<-mids-(mids[2]-mids[1])/2
		step<-locs[2]-locs[1]
		at.vect<-c(locs,locs[length(locs)]+step)
		locs3<-at.vect[seq(1,length(at.vect),2)]
		axis(1,at=locs3,labels=seq(-1,1,0.2),tick=T,las=2,cex.axis=1.2,line=0.25,mgp=c(3,0.5,0.5),tcl=-0.25)	
	}
}

# plot.taph.hists.emp____________________________________________________________
plot.taph.hists.emp<-function(Xct,Xdc,Xrs,Xnc,indicator,yaxis,mains,ymax,labs2,letters,title,type.label,taph.ind,empiricalK){
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
		B<-hist(Xdc[[kendalls]][,i],plot=FALSE,breaks=bins)
		Bprops<-B$counts/nreps
		C<-hist(Xrs[[kendalls]][,i],plot=FALSE,breaks=bins)
		Cprops<-C$counts/nreps
		D<-hist(Xnc[[kendalls]][,i],plot=FALSE,breaks=bins)
		Dprops<-D$counts/nreps
		counts<-cbind(A$counts,B$counts,C$counts,D$counts)
		
		RS.interp<-approx(c(-1,1),c(0,24),xout=empiricalK)
				
		barplot(Aprops,ylim=c(0,ymax),main="",col=c(colorCT),cex.axis=0.8,cex.main=1,las=1,yaxt="n")
		mtext(title,3,line=0.5,cex=0.7)
		mtext(mains[1],3,line=-2,cex=0.7)
		text(-0.9,ymax*0.8,letters[1],pos=4,cex=1.2)
		
		abline(v=RS.interp$y,col="red")
		
		if (yaxis==T) { axis(2,las=1,hadj=0.75,tcl=-0.25) 
		} else { axis(2,las=1,tcl=-0.25,labels=F)}
		if (type.label==T){ mtext("Int. Regime Shift",4,outer=F,cex=0.6) }
				
		IN<-Xct[[kendalls]][,i]
		OUT<-c(Xrs[[kendalls]][,i],Xdc[[kendalls]][,i],Xnc[[kendalls]][,i])
		if (min(IN)>max(OUT)) {
			text(-1,ymax*0.6,"No overlap",pos=4,cex=0.8)
		} else {
			roc<-calcROC(IN, OUT, max.len = 10000)
			if (roc$AUC<0.5){ 
				AUC<-1-round(roc$AUC,3)
			} else { AUC<-round(roc$AUC,3) }
			opt.interp<-approx(c(-1,1),c(0,24),xout=roc$optimal)
			auc<-paste(AUC,"+-",round(roc$se.fit,2))
			abline(v=opt.interp$y,lty=1,lwd=1.5,col="gray40")
			if (round(roc$optimal,3)<0){ xloc=-1
			} else { xloc=-1 }
			text(xloc,ymax*0.6,paste("AUC =", auc),pos=4,cex=0.8)
			text(xloc,ymax*0.5,paste("optimum =", round(roc$optimal,3)),pos=4,cex=0.8)
		}

		
		barplot(Bprops,ylim=c(0,ymax),main="",col=c(colorDC),cex.axis=0.8,cex.main=1,las=1,yaxt="n")
		mtext(mains[2],3,line=-1,cex=0.7)	
		text(-0.9,ymax*0.8,letters[2],pos=4,cex=1.2)
		
		abline(v=RS.interp$y,col="red")		
		if (yaxis==T) { axis(2,las=1,hadj=0.75,tcl=-0.25) 
		} else { axis(2,las=1,tcl=-0.25,labels=F)}
		if (type.label==T){ mtext("Non-hysteretic",4,outer=F,cex=0.6) }
		
		OUT<-Xdc[[kendalls]][,i]
		if (min(IN)>max(OUT)) {
			text(-1,ymax*0.6,"No overlap",pos=4,cex=0.8)
		} else {
			roc<-calcROC(IN, OUT, max.len = 10000)
			if (roc$AUC<0.5){ 
				AUC<-1-round(roc$AUC,3)
			} else { AUC<-round(roc$AUC,3) }
			opt.interp<-approx(c(-1,1),c(0,24),xout=roc$optimal)
			auc<-paste(AUC,"+-",round(roc$se.fit,2))
			abline(v=opt.interp$y,lty=3,lwd=1.5,col="gray40")
			if (round(roc$optimal,3)<0){ xloc=-1
			} else { xloc=-1 }
			text(xloc,ymax*0.6,paste("AUC =", auc),pos=4,cex=0.8)
			text(xloc,ymax*0.5,paste("optimum =", round(roc$optimal,3)),pos=4,cex=0.8)
		}
	
		barplot(Cprops,ylim=c(0,ymax),main="",col=c(colorRS),cex.axis=0.8,cex.main=1,las=1,yaxt="n")
		mtext(mains[3],3,line=-1,cex=0.7)	
		text(-0.8,ymax*0.8,letters[3],pos=4,cex=1.2)	
		
		abline(v=RS.interp$y,col="red")
		if (yaxis==T) { axis(2,las=1,hadj=0.75,tcl=-0.25) 
		} else { axis(2,las=1,tcl=-0.25,labels=F) }
		if (type.label==T){ mtext("Ex. Regime Shift",4,outer=F,cex=0.6) }
		
		OUT<-Xrs[[kendalls]][,i]
		if (min(IN)>max(OUT)) {
			text(-1,ymax*0.6,"No overlap",pos=4,cex=0.8)
		} else {
			roc<-calcROC(IN, OUT, max.len = 10000)
			if (roc$AUC<0.5){ 
				AUC<-1-round(roc$AUC,3)
			} else { AUC<-round(roc$AUC,3) }
			opt.interp<-approx(c(-1,1),c(0,24),xout=roc$optimal)
			auc<-paste(AUC,"+-",round(roc$se.fit,2))
			abline(v=opt.interp$y,lty=3,lwd=1.5,col="gray40")
			if (round(roc$optimal,3)<0){ xloc=-1
			} else { xloc=-1 }
			text(xloc,ymax*0.6,paste("AUC =", auc),pos=4,cex=0.8)
			text(xloc,ymax*0.5,paste("optimum =", round(roc$optimal,3)),pos=4,cex=0.8)
		}
		
		barplot(Dprops,ylim=c(0,ymax),main="",col=c(colorNC),cex.axis=0.8,cex.main=1,las=1,yaxt="n")	
		mtext(mains[4],3,line=-1,cex=0.7)	
		text(-0.9,ymax*0.8,letters[4],pos=4,cex=1.2)	
		
		abline(v=RS.interp$y,col="red")	
		mtext("Kendall's Tau",1,outer=F,line=2.2,cex=0.7)					
		if (yaxis==T) { axis(2,las=1,hadj=0.75,tcl=-0.25)
		} else { axis(2,las=1,tcl=-0.25,labels=F) }
		if (type.label==T){ mtext("No Change",4,outer=F,cex=0.6) }
		
		OUT<-Xnc[[kendalls]][,i]
		if (min(IN)>max(OUT)) {
			text(-1,ymax*0.6,"No overlap",pos=4,cex=0.8)
		} else {
			roc<-calcROC(IN, OUT, max.len = 10000)
			if (roc$AUC<0.5){ 
				AUC<-1-round(roc$AUC,3)
			} else { AUC<-round(roc$AUC,3) }
			opt.interp<-approx(c(-1,1),c(0,24),xout=roc$optimal)
			auc<-paste(AUC,"+-",round(roc$se.fit,2))
			abline(v=opt.interp$y,lty=3,lwd=1.5,col="gray40")
			if (round(roc$optimal,3)<0){ xloc=-1
			} else { xloc=-1 }
			text(xloc,ymax*0.6,paste("AUC =", auc),pos=4,cex=0.8)
			text(xloc,ymax*0.5,paste("optimum =", round(roc$optimal,3)),pos=4,cex=0.8)
		}

		temp.table<-round(t(apply(counts,1,function(x) x/sum(x))),digits=4)	
		mids<-barplot(t(temp.table),plot=FALSE)
		mstep<-(mids[2]-mids[1])/2
		locs<-mids-(mids[2]-mids[1])/2
		step<-locs[2]-locs[1]
		at.vect<-c(locs,locs[length(locs)]+step)
		locs3<-at.vect[seq(1,length(at.vect),2)]
		axis(1,at=locs3,labels=seq(-1,1,0.2),tick=T,las=2,cex=0.6,line=0.25,mgp=c(3,0.5,0.5),tcl=-0.25)	
	}
}

# plot.taph.hists.sub____________________________________________________________
plot.taph.hists.sub<-function(Xct,Xdc,Xrs,Xnc,indicator,yaxis,mains,ymax,labs2,letters,title,type.label,taph.ind,empiricalK){
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
		B<-hist(Xdc[[kendalls]][,i],plot=FALSE,breaks=bins)
		Bprops<-B$counts/nreps
		C<-hist(Xrs[[kendalls]][,i],plot=FALSE,breaks=bins)
		Cprops<-C$counts/nreps
		D<-hist(Xnc[[kendalls]][,i],plot=FALSE,breaks=bins)
		Dprops<-D$counts/nreps
		counts<-cbind(A$counts,B$counts,C$counts,D$counts)
		
		RS.interp<-approx(c(-1,1),c(0,24),xout=empiricalK)
				
		barplot(Aprops,ylim=c(0,ymax),main="",col=c(colorCT),cex.axis=0.8,cex.main=1,las=1,yaxt="n")
		mtext(title,3,line=0.5,cex=0.7)
		mtext(mains[1],3,line=-2,cex=0.7)
		text(-1,ymax*0.9,letters[1],pos=4,cex=2)
		
		abline(v=RS.interp$y,col="red",lwd=1.5)
		
		if (yaxis==T) { axis(2,las=1,hadj=0.75,tcl=-0.25) 
		} else { axis(2,las=1,tcl=-0.25,labels=F)}
		if (type.label==T){ mtext("Int. Regime Shift",4,outer=F,cex=0.6) }
				
		IN<-Xct[[kendalls]][,i]
		OUT<-Xrs[[kendalls]][,i]
		if (min(IN)>max(OUT)) {
			text(-1,ymax*0.6,"No overlap",pos=4,cex=0.8)
		} else {
			roc<-calcROC(IN, OUT, max.len = 10000)
			if (roc$AUC<0.5){ 
				AUC<-1-round(roc$AUC,3)
			} else { AUC<-round(roc$AUC,3) }
			opt.interp<-approx(c(-1,1),c(0,24),xout=roc$optimal)
			auc<-paste(AUC,"+-",round(roc$se.fit,2))
			abline(v=opt.interp$y,lty=3,lwd=2,col="gray40")
			if (round(roc$optimal,3)<0){ xloc=-1
			} else { xloc=-1 }
			
		}
		polygon(x=c(2,2,25,25),y=c(0.35,0.5,0.5,0.35),col="white",border=NA)
		text(1,ymax*0.85,"Critical Transition",pos=4,cex=1)
						
		barplot(Cprops,ylim=c(0,ymax),main="",col=c(colorRS),cex.axis=0.8,cex.main=1,las=1,yaxt="n")
		mtext(mains[3],3,line=-1,cex=0.7)	
		text(-1,ymax*0.9,letters[3],pos=4,cex=2)	
		abline(v=RS.interp$y,col="red",lwd=1.5)
		if (yaxis==T) { axis(2,las=1,hadj=0.75,tcl=-0.25) 
		} else { axis(2,las=1,tcl=-0.25,labels=F) }
		if (type.label==T){ mtext("Ex. Regime Shift",4,outer=F,cex=0.6) }
		
		OUT<-Xrs[[kendalls]][,i]
		if (min(IN)>max(OUT)) {
			text(-1,ymax*0.6,"No overlap",pos=4,cex=0.8)
		} else {
			roc<-calcROC(IN, OUT, max.len = 10000)
			if (roc$AUC<0.5){ 
				AUC<-1-round(roc$AUC,3)
			} else { AUC<-round(roc$AUC,3) }
			opt.interp<-approx(c(-1,1),c(0,24),xout=roc$optimal)
			auc<-paste(AUC,"+-",round(roc$se.fit,2))
			abline(v=opt.interp$y,lty=3,lwd=2,col="gray40")
			if (round(roc$optimal,3)<0){ xloc=-1
			} else { xloc=-1 }
			#text(xloc,ymax*0.6,paste("AUC =", auc),pos=4,cex=0.8)
			#text(xloc,ymax*0.5,paste("optimum =", round(roc$optimal,3)),pos=4,cex=0.8)
		}
		polygon(x=c(2,2,25,25),y=c(0.35,0.5,0.5,0.35),col="white",border=NA)
		text(1,ymax*0.85,"Extrinsic Regime Shift",pos=4,cex=1)
	
		text(-1,ymax*0.5,paste("LWO Kt=", empiricalK,sep=""),pos=4,cex=1.1,col="red")
		text(xloc,ymax*0.7,paste("AUC=", auc,sep=""),pos=4,cex=1.1)
		text(xloc,ymax*0.6,paste("optimum=", round(roc$optimal,3),sep=""),pos=4,cex=1.1)
	
		temp.table<-round(t(apply(counts,1,function(x) x/sum(x))),digits=4)	
		mids<-barplot(t(temp.table),plot=FALSE)
		mstep<-(mids[2]-mids[1])/2
		locs<-mids-(mids[2]-mids[1])/2
		step<-locs[2]-locs[1]
		at.vect<-c(locs,locs[length(locs)]+step)
		locs3<-at.vect[seq(1,length(at.vect),2)]
		axis(1,at=locs3,labels=seq(-1,1,0.2),tick=T,las=2,cex=0.6,line=0.25,mgp=c(3,0.5,0.5),tcl=-0.25)	
	}
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
	
# Kt.summary.stats_________________________________________
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

# plot.supp_________________________________________
plot.supp<-function(Xct,Xdc,Xrs,Xnc,ind,main){
	colorCT<-rgb(0.1,0.3,0.4,1)
	colorDC<-rgb(0.1,0.3,0.4,0.7)
	colorRS<-rgb(0.1,0.3,0.4,0.5)
	colorNC<-rgb(0.1,0.3,0.4,0.3)
	dev.new(width=8,height=4.5)
	par(oma=c(6,4,4,1),mar=c(0.75,0.5,0,0.5))
	nf<-layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),nrow=4,ncol=4,byrow=F))

plot.taph.hists(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=T,mains=mains2,ymax=1.05,labs2=NULL,letters=c("a)","b)","c)","d)"),title="Untransformed",type.label=F,taph.ind=1)
plot.taph.hists(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=c("e)","f)","g)","h)"),title="Age-Depth",type.label=F,taph.ind=2)
plot.taph.hists(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=c("i)","j)","k)","l)"),title="AD+Even Samp",type.label=F,taph.ind=3)
plot.taph.hists(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=c("m)","n)","o)","p)"),title="AD+Targeted Samp",type.label=F,taph.ind=4)

	mtext(main,3,line=1.75,outer=T)
	mtext("Frequency",2,line=2,outer=T,cex=1.1)
	mtext("Kendall's tau",1,line=2.25,outer=T,cex=1.1)
}

# END FUNCTIONS ###########################################################


# RUN SIMULATIONS #########################################################

# Set Grass-Wood model parameters________________________
# parameters as of 30 Jan 2018
h = 0.5 # half saturation constant
r=0.25 #maximum tree growth rate, largely unused by the loops below. Values above about 1.5 are probably unrealistic.
delta_t = 1 # dt used for numerical simulations
gens = 10000 # length of the simulation, in non-dimensional time-steps
K_Start = 1 #initial value for the r parameter
K_Pulse_amt = -0.4 #amount that r changes during the driver pulse applied to the system
pulse_time = 1000 #length of time that the driver pulse is applied for
phi = 0.05
sigma_sd = 0.005
V0 = 1 # initial tree cover
beta_ps<-estBetaParams(mu=0.15, var=0.015)


# Set taphonomic parameters______________________________
exRStime=6000
nreps=100
nsamp=200
AC.buff=0.1
samp.freq2=0.4
steps<-c(1,1,1,1)
a0=2
a1=0.025
sd.pct=0.05
AC.samp=0.4

# Run time series iterations________________________________
# linear (5yr/cm) #
TAtop=5
TAbottom=5
agemodel="linearTA"
windows<-c(2500,600,50,50)

system.time(Xct1<-rep.ews(TStype="TSct",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0,a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xdc1<-rep.ews(TStype="TSdc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=1,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xrs1<-rep.ews(TStype="TSrs",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

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

# broken stick 1 ######
TAtop=5
TAbottom=20
breakpoint=2500
agemodel="brokenstick"
windows<-c(2500,400,50,50)

system.time(Xct3<-rep.ews(TStype="TSct",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0,a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xdc3<-rep.ews(TStype="TSdc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=1,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xrs3<-rep.ews(TStype="TSrs",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xnc3<-rep.ews(TStype="TSnc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

# broken stick 2 ######
TAtop=5
TAbottom=20
breakpoint=4000
agemodel="brokenstick"
windows<-c(2500,300,50,50)

system.time(Xct4<-rep.ews(TStype="TSct",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0,a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xdc4<-rep.ews(TStype="TSdc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=1,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xrs4<-rep.ews(TStype="TSrs",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.RS",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(Xnc4<-rep.ews(TStype="TSnc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff=5000,cutoff2=2000,trim.type="to.set.bounds",start=1000,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

# settings for Figures 4 and 5
colorCT<-rgb(0.1,0.3,0.4,1)
colorDC<-rgb(0.1,0.3,0.4,0.7)
colorRS<-rgb(0.1,0.3,0.4,0.5)
colorNC<-rgb(0.1,0.3,0.4,0.3)
mains2<-c("","","","")

# Plot Figure 4 (SD for each age model)__________________________________________
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

# Plot Figure 5 (AC for each age model)__________________________________________
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

# Plot Figure 6 (SD for age models and subsampling)____________________________________________________________
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


# Set up for plotting figures 2, 3______________________________
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

# Plot Figure 2______________________________
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
lines(1:gens,Knc,col="red",lty=3,lwd=2)
legend("bottomright",c("Tree cover","K parameter"),lty=c(1,2),col=c("gray40","red"),bty="n")

plot(TSrs,type="l",ylab="",ylim=c(-0.05,1.1),xaxt="n",xlab="",las=1,col="gray40",lwd=1.5,xlim=c(-300,10000),yaxt="n")
text(-350,1,"b)",cex=1.75)
axis(2,at=seq(0,1,0.5),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
lines(1:gens,Krs,col="red",lty=3,lwd=2)

plot(TSdc,type="l",ylab="",ylim=c(-0.05,1.1),xaxt="n",xlab="",las=1,col="gray40",lwd=1.5,xlim=c(-300,10000),yaxt="n")
text(-350,1,"c)",cex=1.75)
axis(2,at=seq(0,1,0.5),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
lines(1:gens,Kdc,col="red",lty=3,lwd=2)

plot(TSct,type="l",ylab="",xlab="time steps",ylim=c(-0.05,1.1),las=1,col="gray40",lwd=1.5,xlim=c(-300,10000),yaxt="n",xaxt="n")
text(-350,1,"d)",cex=1.75)
axis(2,at=seq(0,1,0.5),hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
axis(1,hadj=0.6,las=1,tcl=-0.25,cex.axis=1.5)
mtext("Time Steps (Years)",1,line=2,cex=1)
lines(1:gens,Kct,col="red",lty=3,lwd=2)

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

# Plot Figure 3______________________________
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
# Exponential model_____________________________________________________
# set taphonomic parameters
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

# set model paramters
 h = 0.5
 r=0.25 
delta_t = 1 
gens = 300 
K_Start = 1 
K_Pulse_amt = -0.4 
pulse_time = 3
phi = 0.05
sigma_sd = 0.005
V0 = 1 
beta_ps<-estBetaParams(mu=0.15, var=0.015) 

# Run simulations
# occasionally this x1 simulation fails. Just re-run
system.time(x1<-rep.ews(TStype="TSct",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0,a1,cutoff,cutoff2,trim.type="to.RS",start,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(x2<-rep.ews(TStype="TSdc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff,cutoff2,trim.type="to.set.bounds",start,q=1,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(x3<-rep.ews(TStype="TSrs",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff,cutoff2,trim.type="to.RS",start,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

system.time(x4<-rep.ews(TStype="TSnc",nreps=nreps, TAbottom=TAbottom,TAtop=TAtop,sample.by="distribute.samples",nsamp,samp.freq1=NULL,samp.freq2=NULL,AC.buff,windows=windows,steps=steps,agemodel=agemodel,breakpoint=breakpoint,det_method="gaussian",a0=a0,a1=a1,cutoff,cutoff2,trim.type="to.set.bounds",start,q=5,sd.pct=sd.pct,AC.samp=AC.samp))

Xct<-x1
Xdc<-x2
Xrs<-x3
Xnc<-x4

# plot Figure 7___________________________________________________________
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

############# WEST LAKE OKOBOJI ######## 
# Note that the Van Zant paper calls the site "Lake West Okoboji"
# But other sources, including Neotoma, use the name "West Lake Okoboji"

# download West Lake Okoboji data from Neotoma
siteX<-get_site("%Okoboji %")
siteX.dl<-get_download(siteX)
siteX.pol<-siteX.dl[[1]]$counts

# generate new bayesian age model for WLO
# skip this step if using previously-generated chronology
chronconts<-get_chroncontrol(siteX.dl)
chroncont.table<-chronconts[[1]]$chron.control
ages1<-chroncont.table$age
ageSds1<-(chroncont.table$age.old-chroncont.table$age)
positions1<-chroncont.table$depth
# Van Zant 1979 doesn't specify sample thickness
# We assume samples were 1cm thick
thickness1<-chroncont.table$thickness 
thickness1[is.na(thickness1)]<-1
calibcurves1<-(c("normal","normal",rep("intcal13",(length(ages1)-2))))
predictLocs<-siteX.dl[[1]]$sample.meta$depth

# agemodelk<-Bchronology(ages=ages1,ageSds=ageSds1,positions=positions1, positionThicknesses=thickness1,calCurves=calibcurves1,predictPositions=predictLocs,jitterPositions=T)

# chron<-summary(agemodelk)[,4]
# write.csv(summary(agemodelk),"WLOchronology.csv",rownames=F)

# read in WLO bayesian age model
chronfull<-read.csv("WLOchronology.csv")
chron<-chronfull[,4]

# limit to trees and upland herbs, than transform to percents
all_taxa <- do.call(rbind.data.frame, lapply(siteX.dl, function(x)x$taxon.list[,1:6]))
all_taxa <- all_taxa[!duplicated(all_taxa),]
good_cols<-c(which(colnames(siteX.pol) %in% all_taxa[all_taxa$ecological.group %in% c("TRSH","UPHE"),1]))
siteX.pol<-siteX.pol[,good_cols]
siteX.pct<-siteX.pol[,1:ncol(siteX.pol)]/rowSums(siteX.pol[,1:ncol(siteX.pol)],na.rm=TRUE)
	
# calculate percent arboreal pollen	
tree_cols<-c(which(colnames(siteX.pct) %in% all_taxa[all_taxa$ecological.group %in% c("TRSH"),1]))	
tree<-rowSums(siteX.pct[,tree_cols])

# Calculate depostion rates
depth<-siteX.dl[[1]]$sample.meta$depth
depth.difs<-depth[-1]-depth[-length(depth)]
age.difs<-chron[-1]-chron[-length(chron)]
dep.rates<-age.difs/depth.difs
plot(chron[-length(chron)],dep.rates)
dep.tab<-cbind(depth[-length(depth)],dep.rates)

# create a data table for LWO
data.tab<-cbind(rev(chron[-length(chron)]),rev(depth[-length(depth)]),rev(dep.rates),rev(tree[-length(tree)]))
data.tab<-data.tab[data.tab[,3]>0,]
colnames(data.tab)<-c("ages","depth","dep.rates","AP")

# find time of regime shift
bp.out<-CE.Normal.Mean(as.data.frame(tree),Nmax=1)
timeCT<-chron[bp.out$BP.Loc]

# find length of time prior to regime shift at LWO
yrs.preCT<-max(chron)-timeCT

# run rolling sd
win<-sum(chron>timeCT)/2
det.pol<-detrendTS(data.tab[,c(1,4)],method="gaussian")
tree.sd<-std.sd(det.pol,winlen=win,step.size=1)
preCTinds<-which(tree.sd$windows2[,2]>timeCT)
tree.sd.trim<-tree.sd$std.SD[preCTinds]
time.trim<-tree.sd$windows2[preCTinds,1]
time.temp<-max(time.trim)-time.trim #need to reverse order, so that we are moving forward in time
empiricalK<-cor(tree.sd.trim,time.temp,method="kendall")

# Simulations_____________________________________________________________
# set number of gens, and amt. of time prior to the regime shift based on observations at WLO
gens<-ceiling(max(data.tab[,1]))
amt.preCTtime<-max(chron)-timeCT

# set up vector of time averaging
interpTA<-approx(data.tab[,2],data.tab[,3],xout=seq(min(data.tab[,2]),max(data.tab[,2]),1))
TAbins<-interpTA$y
cTAbins<-c(1,cumsum(TAbins))
cTAbins<-cTAbins[which(cTAbins<gens)]

# set up sampling
sample.at<-data.tab[,2] #observed sampling at WLO

# set up 2x smpling
sample.at.mid<-c()
for(k in 1:(length(sample.at)-1)){
	step<-(sample.at[k]-sample.at[k+1])/2
	sample.at.mid[k]<-sample.at[k]-step
}
sample.at2<-sort(round(c(as.vector(sample.at),sample.at.mid)),decreasing=T) 

# Set model parameters__________________________
 h = 0.5 
 r=0.25 
 q=5
 c=1  
delta_t = 1 
r_Start = 1 
r_Pulse_amt = -0.3 
pulse_time = 1000 
sigma_sd = 0.005
V0 = 1 
FRI = 1
beta_ps<-estBetaParams(mu=0.15, var=0.015) 

# Set taphonomic parameters__________________________
exRStime=timeCT
nreps=100
steps<-c(1,1,1,1)

# run simulations for each time series type
Xct<-rep.ews.emp(TStype="TSct",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT)
Xdc<-rep.ews.emp(TStype="TSdc",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT)
Xrs<-rep.ews.emp(TStype="TSrs",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT)
Xnc<-rep.ews.emp(TStype="TSnc",nreps,TAbins,sample.at,sample.at2,detrend_method="gaussian",yrs.preCT)

# plot Figure 8_____________________________________________
ind="sd"
colorCT<-rgb(0.1,0.3,0.4,1)
colorDC<-rgb(0.1,0.3,0.4,0.7)
colorRS<-rgb(0.1,0.3,0.4,0.5)
colorNC<-rgb(0.1,0.3,0.4,0.3)
mains2<-c("","","","")

dev.new(width=5,height=7)
matx<-c(1,1,2,2,3,4,4,5,5,
1,1,2,2,3,6,6,7,7)
nf<-layout(matrix(matx,nrow=9,ncol=2,byrow=F))
par(oma=c(4,5,1,2),mar=c(1.5,1,0,1))

empiricalK<-round(empiricalK,digits=2)

preRS<-data.tab[data.tab[,1]>timeCT,]
plot(data.tab[,1],data.tab[,4],pch=16,xlim=c(16000,0),xaxt="n",las=1,cex.axis=1.5)
lines(data.tab[,1],data.tab[,4])
mtext("Prop. AP",2,line=3.25,cex=1.1)
mtext("a)",3,line=-1.6,cex=1.2,adj=0.015)
abline(v=timeCT,lty=3)

preRSsd<-tree.sd$std.SD[tree.sd$windows2[,2]>timeCT]
plot(tree.sd$windows2[tree.sd$windows2[,2]>timeCT,2],preRSsd,xlim=c(16000,0),pch=16,las=1,cex.axis=1.5)
lines(tree.sd$windows2[tree.sd$windows2[,2]>timeCT,2],preRSsd)
mtext("Std. SD",2,line=3.25,cex=1.1)
mtext("Calendar YBP",1,line=2.75,cex=1.1)
mtext("b)",3,line=-1.6,cex=1.2,adj=0.015)
abline(v=timeCT,lty=3)

plot(c(0,1),c(0,1),type="n",bty="n",xaxt="n",yaxt="n",ylab="",xlab="")

letters.list<-c("","","","")
plot.taph.hists.sub(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=T,mains=mains2,ymax=0.45,labs2=NULL,letters=letters.list,title="",type.label=F,taph.ind=3,empiricalK)
mtext("Frequency",2,line=2.75,outer=F,cex=1.1)   
mtext("Frequency",2,line=2.75,outer=F,cex=1.1,adj=7)   
mtext("Kendall's tau",1,line=3.5,cex=1.1)
mtext("c)",3,line=7.5,cex=1.2,adj=0.05)
mtext("Grad. CT",3,line=7.6,cex=0.8,adj=0.22)
mtext("d)",3,line=-1.9,cex=1.2,adj=0.05)
mtext("Abrupt CT",3,line=-1.8,cex=0.8,adj=0.25)
mtext("Observed Sampling",3,line=10,cex=1.2)

plot.taph.hists.sub(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=0.45,labs2=NULL,letters=letters.list,title="",type.label=F,taph.ind=4,empiricalK)
mtext("Kendall's tau",1,line=3.5,cex=1.1)
mtext("e)",3,line=7.5,cex=1.2,adj=0.05)
mtext("f)",3,line=-1.9,cex=1.2,adj=0.05)
mtext("Grad. CT",3,line=7.6,cex=0.8,adj=0.22)
mtext("Abrupt CT",3,line=-1.8,cex=0.8,adj=0.2)
mtext("2x Observed Sampling",3,line=10,cex=1.2)

# plot supplemental data: historgams of Kendall's tau values________________________________
dev.new(width=7,height=4.5)
par(oma=c(6,4,3,1),mar=c(0.75,0.5,0,0.5))
nf<-layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),nrow=4,ncol=4,byrow=F))

ind="sd"
letters.list<-c("a)","b)","c)","d)")
plot.taph.hists.emp(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=T,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Untransformed",type.label=F,taph.ind=1,empiricalK)
mtext("Frequency",2,line=2,outer=T,cex=0.7)  
letters.list<-c("e)","f)","g)","h)")
plot.taph.hists.emp(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Age-Depth",type.label=F,taph.ind=2,empiricalK)
mtext("Frequency",2,line=2,outer=T,cex=0.7)  
letters.list<--c("i)","j)","k)","l)")
plot.taph.hists.emp(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Age-Depth and Sampled",type.label=F,taph.ind=3,empiricalK)
mtext("Frequency",2,line=2,outer=T,cex=0.7)  
letters.list<-c("m)","n)","o)","p)")
plot.taph.hists.emp(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Age-Depth and 2x Sampled",type.label=F,taph.ind=4,empiricalK)
mtext("Frequency",2,line=2,outer=T,cex=0.7)  

dev.new(width=7,height=4.5)
par(oma=c(6,4,3,1),mar=c(0.75,0.5,0,0.5))
nf<-layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),nrow=4,ncol=4,byrow=F))

ind="ac"
letters.list<-c("a)","b)","c)","d)")
plot.taph.hists.emp(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=T,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Untransformed",type.label=F,taph.ind=1,empiricalK)
mtext("Frequency",2,line=2,outer=T,cex=0.7)  
letters.list<-c("e)","f)","g)","h)")
plot.taph.hists.emp(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Age-Depth",type.label=F,taph.ind=2,empiricalK)
mtext("Frequency",2,line=2,outer=T,cex=0.7)  
letters.list<-c("i)","j)","k)","l)")
plot.taph.hists.emp(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Age-Depth and Sampled",type.label=F,taph.ind=3,empiricalK)
mtext("Frequency",2,line=2,outer=T,cex=0.7)  
letters.list<-c("m)","n)","o)","p)")
plot.taph.hists.emp(Xct,Xdc,Xrs,Xnc,indicator=ind,yaxis=F,mains=mains2,ymax=1.05,labs2=NULL,letters=letters.list,title="Age-Depth and 2x Sampled",type.label=F,taph.ind=4,empiricalK)
mtext("Frequency",2,line=2,outer=T,cex=0.7)  

SuppWLO1<-Kt.summary.stats(Xct,Xdc,Xrs,Xnc,3,"sd")
SuppWLO2<-Kt.summary.stats(Xct,Xdc,Xrs,Xnc,3,"ac")

SuppWLO<-rbind(SuppWLO1[,c(1,2,3,5,6,8,9,10)],SuppWLO2[,c(1,2,3,5,6,8,9,10)])

#######################################
# Paleoecological Resilience Indicator functions
# Stegner et al. 
# updated 08 July 2018
#######################################

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

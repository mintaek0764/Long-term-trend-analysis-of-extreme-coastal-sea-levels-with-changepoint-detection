##########################################################################################
##########################################################################################
#####                                                                               ######
##### R Functions Used for the Analysis of the GESLA Data Set                         ######
#####                                                                               ######
##### Supplementary material to the "Long-term trend analysis of extreme coastal    ######
##### sea levels with changepoint detection"                                        ######
#####                                                                               ######
##### Last Updated on Monday, 7 December 2020                                       ######
#####                                                                               ######
##########################################################################################
##########################################################################################



##########################################################################################
##### Function that calculates number of non-missing observations between           ######
##### two locations (which could be changepoints)                                   ######
##########################################################################################



get.distance <- function(a, b, t.valid) {
  length(t.valid[which(t.valid >= a & t.valid < b)])
}



##########################################################################################
##### Function that calculates the penalty term in the MDL for the specified        ######
##### changepoint configuration                                                     ######
##########################################################################################



get.penalty <- function(cp, n, t.valid) {
  # tau_0 = 1 and tau_(k+1) = N+1
  tau <- c(1, cp, n+1)
  # k is number of changepoints
  k <- length(tau) - 2
  # no penalty assigned when there is no changepoint
  if(k == 0) {
	Penalty <- 0
  # penalty for changepoints
  } else {
  # calculating tau_i - tau_{i-1} where i = 2,3,..,k+1 by only considering non-missing values
	diff.tau <- numeric()
	# index here is 3:k+2 because R's index starts from 1 (it is in fact working from i=2 to k+1)
	for (i in 3:(k+2)) {
		diff.tau[i-2] <- get.distance(tau[i-1], tau[i], t.valid)
	}
	# penalty term for changepoints
	Penalty <- sum(log2(diff.tau))/2 + log2(k+1) + sum(log2(tau[3:(k+2)]))
  }
	
  # add log2(N)/2 penalty term to account for the long-term linear trend parameter, alpha
  Penalty <- Penalty + log2(length(t.valid))/2

  # returning the penalty term
  Penalty
}



##########################################################################################
##### Function that fits the GEV model with specified changepoints                  ######
##### for use in GESLA data analysis                                                ######
##########################################################################################



# gev optimization function (using periodic location parameters)
fit.gev.changepoint <- function(data, cp) {
  # Count how many observations are there
  t <- 1:length(data)
  # Periodicity
  T <- 12
  # total number of observations
  n <- length(t)

  # times that observations were missing
  na.locations <- t[is.na(data)]

  # times that can be chosen as potential changepoints
  # disallowing first year and last year of data to be chosen as changepoints
  t.valid <- t
  t.valid[c(1:T,(n-T+1):n)] <- 0
  # disallowing times when data are not available to be chosen as changepoints
  t.valid <- t.valid[!is.na(data)]
  t.valid <- t.valid[t.valid > 0]

  # calculating cosine and sine values to be used in GEV fitting in advance
  cos.vals <- cos(2*pi*t/T)
  sin.vals <- sin(2*pi*t/T)
  
  # calculating cosine and sine values to be used in GEV fitting in advance
  cos.vals2 <- cos(2*pi*t/(T/2))
  sin.vals2 <- sin(2*pi*t/(T/2))
  
  # target function for the GEV fitting
  gev.lik.optim <- function(a) {
    # a is a vector containing all parameters (mu, sigma, xi, alpha)
    m <- a[1:5]
    c <- a[6:10]
    xi <- a[11]
	alpha <- a[12]
	min_num_par <- 12
    
    # periodic scale parameter (sigma)
    sig <- c[1] + c[2]*cos.vals + c[3]*sin.vals + c[4]*cos.vals2 + c[5]*sin.vals2
    
    # mean shifts
    delta <- rep(0,length(data))
    # changepoint(s), if any
    if(length(a) > min_num_par) {
      for (i in 1:length(cp)) {
        delta[cp[i]:n] <- a[min_num_par +i]
      }
    }
    
    # periodic location parameter base
    mu0 <- m[1] + m[2]*cos.vals + m[3]*sin.vals + m[4]*cos.vals2 + m[5]*sin.vals2
        
    # location parameter with all necessary components
    mu <- mu0 + alpha*t/T/10 + delta
    y <- 1 + xi * {data - mu}/sig
    
    # exclude values when data is missing
    sig <- sig[!is.na(data)]
    y <- y[!is.na(data)]
    
    # calculating negative log-likelihood value for GEV models
    if (any(y <= 0) || any(sig <= 0)) {
      negloglik <- 10^12
    } else {
      negloglik <- sum(log(sig)) + sum(y^{-1/xi}) + sum(log(y) * {1/xi + 1})
    }
    # return the negative log-likelihood value
    negloglik
  }
  
  # a series containing corresponding months for each observation
  months <- rep(1:T,length.out=n)
  
  # intiial values for mean shifts
  deltainit <- rep(0,length.out=length(cp))
  cp_full <- c(1,cp,length(data)+1)
  if(length(deltainit) > 0) {
    for(i in 1:length(cp)) {
    	# number of months available before the changepoint
    	month_length_front <- cp_full[i+1]-cp_full[i]
    	# number of months available after the changepoint
    	month_length_back <- cp_full[i+2]-cp_full[i+1]
    	# determine how many months to use for estimating delta initial value
    	month_length <- min(24,month_length_front,month_length_back)
    	# if earlier batch has less months, use immediately preceding 24 months from that batch
    	if(month_length_front <= month_length_back) {
    		# take 24 (or less depending on availability) observations directly preceding the changepoint
    		front_batch_index <- (cp_full[i+1]-month_length):(cp_full[i+1]-1)
    		front_batch_months <- months[(cp_full[i+1]-month_length):(cp_full[i+1]-1)]
    		# find out how many months that need to be added to align the front batch and back batch
    		back_batch_alignment <- (front_batch_months[1] - months[cp_full[i+1]]) %% 12
    		back_batch_index <- (cp_full[i+1]+back_batch_alignment):(cp_full[i+1]+back_batch_alignment+month_length-1)
    		back_batch_index_corrected <- back_batch_index[back_batch_index %in% t]
    		back_batch <- data[back_batch_index_corrected]
       		# resample front batch using corrected indeces
       		front_batch_index_corrected <- front_batch_index[back_batch_index %in% t]
       		front_batch <- data[front_batch_index_corrected]
    	} else {
			# take 24 (or less depending on availability) observations directly following the changepoint
			back_batch_index <- (cp_full[i+1]):(cp_full[i+1]+month_length-1)
    		back_batch_months <- months[(cp_full[i+1]):(cp_full[i+1]+month_length-1)]
    		# find out how many months that need to be added to align the front batch and back batch
    		front_batch_alignment <- (months[cp_full[i+1]-month_length] - back_batch_months[1]) %% 12
    		front_batch_index <- (cp_full[i+1]-front_batch_alignment-month_length):(cp_full[i+1]-1-front_batch_alignment)
    		front_batch_index_corrected <- front_batch_index[front_batch_index %in% t]
    		front_batch <- data[front_batch_index_corrected]
       		# resample back batch using corrected indeces
       		back_batch_index_corrected <- back_batch_index[front_batch_index %in% t]
       		back_batch <- data[back_batch_index_corrected]
    	}
    	deviations <- back_batch - front_batch
    	avg_deviation <- mean(deviations,na.rm=TRUE)
    	avg_deviation[is.na(avg_deviation)] <- 0
    	deltainit[i] <- avg_deviation
    }
    deltainit <- cumsum(deltainit)
  }

  data_shifted <- data
  if(length(deltainit) > 0) {
    for(i in 2:(length(cp)+1)) {
      data_shifted[(cp_full[i]):(cp_full[i+1]-1)] <- data[(cp_full[i]):(cp_full[i+1]-1)] - deltainit[i-1]
    }
  }
  
  # initial value for shape parameter, xi
  # four initial values will be used
  xiinit_neg2 <- -0.15
  xiinit_neg1 <- -0.1
  xiinit_pos1 <- 0.1
  xiinit_pos2 <- 0.15
  
  # initial values for periodic scale parameter, sigma
  cinit <- c(sqrt(6*max(tapply(data_shifted,months,var,na.rm=TRUE)))/pi, 0, 0, 0, 0)
  # initial value for mu0
  minit <- c(mean(data_shifted,na.rm=TRUE)-0.57722*sqrt(6*var(data_shifted,na.rm=TRUE))/pi, 0, 0, 0, 0)
  
  # initiap value for the linear trend term
  alphainit <- 0
  
  # all initial values
  init_xineg2 <- c(minit, cinit, xiinit_neg2, alphainit, deltainit)
  init_xineg1 <- c(minit, cinit, xiinit_neg1, alphainit, deltainit)
  init_xipos1 <- c(minit, cinit, xiinit_pos1, alphainit, deltainit)
  init_xipos2 <- c(minit, cinit, xiinit_pos2, alphainit, deltainit)
  
  # fitting the GEV model
  fit_xineg2 <- optimr(init_xineg2, gev.lik.optim, hessian = TRUE, method = "BFGS", control = list(maxit = 10000))
  fit_xineg1 <- optimr(init_xineg1, gev.lik.optim, hessian = TRUE, method = "BFGS", control = list(maxit = 10000))
  fit_xipos1 <- optimr(init_xipos1, gev.lik.optim, hessian = TRUE, method = "BFGS", control = list(maxit = 10000))
  fit_xipos2 <- optimr(init_xipos2, gev.lik.optim, hessian = TRUE, method = "BFGS", control = list(maxit = 10000))

  # picking the GEV model with lower neg log lik value
  nnloglik <- c(fit_xineg2$value, fit_xineg1$value, fit_xipos1$value, fit_xipos2$value)

  best_fit <- which(nnloglik == min(nnloglik))
  
  if(best_fit == 1) {
  	fit <- fit_xineg2
  } else if (best_fit == 2) {
  	fit <- fit_xineg1
  } else if (best_fit == 3) {
  	fit <- fit_xipos1
  } else {
  	fit <- fit_xipos2
  }

  fit$MDL <- fit$value/log(2) + get.penalty(cp, n, t.valid)
  
  # returning the MLE model
  fit
}



##########################################################################################
##### Function that fits the GEV model with specified changepoints                  ######
##### for use in bootstrap                                                          ######
##########################################################################################



fit.gev.changepoint.bootstrap <- function(resample, cp_original, n_original, delta_original, boot_index) {
  # Count how many observations are there
  t <- boot_index
  # Periodicity
  T <- 12
  # total number of observations
  n <- length(t)

  # calculating cosine and sine values to be used in GEV fitting in advance
  cos.vals <- cos(2*pi*t/T)
  sin.vals <- sin(2*pi*t/T)
  
  # calculating cosine and sine values to be used in GEV fitting in advance
  cos.vals2 <- cos(2*pi*t/(T/2))
  sin.vals2 <- sin(2*pi*t/(T/2))

  # target function for the GEV fitting
  gev.lik.optim <- function(a) {
    # a is a vector containing all parameters
    m <- a[1:5]
    c <- a[6:10]
    xi <- a[11]
    alpha <- a[12]
    
    # periodic scale parameter (sigma)
    sig <- c[1] + c[2]*cos.vals + c[3]*sin.vals + c[4]*cos.vals2 + c[5]*sin.vals2
    
    # mean shifts
    delta <- rep(0,length(resample))
    # changepoint(s), if any
    if(length(a) > 12) {
      for (i in 1:length(cp_original)) {
      	cp_index <- cp_original[i]:n_original
        delta[boot_index %in% cp_index] <- a[12+i]
      }
    }
    
    # periodic location parameter base
    mu0 <- m[1] + m[2]*cos.vals + m[3]*sin.vals + m[4]*cos.vals2 + m[5]*sin.vals2

    # location parameter with all necessary components
    mu <- mu0 + alpha*t/T/10 + delta
    y <- 1 + xi * {resample - mu}/sig
    
    # exclude values when data is missing
    sig <- sig[!is.na(resample)]
    y <- y[!is.na(resample)]
    
    # calculating negative log-likelihood value for GEV models
    if (any(y <= 0) || any(sig <= 0)) {
      negloglik <- 10^12
    } else {
      negloglik <- sum(log(sig)) + sum(y^{-1/xi}) + sum(log(y) * {1/xi + 1})
    }
    # return the negative log-likelihood value
    negloglik
  }
  
  # a series containing corresponding months for each observation
  months <- rep(1:T,length.out=n_original)
  months <- months[boot_index]
  
  # intiial values for mean shifts
  deltainit <- unique(delta_original)[-1]

  # remove the mean shifts from the resampled data
  resample_shifted <- resample - delta_original[boot_index]
  
  # initial value for shape parameter (xi)
  xiinit_neg2 <- -0.15
  xiinit_neg1 <- -0.1
  xiinit_pos1 <- 0.1
  xiinit_pos2 <- 0.15
  
  # initial values for periodic scale parameter (sigma)
  cinit <- c(sqrt(6*max(tapply(resample_shifted,months,var,na.rm=TRUE)))/pi, 0, 0, 0, 0)
  # initial value for mu0
  minit <- c(mean(resample_shifted,na.rm=TRUE)-0.57722*sqrt(6*var(resample_shifted,na.rm=TRUE))/pi, 0, 0, 0, 0)

  # initial value for the linear trend term
  alphainit <- 0
  
  # all initial values
  init_xineg2 <- c(minit, cinit, xiinit_neg2, alphainit, deltainit)
  init_xineg1 <- c(minit, cinit, xiinit_neg1, alphainit, deltainit)
  init_xipos1 <- c(minit, cinit, xiinit_pos1, alphainit, deltainit)
  init_xipos2 <- c(minit, cinit, xiinit_pos2, alphainit, deltainit)
  
  # fitting the GEV model
  fit_xineg2 <- optimr(init_xineg2, gev.lik.optim, hessian = TRUE, method = "BFGS", control = list(maxit = 10000))
  fit_xineg1 <- optimr(init_xineg1, gev.lik.optim, hessian = TRUE, method = "BFGS", control = list(maxit = 10000))
  fit_xipos1 <- optimr(init_xipos1, gev.lik.optim, hessian = TRUE, method = "BFGS", control = list(maxit = 10000))
  fit_xipos2 <- optimr(init_xipos2, gev.lik.optim, hessian = TRUE, method = "BFGS", control = list(maxit = 10000))

  # picking the GEV model with lower neg log lik value
  nnloglik <- c(fit_xineg2$value, fit_xineg1$value, fit_xipos1$value, fit_xipos2$value)

  best_fit <- which(nnloglik == min(nnloglik))
  
  if(best_fit == 1) {
  	fit <- fit_xineg2
  } else if (best_fit == 2) {
  	fit <- fit_xineg1
  } else if (best_fit == 3) {
  	fit <- fit_xipos1
  } else {
  	fit <- fit_xipos2
  }
  
  # returning the MLE model
  fit
}



##########################################################################################
##### Calculate return exceedances from the median sea level.                       ######
##########################################################################################



GEV_return_exceedances_extremal <- function(data, extremal, cp_station, median_SL, last_month, years = c(25,50,75,100), boot_index = numeric(0)) {

	# time index (months)
	if (length(boot_index) == 0) {
		t <- 1:length(data)
		gev_MLE <- fit.gev.changepoint(data, cp_station)
	} else {
		t <- boot_index
		resample <- data[t]
		
		if(length(cp_station) > 0) {
			mean_shifts_original <- fit.gev.changepoint(data, cp_station)$par[13:(12 + length(cp_station))]
			# mean shifts
			delta_original <- rep(0, length(data))
			# changepoint(s), if any
			if(length(mean_shifts_original) > 0) {
			  for (i in 1:length(cp_station)) {
			    delta_original[cp_station[i]:length(data)] <- mean_shifts_original[i]
			  }
			}
			
			gev_MLE <- fit.gev.changepoint.bootstrap(resample, cp_station, length(data), delta_original, boot_index)
		} else {
			delta_original <- rep(0,length(data))
			gev_MLE <- fit.gev.changepoint.bootstrap(resample, cp_station, length(data), delta_original, boot_index)
		}
		
	}
	
	# save estimated GEV parameters
	params <- gev_MLE$par
	
	# index for January 2020 (used in return level calculations)
	t_2020_01 <- length(data) + (2641 - last_month)
	
	### Calculating return levels with extremal index
	rl_extremal <- numeric(0)
	for(i in 1:length(years)) {
		rl_year <- years[i]
		value <- round(min(data, na.rm = TRUE),3) - 0.001
		cumprob <- 100
		while(cumprob >= 1) {
			value <- value + 0.001
			probs <- numeric(0)
			for(a in 1:(rl_year*12))	{
				probs[a] <- GEV_CCDF(value, t_2020_01 + a - 1, params, extremal)
			}
			cumprob <- sum(probs)
		}
		rl_extremal[i] <- value
	}
	re_extremal <- rl_extremal - median_SL
	
	### Calculating return levels without extremal index
	rl_no_extremal <- numeric(0)
	for(i in 1:length(years)) {
		rl_year <- years[i]
		value <- round(min(data, na.rm = TRUE),3) - 0.001
		cumprob <- 100
		while(cumprob >= 1) {
			value <- value + 0.001
			probs <- numeric(0)
			for(a in 1:(rl_year*12))	{
				probs[a] <- GEV_CCDF(value, t_2020_01 + a - 1, params, 1)
			}
			cumprob <- sum(probs)/12
		}
		rl_no_extremal[i] <- value
	}
	re_no_extremal <- rl_no_extremal - median_SL

	# returning estimated return exceedances, linear trend, and others
	list(return_level_no_extremal=rl_no_extremal, return_level_extremal=rl_extremal, return_exceedance_no_extremal=re_no_extremal, return_exceedance_extremal=re_extremal, trend=params[12])
}



##########################################################################################
##### Calculate the GEV complementary CDF value for a given value.                  ######
##########################################################################################



GEV_CCDF <- function(value, time, params, extremal) {
	# retrieve parameter estimates
	m <- params[1:5]
    c <- params[6:10]
    xi <- params[11]
    alpha <- params[12]
    if(length(params) > 12) {
    	delta <- params[length(params)]
    } else {
    	delta <- 0
    }
    
    # Periodicity
    T <- 12
    
    # calculating GEV parameters before extremal index is applied
    sig <- c[1] + c[2]*cos(2*pi*time/T) + c[3]*sin(2*pi*time/T) + c[4]*cos(2*pi*time/(T/2)) + c[5]*sin(2*pi*time/(T/2))
    mu0 <- m[1] + m[2]*cos(2*pi*time/T) + m[3]*sin(2*pi*time/T) + m[4]*cos(2*pi*time/(T/2)) + m[5]*sin(2*pi*time/(T/2))
    mu <- mu0 + alpha*time/T/10 + delta
    
    # calculating GEV parameters with extremal index
    mu_ast <- mu - (sig/xi)*(1 - extremal^(xi))
    sig_ast <- sig*(extremal^xi)
    
    # calculate the cdf value for the corresponding value
    yz <- max(0, 1 + xi*((value - mu_ast)/sig_ast))
    ccdf <- 1 - exp(-yz^(-1/xi))
    
    # return the cdf value
    ccdf
}



##########################################################################################
##### A function that makes the quantile-quantile plot using Gumbel scale           ######
##########################################################################################



quantile_gumbel <- function(data, cp_station) {
	# the user can supply "numeric(0)" for the cp_station if they wish to obtain the QQ plot without changepoints

	gev_MLE <- fit.gev.changepoint(data, cp_station)
	params <- gev_MLE$par

	### QQ Plot ###

	sorted_data <- sort(data)


	# Count how many observations are there
	t <- 1:length(data)
	# Periodicity
	T <- 12
	# total number of observations
	n <- length(t)

	# calculating cosine and sine values to be used in GEV fitting in advance
	cos.vals <- cos(2*pi*t/T)
	sin.vals <- sin(2*pi*t/T)

	# calculating cosine and sine values to be used in GEV fitting in advance
	cos.vals2 <- cos(2*pi*t/(T/2))
	sin.vals2 <- sin(2*pi*t/(T/2))

	# sorting out parameters
	
	m <- params[1:5]
	c <- params[6:10]
	xi <- params[11]
	alpha <- params[12]

	# defining sigma
	sig <- c[1] + c[2]*cos.vals + c[3]*sin.vals + c[4]*cos.vals2 + c[5]*sin.vals2

	# mean shifts
	delta <- rep(0,length(data))
	# changepoint(s), if any
	if(length(params) > 12) {
	  for (i in 1:length(cp_station)) {
	    delta[cp_station[i]:n] <- params[12+i]
	  }
	}
    
	# periodic location parameter base
	mu0 <- m[1] + m[2]*cos.vals + m[3]*sin.vals + m[4]*cos.vals2 + m[5]*sin.vals2

	# location parameter with all necessary components
	mu <- mu0 + alpha*t/T/10 + delta	
	
	##### Transform data to Gumbel distribution

	gumbel <- xi^(-1)*log(1+xi*(data-mu)/sig)
	gumbel <- sort(gumbel)
	
	quantiles <- 1:length(data)/(length(data) + 1)
	quantiles <- quantiles[!is.na(data)]
	
	x_qq <- -log(-log(quantiles))
	
	##### Drawing a quantile Plot
	if(length(cp_station) == 0) {
		displayname <- "Quantile-Quantile Plot (without changepoints)"
	} else {
		displayname <- "Quantile-Quantile Plot (with changepoints)"
	}
	plot(x_qq, gumbel,ylab = "Sample Quantiles", xlab = "Theoretical Quantiles", main = displayname)
	abline(0, 1, col = 4)
}



##########################################################################################
##### Implement the moving Block Bootstrap resampling.                              ######
##########################################################################################



moving.block <- function(size,num.total,seed) {
	set.seed(seed)
	num.blocks <- num.total - size + 1
	blocks.out <- floor(num.total/size)
	block.index <- sample(1:num.blocks,size=blocks.out,replace=TRUE)
	index <- numeric()
	for(i in 1:blocks.out) {
		index[(1+size*(i-1)):(size*i)] <- seq(from=block.index[i],to=block.index[i]+size-1)
	}
	index
}



##########################################################################################
##### End of File                                                                   ######
##########################################################################################
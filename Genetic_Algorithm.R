##########################################################################################
##########################################################################################
#####                                                                               ######
##### Genetic Algorithm for Detecting Changepoints in the GESLA Monthly Maximum     ######
##### Sea Level Series                                                              ######
#####                                                                               ######
##### Supplementary material to the "Long-term trend analysis of extreme coastal    ######
##### sea levels with changepoint detection"                                        ######
#####                                                                               ######
##### Last Updated on Monday, 7 December 2020                                       ######
#####                                                                               ######
##########################################################################################
##########################################################################################



##########################################################################################
##### Reading Data to Run the GA Algorithm                                          ######
##########################################################################################



# Taking the argument passed from the command line
args <- commandArgs(trailingOnly=TRUE)
# First (and the only) argument is GESLA tide gauge id, ranging from 1 to 300
tide_gauge_ID <- as.numeric(args[1])

### Loading pre-defined custom R functions
source("Functions.R")

# Read the GESLA data set
GESLA_dataset_name <- paste('MonthlyMax/SLstation', tide_gauge_ID, '.csv', sep='')
data_read <- read.csv(GESLA_dataset_name)
months <- data_read[,1]
data <- data_read[,2]
rm(data_read)

message("Processing for GESLA Station #", tide_gauge_ID)
message("")

# seed value to be used for random number generation
seed_value <- 20200805

## Loading required libraries
library(foreach)
library(doParallel)
library(optimr)

## setup parallel backend to use many processors
cores <- detectCores()
cluster <- makeCluster(cores[1] - 1)
registerDoParallel(cluster)

message("Using ", cores[1], " cores.")
message("")



##########################################################################################
##### Setting Up Variables to be Used for GEV Estimation                            ######
##########################################################################################



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
cos.vals2 <- cos(2*pi*t/(T/2))
sin.vals2 <- sin(2*pi*t/(T/2))



##########################################################################################
##### Setting Up Variables to be Used for Genetic Algorithm                         ######
##########################################################################################



# how many children gets produced in each generation
gen.size <- 200

# for the initial generation only, limit how many changepoints can be generated per each parent
min.initial.cp <- 0
max.initial.cp <- 5

# maximum number of generations to run
max.gen <- 300

# deciding how far each changepoint have to be to be considered to be distinctive
# if more than one changepoints are chosen within that distance, they will be merged
merge.distance <- T # equivalent to one year apart
# probability of mutation
mutation.prob <- 0.002

# runs without improvements on MDL before being considered to have converged (not currently used)
max.runs <- max.gen



##########################################################################################
##### Start of Genetic Algorithm for Changepoint Detection                          ######
##### The algorithm will be run considering the linear trend.                       ######
##########################################################################################



# to calculate how long it took to process the entire algorithm
beginning.time <- proc.time()
message("Genetic algorithm for changepoint detection has started running")
message("")
message("Each generation will have ", gen.size, " children")
message("Genetic algorithm will run for ", max.gen, " generations")
message("Mutation probability is ", mutation.prob)
message("")

# counter for the runs without improvements
run <- 0

# declaring/emptying storage variables
# This matrix stores all MDL values
MDL <- matrix(0,max.gen,gen.size)
# This list stores all children
history <- list()
# This list stores the best from each generation only
best.fit <- list()
best.MDL <- numeric(max.gen)



##########################################################################################
##### Initializing the First Generation                                             ######
##########################################################################################



# setting beginning time
initial.start <- proc.time()
# setting seeds for random number generations
set.seed(seed_value)
# temporary list to store current changepoints generated
cp_list <- list()
# temporary list to store current configurations
configs <- list()

# generating one child at a time for the initial generation
i <- 1
while(i <= gen.size) {	
	# randomly generating changepoint occurence times
	cp <- sort(sample(t.valid, sample(min.initial.cp:max.initial.cp, 1)), method="quick")
	# merging changepoints that are within 20 observations
	if(length(cp) > 1) {
		compare <- cp[1]
		merged <- numeric(0)
		cp.to.merge <- cp[1]
		for(p in 2:length(cp)) {
			if(get.distance(compare, cp[p], t.valid) <= merge.distance) {
				cp.to.merge <- c(cp.to.merge,cp[p])
			} else {
				merged <- c(merged,sum(cp.to.merge)/length(cp.to.merge))
				compare <- cp[p]
				cp.to.merge <- cp[p]
			}
		}
		cp <- round(c(merged,sum(cp.to.merge)/length(cp.to.merge)))
	}
	# checking if the produced child has already been produced in the current generation
	duplicate <- FALSE
	if(i > 1){
		for (k in 1:(i-1)) {
			if(all(cp %in% configs[[k]])) {
				duplicate <- TRUE
				break
			}
		}
	}
	# repeating the process if produced child is a duplicate
	if(duplicate) {
		# discard and repeat this iteration
		i <- i
	# proceed otherwise
	} else {
		configs[[i]] <- cp
		i <- i + 1
	}
}

# calculating MDL value using foreach function
out <- foreach(i = 1:gen.size,.packages="optimr") %dopar% {
	fit.gev.changepoint(data, configs[[i]])$MDL
}
out <- unlist(out)
# ensures all MDL values are finite
out[is.infinite(out)] <- 10^12
# saving MDL values to a matrix
MDL[1,] <- out

# archiving current configuration
history[[1]] <- configs
# determining which child was the fittest in terms of MDL scores
best.where <- which(MDL[1,] == min(MDL[1,]))
# if there were ties in the fittest child, randomly select one of them
if (length(best.where) > 1) {
	best.where <- sample(best.where, 1)
}
# archiving best fit and corresponding MDL and trend
best.fit[[1]] <- configs[[best.where]]
best.MDL[1] <- MDL[1,best.where]

# calculating how long it took to finish to intiate the first generation
initial.end <- proc.time() - initial.start
message("Took ", as.numeric(initial.end[3]), " seconds for generation ", 1)
message("Best fit for this generation was ")
print(best.fit[[1]])
message("with MDL of ")
print(best.MDL[[1]])
message(" ")



##########################################################################################
##### Initializing the Subsequent Generations                                       ######
##########################################################################################



# setting seeds for random number generations
set.seed(seed_value)

# running the algorithm from the second generation until maximum generation number allowed
for (g in 2:max.gen) {
	# setting beginning time
	ptm <- proc.time()
	message("Started processing generation ", g)
	# temporary list to store current configurations
	configs <- list()
 	# Assigning which potential parent corresponds to rank 1 through rank 'gen.size'
	# rank 'gen.size' corresponds to the fittest and rank 1 corresponds to the least fittest
	index <- sort(MDL[g-1,], decreasing=TRUE, method="quick", index=TRUE)$ix

	##########################################################################################
	##### Initializing distinctive child using different parents begins here            ######
	##########################################################################################
	
	# The best changepoint configuration from the previous generation is passed over unaltered
	configs[[1]] <- best.fit[[g-1]]
	# Second best changepoint configuration from the previous generation is also passed over unaltered
	configs[[2]] <- unlist(history[[g-1]][order(MDL[g-1,])[2]])
	
	# crossing the two best changepoint configurations to get the third child
	child.produced <- sort(unique(c(configs[[1]],configs[[2]])), method="quick")
	thinned <- child.produced[runif(length(child.produced)) > 0.5]
	thinned <- thinned + sample(c(-1,0,1), size=length(thinned), prob=c(0.3,0.4,0.3), replace=TRUE)
	mutation.pool <- t.valid[!t.valid %in% thinned]
	mutation <- mutation.pool * sample(c(0,1), length(mutation.pool), prob=c(1-mutation.prob,mutation.prob), replace=TRUE)
	thinned <- c(thinned,mutation[mutation > 0])
	# Disallowing changepoints that are not possible
	# if resulted in one of missing areas
	if(any(thinned %in% na.locations)) {
		for(j in 1:length(thinned)) {
			if(thinned[j] %in% na.locations) {
				thinned[j] <- min(t.valid[t.valid > thinned[j]],n)
			}
		}
	}
	# disallowing first or last observation as possible changepoints
	child <- thinned[which(thinned > 1 & thinned < n)]
    
    # sorting changepoints
	child <- sort(unique(child), method="quick")

	# merging changepoints that are within a pre-specified number of observations
	if(length(child) > 1) {
		compare <- child[1]
		merged <- numeric(0)
		cp.to.merge <- child[1]
		for(p in 2:length(child)) {
			if(get.distance(compare, child[p], t.valid) <= merge.distance) {
				cp.to.merge <- c(cp.to.merge,child[p])
			} else {
				merged <- c(merged,sum(cp.to.merge)/length(cp.to.merge))
				compare <- child[p]
				cp.to.merge <- child[p]
			}
		}
		child <- round(c(merged,sum(cp.to.merge)/length(cp.to.merge)))
	}
	configs[[3]] <- child
	
	# generate the remaining children
	i <- 4
	while(i <= gen.size) {
		# Choosing father from the previous generation
		father.rank <- numeric(gen.size)
		father.rank[index] <- 1:gen.size
		# assigning probability of getting chosen based on their rank
		father.prob <- father.rank/sum(father.rank)
		# selecting the father based on that probability
		father.index <- sample(1:gen.size, size=1, prob=father.prob)
		father <- history[[g-1]][[father.index]]

		# Choosing mother from the previous generation excluding father
		mother.rank <- father.rank[-father.index]
		mother.rank[mother.rank > father.rank[father.index]] <- mother.rank[mother.rank > father.rank[father.index]] - 1
		# assigning probability of getting chosen based on their rank
		mother.prob <- mother.rank/sum(mother.rank)
		# selecting the mother based on that probability
		mother.index <- sample(1:c(gen.size-1), size=1, prob=mother.prob)
		if(mother.index >= father.index) {
			mother.index <- mother.index + 1
		}
		mother <- history[[g-1]][[mother.index]]

		# Choosing child
		child.produced <- sort(unique(c(father,mother)), method="quick")

		# thinning each changepoint with 0.5 probability
		thinned <- child.produced[runif(length(child.produced)) > 0.5]
		# if resulted in empty set, assign 0 instead
		if(length(thinned) < 1) {
			thinned <- 1
		} else {
			thinned <- thinned + sample(c(-1,0,1), size=length(thinned), prob=c(0.3,0.4,0.3), replace=TRUE)
		}

		# Allowing mutations with a preset mutation probability for each changepoint that is not present
		mutation.pool <- t.valid[!t.valid %in% thinned]
		mutation <- mutation.pool * sample(c(0,1), length(mutation.pool), prob=c(1-mutation.prob,mutation.prob), replace=TRUE)
		thinned <- c(thinned,mutation[mutation > 0])

		# Disallowing changepoints that are not possible
		# if resulted in one of missing areas
		if(any(thinned %in% na.locations)) {
			for(j in 1:length(thinned)) {
				if(thinned[j] %in% na.locations) {
					thinned[j] <- min(t.valid[t.valid > thinned[j]], n)
				}
			}
		}
		# disallowing first or last observation as possible changepoints
		child <- thinned[which(thinned > 1 & thinned < n)]
    
		# sorting changepoints
		child <- sort(unique(child), method="quick")

		# merging changepoints that are within a pre-specified number of observations
		if(length(child) > 1) {
			compare <- child[1]
			merged <- numeric(0)
			cp.to.merge <- child[1]
			for(p in 2:length(child)) {
				if(get.distance(compare, child[p], t.valid) <= merge.distance) {
					cp.to.merge <- c(cp.to.merge,child[p])
				} else {
					merged <- c(merged,sum(cp.to.merge)/length(cp.to.merge))
					compare <- child[p]
					cp.to.merge <- child[p]
				}
			}
			child <- round(c(merged,sum(cp.to.merge)/length(cp.to.merge)))
		}

		# checking if the produced child has already been produced in the current generation
		duplicate <- FALSE
		if(i > 1) {
			for (k in 1:(i-1)) {
				if(all(child %in% configs[[k]])) {
					duplicate <- TRUE
					break
				}
			}
		}

		# repeating the process if produced child is a duplicate
		if(duplicate) {
			# discard and repeat the iteration
			i <- i
		# proceed otherwise
		} else {
			# archive the child
			configs[[i]] <- child
			i <- i + 1
		}
	}

	# calculating MDL value in parallel using foreach function
	out <- foreach(i = 1:gen.size) %dopar% {
		fit.gev.changepoint(data, configs[[i]])$MDL
	}
	out <- unlist(out)
	# this ensures all MDL values are finite
	out[is.infinite(out)] <- 10^12
	# saving MDL values to a matrix
	MDL[g,] <- out



	##########################################################################################
	##### Initializing distinctive child using different parents ends here              ######
	##########################################################################################



	# archiving current configuration
	history[[g]] <- configs
	# determining which child was the fittest in terms of MDL scores
	best.where <- which(MDL[g,] == min(MDL[g,]))
	# if there were ties in the fittest child, randomly select one of them
	if (length(best.where) > 1) {
		best.where <- sample(best.where,1)
	}
	# storing the fittest child from the generation
	best.fit[[g]] <- configs[[best.where]]
	best.MDL[g] <- MDL[g,best.where]
   
	# count number of consecutive runs that did not see an improvement in MDL
	if(best.MDL[g] > min(best.MDL[1:(g-1)])) {
		run <- run + 1
	# reset the counter if there is an improvement in MDL
	} else {
		run <- 0
	}
	message("Number of runs without improvements on MDL: ", run)
 
	# calculating elapsed time to process each generation
	elapsed <- proc.time() - ptm
	message("It took ", as.numeric(elapsed[3]), " seconds to process generation ", g)
	message("Best fit for this generation was ")
	print(best.fit[[g]])
	message("with MDL of ")
	print(best.MDL[[g]])
	message(" ")
 	
	# If we have certain number of consecutive runs without an improvement in MDL, GA can be considered to have converged.
	if(run >= max.runs) {
		break
	}
}

# this stops the cluster that was created to handle parallel for loops
stopCluster(cluster)



##########################################################################################
##### End of Genetic Algorithm for Changepoint Detection                            ######
##########################################################################################



# Reporting total elapsed time
total.elapsed <- proc.time() - beginning.time
message("It took ", round(as.numeric(total.elapsed[3])/3600,2), " hours or ", round(as.numeric(total.elapsed[3])/60,2), " minutes total")

# determining which child was the fittest in terms of MDL scores
best.where <- which(best.MDL == min(best.MDL))
# if there were ties in the fittest child, randomly select one of them
if (length(best.where) > 1) {
	best.where <- sample(best.where, 1)
}

# Reporting the best fit
message("The fittest solution was the following")
print(best.fit[[best.where]])
message("MDL: ", round(best.MDL[best.where], 4))
message("From generation ", best.where)

# save the overall best fit
best.fit <- best.fit[[best.where]]

# save the estimated long-term trend parameter, alpha
trend.fitted <- fit.gev.changepoint(data,best.fit)$par[8]

# Exporting outputs
outputfilename1 <- paste("GA_Results_", tide_gauge_ID, ".RData", sep="")
save(best.fit, trend.fitted, data, file=outputfilename1)

message("All variables are exported to ", outputfilename1, " in the current working directory.")



##########################################################################################
##### End of File                                                                   ######
##########################################################################################
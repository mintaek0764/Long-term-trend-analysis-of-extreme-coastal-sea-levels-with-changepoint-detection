##########################################################################################
##########################################################################################
#####                                                                               ######
##### Compute BCa Bootstrap Confidence Intervals for Return Sea Levels and          ######
##### Long-term Linear Trend Parameter for the GESLA Data Set                       ######
#####                                                                               ######
##### Supplementary material to the "Long-term trend analysis of extreme coastal    ######
##### sea levels with changepoint detection"                                        ######
#####                                                                               ######
##### Last Updated on Monday, 7 December 2020                                       ######
#####                                                                               ######
##########################################################################################
##########################################################################################




##########################################################################################
##### Initial set up                                                                ######
##########################################################################################



## Start the clock
ptm <- proc.time()

# Taking the argument passed from the command line
args <- commandArgs(trailingOnly=TRUE)
# First (and the only) argument is GESLA tide gauge id, ranging from 1 to 300
index_to_check <- as.numeric(args[1])

## Loading required libraries
library(foreach)
library(doParallel)
library(optimr)

### Loading pre-defined custom R functions
source("Functions.R")

## Read the GESLA monthly maximum series
GESLA_dataset_name <- paste('MonthlyMax/SLstation', index_to_check, '.csv', sep='')
data_read <- read.csv(GESLA_dataset_name)
months <- data_read[,1]
data <- data_read[,2]

message("Processing for GESLA Station #", index_to_check)
message("")

# set the seed value to be used for random number generation
seed_value <- 20200821
# set the number of bootstrap samples to obtain
iteration <- 10
# set the size of the block for the moving block bootstrap
block.size <- 12
# set the confidence level for the confidence interval
confidence <- 0.95

### Loading the estimated extremal index
extremal_index <- scan("SavedData/GESLA_Extremal_Index.txt")
extremal_estimated <- extremal_index[[index_to_check]]

### Loading median sea level of the most recent 12 months
median_sea_level <- scan("SavedData/GESLA_Median_Sea_Level.txt")
median_station <- median_sea_level[index_to_check]

### Loading the GA estimated changepoints for the GESLA data set
load('SavedData/GA_Changepoints_GESLA.RData')
changepoint_estimated <- cp[[index_to_check]]



##########################################################################################
##### Begin Bootstrap CI for GESLA                                                  ######
##########################################################################################



## setup parallel backend to use many processors
cores <- detectCores()
cluster <- makeCluster(cores[1] - 1)
registerDoParallel(cluster)

message("Using ", cores[1], " cores.")
message("")

### Fitting the GEV model for original observations
MLE_GEV_cp <- GEV_return_exceedances_extremal(data, extremal = extremal_estimated, cp_station = changepoint_estimated, median_SL = median_station, last_month = max(months))
MLE_GEV_nocp <- GEV_return_exceedances_extremal(data, extremal = extremal_estimated, cp_station = numeric(0), median_SL = median_station, last_month = max(months))

### Performing pairs bootstrap

bootstrap_output <- foreach(i = 1:iteration,.packages="optimr") %dopar% {
	index <- moving.block(block.size, length(data),i)
	
	gev_fit_boot_cp <- GEV_return_exceedances_extremal(data, extremal = extremal_estimated, cp_station = changepoint_estimated, median_SL = median_station, last_month = max(months), boot_index = index)
	gev_fit_boot_nocp <- GEV_return_exceedances_extremal(data, extremal = extremal_estimated, cp_station = numeric(0), median_SL = median_station, last_month = max(months), boot_index = index)

	c(gev_fit_boot_cp$trend, gev_fit_boot_cp$return_exceedance_no_extremal, gev_fit_boot_cp$return_exceedance_extremal, gev_fit_boot_nocp$trend, gev_fit_boot_nocp$return_exceedance_no_extremal, gev_fit_boot_nocp$return_exceedance_extremal)
}

trends_cp_boot <- numeric()
re_cp_no_extremal_boot <- matrix(0, iteration, 4)
re_cp_extremal_boot <- matrix(0, iteration, 4)

trends_nocp_boot <- numeric()
re_nocp_no_extremal_boot <- matrix(0, iteration, 4)
re_nocp_extremal_boot <- matrix(0, iteration, 4)

for (i in 1:iteration) {
	# estimated linear trends with changepoints
	trends_cp_boot[i] <- bootstrap_output[[i]][1]
	# estimated return levels with changepoints
	re_cp_no_extremal_boot[i,] <- bootstrap_output[[i]][2:5]
	re_cp_extremal_boot[i,] <- bootstrap_output[[i]][6:9]
	# estimated linear trends without changepoints
	trends_nocp_boot[i] <- bootstrap_output[[i]][10]
	# estimated return levels without changepoints
	re_nocp_no_extremal_boot[i,] <- bootstrap_output[[i]][11:14]
	re_nocp_extremal_boot[i,] <- bootstrap_output[[i]][15:18]
}

### Calculating bias-correction constant, zBC

zBC <- numeric(18)
zBC[1] <- qnorm(sum(trends_cp_boot < MLE_GEV_cp$trend)/iteration, mean=0, sd=1)
for(i in 1:4) {
	zBC[1+i] <- qnorm(sum(re_cp_no_extremal_boot[,i] < MLE_GEV_cp$return_exceedance_no_extremal[i])/iteration, mean=0, sd=1)
	zBC[5+i] <- qnorm(sum(re_cp_extremal_boot[,i] < MLE_GEV_cp$return_exceedance_extremal[i])/iteration, mean=0, sd=1)
}
zBC[10] <- qnorm(sum(trends_nocp_boot < MLE_GEV_nocp$trend)/iteration,mean=0,sd=1)
for(i in 1:4) {
	zBC[10+i] <- qnorm(sum(re_nocp_no_extremal_boot[,i] < MLE_GEV_nocp$return_exceedance_no_extremal[i])/iteration,mean=0, sd=1)
	zBC[14+i] <- qnorm(sum(re_nocp_extremal_boot[,i] < MLE_GEV_nocp$return_exceedance_extremal[i])/iteration,mean=0, sd=1)
}


### Computing delete-1 Jackknife Return Level Samples
jackknife_output <- foreach(i = 1:length(data),.packages="optimr") %dopar% {
	index <- (1:length(data))[-i]
	gev_fit_jack_cp <- GEV_return_exceedances_extremal(data, extremal_estimated, cp_station = changepoint_estimated, median_SL = median_station, last_month = max(months), boot_index = index)
	gev_fit_jack_nocp <- GEV_return_exceedances_extremal(data, extremal_estimated, cp_station = numeric(0), median_SL = median_station, last_month = max(months), boot_index = index)

	c(gev_fit_jack_cp$trend, gev_fit_jack_cp$return_exceedance_no_extremal, gev_fit_jack_cp$return_exceedance_extremal, gev_fit_jack_nocp$trend, gev_fit_jack_nocp$return_exceedance_no_extremal, gev_fit_jack_nocp$return_exceedance_extremal)
}

trends_cp_jack <- numeric()
re_cp_no_extremal_jack <- matrix(0, length(data), 4)
re_cp_extremal_jack <- matrix(0, length(data), 4)

trends_nocp_jack <- numeric()
re_nocp_no_extremal_jack <- matrix(0, length(data), 4)
re_nocp_extremal_jack <- matrix(0, length(data), 4)

for (i in 1:length(data)) {
	# estimated linear trends with changepoints
	trends_cp_jack[i] <- jackknife_output[[i]][1]
	# estimated return levels with changepoints
	re_cp_no_extremal_jack[i,] <- jackknife_output[[i]][2:5]
	re_cp_extremal_jack[i,] <- jackknife_output[[i]][6:9]
	# estimated linear trends without changepoints
	trends_nocp_jack[i] <- jackknife_output[[i]][10]
	# estimated return levels without changepoints
	re_nocp_no_extremal_jack[i,] <- jackknife_output[[i]][11:14]
	re_nocp_extremal_jack[i,] <- jackknife_output[[i]][15:18]
}


### Calculating acceleration constant cA
cA <- numeric(18)
summation_trend_cp <- trends_cp_jack - mean(trends_cp_jack)
cA[1] <- (1/6)*(sum(summation_trend_cp^3)/(sum(summation_trend_cp^2)^(3/2)))
for(i in 1:4) {
	summation_re_cp_no_extremal <- re_cp_no_extremal_jack[,i] - mean(re_cp_no_extremal_jack[,i])
	cA[1+i] <- (1/6)*(sum(summation_re_cp_no_extremal^3)/(sum(summation_re_cp_no_extremal^2)^(3/2)))
	summation_re_cp_extremal <- re_cp_extremal_jack[,i] - mean(re_cp_extremal_jack[,i])
	cA[5+i] <- (1/6)*(sum(summation_re_cp_no_extremal^3)/(sum(summation_re_cp_no_extremal^2)^(3/2)))
}
summation_trend_nocp <- trends_nocp_jack - mean(trends_nocp_jack)
cA[10] <- (1/6)*(sum(summation_trend_nocp^3)/(sum(summation_trend_nocp^2)^(3/2)))
for(i in 1:4) {
	summation_re_nocp_no_extremal <- re_nocp_no_extremal_jack[,i] - mean(re_nocp_no_extremal_jack[,i])
	cA[10+i] <- (1/6)*(sum(summation_re_nocp_no_extremal^3)/(sum(summation_re_nocp_no_extremal^2)^(3/2)))
	summation_re_nocp_extremal <- re_nocp_extremal_jack[,i] - mean(re_nocp_extremal_jack[,i])
	cA[14+i] <- (1/6)*(sum(summation_re_nocp_no_extremal^3)/(sum(summation_re_nocp_no_extremal^2)^(3/2)))
}


### Calculating Quantiles for BCa
Z <- qnorm(1-(1-confidence)/2)

lower <- zBC + (zBC-Z)/(1-(cA*(zBC-Z)))
upper <- zBC + (zBC+Z)/(1-(cA*(zBC+Z)))

quantile_BCa <- pnorm(cbind(lower, upper))


### Computing Confidence Intervals with changepoints
CI_cp <- matrix(0, 10, 7)
CI_cp[,1] <- c("", "Trend","25-yr RL w/o theta", "50-yr RL w/o theta", "75-yr RL w/o theta", "100-yr RL w/o theta",  "25-yr RL w/ theta", "50-yr RL w/ theta", "75-yr RL w/ theta", "100-yr RL w/ theta")
CI_cp[1,] <- c("", "Percentile Lower", "Percentile Upper", "BCa Lower", "BCa Upper", "MLE", "Median")

# Assigning conventional percentile bootstrap CI on the left and BCa CI on the right
# For each CI, left: (alpha/2)%, middle: median, right: (1-alpha/2)%
CI_cp[2, 2:5] <- round(c(quantile(trends_cp_boot, (1-confidence)/2),quantile(trends_cp_boot, 1-(1-confidence)/2), quantile(trends_cp_boot, quantile_BCa[1,1]), quantile(trends_cp_boot, quantile_BCa[1,2])),4)
for(i in 1:4) {
	CI_cp[2+i,2:5] <- round(c(quantile(re_cp_no_extremal_boot[,i], (1-confidence)/2), quantile(re_cp_no_extremal_boot[,i], 1-(1-confidence)/2), quantile(re_cp_no_extremal_boot[,i], quantile_BCa[1+i,1]), quantile(re_cp_no_extremal_boot[,i], quantile_BCa[1+i,2])),4)
	CI_cp[6+i,2:5] <- round(c(quantile(re_cp_extremal_boot[,i], (1-confidence)/2),quantile(re_cp_extremal_boot[,i], 1-(1-confidence)/2), quantile(re_cp_extremal_boot[,i], quantile_BCa[5+i,1]), quantile(re_cp_extremal_boot[,i], quantile_BCa[5+i,2])),4)
}
CI_cp[,6] <- c("MLE", round(c(MLE_GEV_cp$trend,MLE_GEV_cp$return_exceedance_no_extremal, MLE_GEV_cp$return_exceedance_extremal), 4))
CI_cp[,7] <- c("Median", round(c(median(trends_cp_boot), median(re_cp_no_extremal_boot[,1]), median(re_cp_no_extremal_boot[,2]), median(re_cp_no_extremal_boot[,3]), median(re_cp_no_extremal_boot[,4]), median(re_cp_extremal_boot[,1]), median(re_cp_extremal_boot[,2]), median(re_cp_extremal_boot[,3]), median(re_cp_extremal_boot[,4])), 4))


### Computing Confidence Intervals without changepoints
CI_nocp <- matrix(0, 10, 7)
CI_nocp[,1] <- c("", "Trend", "25-yr RL w/o theta", "50-yr RL w/o theta", "75-yr RL w/o theta", "100-yr RL w/o theta", "25-yr RL w/ theta", "50-yr RL w/ theta", "75-yr RL w/ theta","100-yr RL w/ theta")
CI_nocp[1,] <- c("", "Percentile Lower", "Percentile Upper", "BCa Lower", "BCa Upper", "MLE", "Median")

# Assigning conventional percentile bootstrap CI on the left and BCa CI on the right
# For each CI, left: (alpha/2)%, middle: median, right: (1-alpha/2)%
CI_nocp[2,2:5] <- round(c(quantile(trends_nocp_boot, (1-confidence)/2), quantile(trends_nocp_boot, 1-(1-confidence)/2), quantile(trends_nocp_boot, quantile_BCa[10,1]), quantile(trends_nocp_boot, quantile_BCa[10,2])), 4)
for(i in 1:4) {
	CI_nocp[2+i,2:5] <- round(c(quantile(re_nocp_no_extremal_boot[,i], (1-confidence)/2), quantile(re_nocp_no_extremal_boot[,i], 1-(1-confidence)/2), quantile(re_nocp_no_extremal_boot[,i], quantile_BCa[10+i,1]), quantile(re_nocp_no_extremal_boot[,i], quantile_BCa[10+i,2])), 4)
	CI_nocp[6+i,2:5] <- round(c(quantile(re_nocp_extremal_boot[,i], (1-confidence)/2), quantile(re_nocp_extremal_boot[,i], 1-(1-confidence)/2), quantile(re_nocp_extremal_boot[,i], quantile_BCa[14+i,1]), quantile(re_nocp_extremal_boot[,i], quantile_BCa[14+i,2])), 4)
}
CI_nocp[,6] <- c("MLE",round(c(MLE_GEV_nocp$trend,MLE_GEV_nocp$return_exceedance_no_extremal,MLE_GEV_nocp$return_exceedance_extremal), 4))
	CI_nocp[,7] <- c("Median",round(c(median(trends_nocp_boot),median(re_nocp_no_extremal_boot[,1]), median(re_nocp_no_extremal_boot[,2]), median(re_nocp_no_extremal_boot[,3]), median(re_nocp_no_extremal_boot[,4]), median(re_nocp_extremal_boot[,1]), median(re_nocp_extremal_boot[,2]), median(re_nocp_extremal_boot[,3]), median(re_nocp_extremal_boot[,4])), 4))

outputfilename <- paste("Bootstrap_Results_", index_to_check, ".RData", sep="")
save(MLE_GEV_cp, MLE_GEV_nocp, bootstrap_output, zBC, jackknife_output, cA, quantile_BCa, CI_cp, CI_nocp, file=outputfilename)

message("All relevant variables are exported to ", outputfilename, " in the current working directory.")


### Stop the cluster when done
stopCluster(cluster)

### Stop the clock
time_elapsed <- proc.time() - ptm
message("It took ", round(as.numeric(time_elapsed[3])/3600, 4), " hours to run this code.")



##########################################################################################
##### End of File                                                                   ######
##########################################################################################
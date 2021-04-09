##########################################################################################
##########################################################################################
#####                                                                               ######
##### Main Analysis for the GESLA Data Set                                          ######
#####                                                                               ######
##### Supplementary material to the "Long-term trend analysis of extreme coastal    ######
##### sea levels with changepoint detection"                                        ######
#####                                                                               ######
##### Last Updated on Monday, 7 December 2020                                       ######
#####                                                                               ######
##########################################################################################
##########################################################################################



##########################################################################################
##### Preparing GESLA data set and changepoints for further analysis.               ######
##########################################################################################



### Loading necessary package
library(optimr)
library(foreach)
library(doParallel)

### Loading pre-defined custom R functions
source("Functions.R")

### reading necessary metadata
names <- readLines("SavedData/GESLA_Station_Names.txt")
coord <- read.table("SavedData/GESLA_Station_Coordinates.txt", head=TRUE)

### Loading the GA estimated changepoints for the GESLA data set
# changepoints are in terms of the observation index (1 through n) and not month index
load('SavedData/GA_Changepoints_GESLA.RData')

### Saving monthly maximum series for each station into one list
num_stn <- length(names)
SLmax <- list()
cp_month <- list()
trends_cp <- numeric(0)
for (i in 1:num_stn) {
	# Reading the GESLA data set
	file_name2 <- paste("MonthlyMax/SLstation", i, ".csv", sep="")
	SLmax[[i]] <- read.csv(file_name2)
	# Giving names to each column
	colnames(SLmax[[i]]) <- c("month", "max", "mean", "median")
	# Extracting the months index for the GA estimated changepoints
	cp_month[[i]] <- SLmax[[i]][,1][cp[[i]]]
}

### Loading the estimated extremal index
extremal_index <- scan("SavedData/GESLA_Extremal_Index.txt")

### Loading median sea level of the most recent 12 months
median_sea_level <- scan("SavedData/GESLA_Median_Sea_Level.txt")



##########################################################################################
##### Fitting GEV models with and without changepoints                              ######
##########################################################################################

### CAUTION: This takes a bit of time



# Create list objects to save the estimation results
gevfit_cp <- list()
gevfit_nocp <- list()

# Fit GEV models to each and every 300 selected GESLA stations
for (i in 1:num_stn) {
	data <- SLmax[[i]][,2]
	gevfit_cp[[i]] <- fit.gev.changepoint(data, cp[[i]])
	gevfit_nocp[[i]] <- fit.gev.changepoint(data, numeric(0))
}

# Extracte the estimated linear trends from the fitted GEV models
trends_cp <- numeric(num_stn)
trends_nocp <- numeric(num_stn)

for (i in 1:num_stn) {
	trends_cp[i] <- gevfit_cp[[i]]$par[12]
	trends_nocp[i] <- gevfit_nocp[[i]]$par[12]
}



##########################################################################################
##### Calculating the Nonstationary Return Exceedances from the Median Sea Level    ######
##########################################################################################

### CAUTION: This takes a quite a lot of time



## setup parallel backend to use many processors (if needed)
cores <- detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)

message("Using ", cores[1], " cores.")
message("")

# Calculate the estimated nonstationary return exceedances from the median sea levels for all 300 selected GESLA stations
# With changepoints considered
re_estimates_cp <- foreach(index_to_check = 1:num_stn,.packages="optimr") %dopar% {
  # read the data for the corresponding station
  read_data <- SLmax[[index_to_check]]
  data <- read_data[,2]
  months <- read_data[,1]
  last_month <- max(months)
  # extract the GA estimated changepoints for the station
  cp_station <- cp[[index_to_check]]
  # extract the estimated extremal index for the station
  extremal_station <- extremal_index[index_to_check]
  # extract the median sea level of the most recent 12 months for the station
  median_SL_station <- median_sea_level[index_to_check]
  # set what return periods to consider in return exceedance estimation
  years = c(25,50,75,100)
  
  # calculate the GEV estimated return exceedances from the median sea level
  # with and without extremal index
  GEV_return_exceedances_extremal(data, extremal_station, cp_station, median_SL_station, last_month, years)
}

# set up empty lists to store reassigned estimated return exceedances 
re_cp_extremal <- list()
re_cp_noextremal <- list()

# reassign estimated return exceedances for each case
for (index_to_check in 1:num_stn) {
  # Case 1: Consider both changepoints and extremal index
  re_cp_extremal[[index_to_check]] <- re_estimates_cp[[index_to_check]]$return_exceedance_extremal
  # Case 2: Consider changepoints but ignore extremal index
  re_cp_noextremal[[index_to_check]] <- re_estimates_cp[[index_to_check]]$return_exceedance_no_extremal
}

# Calculate the estimated nonstationary return exceedances from the median sea levels for all 300 selected GESLA stations
# Without changepoints considered

re_estimates_nocp <- foreach(index_to_check = 1:num_stn,.packages="optimr") %dopar% {
  # read the data for the corresponding station
  read_data <- SLmax[[index_to_check]]
  data <- read_data[,2]
  months <- read_data[,1]
  last_month <- max(months)
  # set no changepoints
  cp_station <- numeric(0)
  # extract the estimated extremal index for the station
  extremal_station <- extremal_index[index_to_check]
  # extract the median sea level of the most recent 12 months for the station
  median_SL_station <- median_sea_level[index_to_check]
  # set what return periods to consider in return exceedance estimation
  years = c(25,50,75,100)
  
  # calculate the GEV estimated return exceedances from the median sea level
  # with and without extremal index
  GEV_return_exceedances_extremal(data, extremal_station, cp_station, median_SL_station, last_month, years)
}

# set up empty lists to store reassigned estimated return exceedances 
re_nocp_extremal <- list()
re_nocp_noextremal <- list()

# reassign estimated return exceedances for each case
for (index_to_check in 1:num_stn) {
  # Case 3: Ignore changepoints but consider extremal index
  re_nocp_extremal[[index_to_check]] <- re_estimates_cp[[index_to_check]]$return_exceedance_extremal
  # Case 3: Ignore both changepoints and extremal index
  re_nocp_noextremal[[index_to_check]] <- re_estimates_cp[[index_to_check]]$return_exceedance_no_extremal
}

# this stops the cluster that was created to handle parallel for loops
stopCluster(cluster)



##########################################################################################
##### Make a quantile-quantile plot using Gumbel scale for a GESLA station          ######
##########################################################################################



# Specify what GESLA station to assess the GEV goodness-of-fit
index_to_check <- 63 # index for Fishguard, UK

read_data <- SLmax[[index_to_check]]
data <- read_data[,2]
cp_station <- cp[[index_to_check]]

# Fit the GEV model with changepoints
cp_gev <- fit.gev.changepoint(data,cp_station)
params <- cp_gev$par

# Fit the GEV model without changepoints
nocp_gev <- fit.gev.changepoint(data,numeric(0))
params <- nocp_gev$par

# QQ Plot for GEV model with changepoints
quantile_gumbel(data, cp[[index_to_check]])
# QQ Plot for GEV model without changepoints
quantile_gumbel(data, numeric(0))



##########################################################################################
##### End of File                                                                   ######
##########################################################################################
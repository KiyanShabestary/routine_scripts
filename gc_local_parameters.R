# 230301 KS
# Script to get local growth parameters of fit obtained from GrowthCurver package

require(growthcurver)
require(dplyr)

get_local_growth_params <- function(OD,time,trim_at_time){
  # Function that computes local growth parameters for a given column of a 
  # dataframe. To be used with apply.
  # Input:
  #       - column representing time in hours
  #       - each other column representing the timecourse of a given well, already blanked
  #       - trim_at_time value in hours
  # Output:
  #       - vector containing local (mu_lin) and exponential (mu_log) growth 
  #         estimates as well as doubling time (from mu_log).
  
  
  
  # Compute fit to get smooth version of the curve
  gc_fit <- SummarizeGrowth(data_t = time, 
                            data_n = OD,
                            t_trim = trim_at_time,
                            bg_correct = "none")
  
  # Extract the fit and binds it to the time column
  fit <- cbind(gc_fit$data$t, predict(gc_fit$model)) %>%
    as.data.frame() %>% dplyr::rename(time = V1, OD_fit = V2)
  rownames(fit)<-NULL
  
  # Compute log of the predicted fitted data
  fit$ln_OD_fit <- log(fit$OD_fit) # In R, log function is (logically) the natural logarithm
  
  # Calculate differences between OD and time at each timestep 
  # (note that in R: diff_i <- step_i+1 - step_i)
  local_rates = diff(as.matrix(fit)) %>% as.data.frame()
  
  local_rates$mu_lin <- local_rates$OD_fit/local_rates$time
  local_rates$mu_log <- local_rates$ln_OD_fit/local_rates$time
  
  max_mu_lin <- max(local_rates$mu_lin)
  max_mu_log <- max(local_rates$mu_log)
  
  min_doubling_time <- log(2)/max_mu_log
  
  c(max_mu_lin,max_mu_log,min_doubling_time)
}

# Download growthcurver example data
d <- growthdata

growth_params <- d %>% 
  select(-time) %>% 
  apply(MARGIN=2, FUN= get_local_growth_params, d$time, 20) %>%
  t() %>% as.data.frame() %>% dplyr::rename(max_mu_lin=V1,max_mu_log=V2,min_td=V3)
growth_params$well <- rownames(growth_params)
rownames(growth_params)<-NULL


# For debugging / example on one well

OD=d$A1
time=d$time
trim_at_time = 20

gc_fit <- SummarizeGrowth(data_t = time, 
                          data_n = OD,
                          t_trim = trim_at_time,
                          bg_correct = "none")

# Extract the fit and binds it to the time column
fit <- cbind(gc_fit$data$t, predict(gc_fit$model)) %>%
  as.data.frame() %>% dplyr::rename(time = V1, OD_fit = V2)
rownames(fit)<-NULL

# Compute log of the predicted fitted data
fit$ln_OD_fit <- log(fit$OD_fit) # In R, log function is (logically) the natural logarithm

# Calculate differences between OD and time at each timestep 
# (note that in R: diff_i <- step_i+1 - step_i)
local_rates = diff(as.matrix(fit)) %>% as.data.frame()

local_rates$mu_lin <- local_rates$OD_fit/local_rates$time
local_rates$mu_log <- local_rates$ln_OD_fit/local_rates$time

max_mu_lin <- max(local_rates$mu_lin)
max_mu_log <- max(local_rates$mu_log)

min_doubling_time <- log(2)/max_mu_log



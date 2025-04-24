library(tidyverse)
library(e1071)

####################################
#####       QUESTION 1        ######
####################################

#################
#  PART A   
#################
#initializing data
finch_data = read_csv("zebrafinches.csv")
further_data = finch_data$further
skew_further = skewness(further_data)
n=nrow(finch_data)
#Calculating t statistic
t_val = mean(further_data)/(sd(further_data)/sqrt(n))

#Caculating potential error
potential_err = pnorm(t_val) + (skew_further/sqrt(n))*((2*t_val^2+1)/6)*dnorm(t_val)

#################
#  PART B
#################
#Initializing
t_vals = seq(from = -10, to = 10, by = 0.01)
R = length(t_vals)
err_tibble = tibble(err=rep(NA,R))











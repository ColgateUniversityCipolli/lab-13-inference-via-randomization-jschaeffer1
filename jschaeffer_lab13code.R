library(tidyverse)
library(e1071)
library(boot)

####################################
#####       QUESTION 1        ######
####################################

#################
#  PART A   
#################
#initializing data
finch_data = read_csv("zebrafinches.csv")
further_data = finch_data$further
closer_data = finch_data$closer
diff_data = finch_data$diff
skew_further = skewness(further_data)
n=nrow(finch_data)
#Calculating t statistic
t_val = mean(further_data)/(sd(further_data)/sqrt(n))

#Caculating potential error
potential_err = (skew_further/sqrt(n))*((2*t_val^2+1)/6)*dnorm(t_val)

#################
#  PART B
#################
#Initializing
t_vals = seq(from = -10, to = 10, by = 0.01)
R = length(t_vals)
err_tibble = tibble(err=rep(NA,R))

#Looping over each t value
for (i in 1:R){
  t_val = t_vals[i]
  #Calculating potential error
  potential_err = (skew_further/sqrt(n))*((2*t_val^2+1)/6)*dnorm(t_val)
  
  #Adding error to tibble
  err_tibble$err[i] = potential_err
}
#view(err_tibble)

#Making plot of values
error_plot = ggplot(data=err_tibble, aes(t_vals, err)) +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept = 0) + 
  xlab("T Value") +
  ylab("Potential Error")


#################
#  PART C
#################
alpha = 0.05
t_val = qnorm(alpha)
n = ((skew_further/(6*0.1*alpha))* (2*t_val^2+1)*dnorm(t_val))^2



####################################
#####       QUESTION 2        ######
####################################

#################
#  PART A
#################

R <- 10000
resamples <- tibble(further = rep(NA, R),
                    closer = rep(NA, R),
                    diff = rep(NA, R))

for(i in 1:R){
  #Making a sample for further
  curr.resample <- sample(further_data,
                          size = length(further_data),
                          replace = T)
  
  resamples$further[i] <- mean(curr.resample)/(sd(further_data)/sqrt(length(further_data)))
  
  
  
  #Making a sample for closer
  curr.resample <- sample(closer_data,
                          size = length(closer_data),
                          replace = T)
  
  resamples$closer[i] <- mean(curr.resample)/(sd(closer_data)/sqrt(length(closer_data)))
  
  #Making a sample for difference
  curr.resample <- sample(diff_data,
                          size = length(diff_data),
                          replace = T)
  
  resamples$diff[i] <- mean(curr.resample)/(sd(diff_data)/sqrt(length(diff_data)))
}
resamples_shifted = resamples |>
  mutate(further = further-mean(further_data)/(sd(further_data)/sqrt(25))) |>
  mutate(closer = closer-mean(closer_data)/(sd(closer_data)/sqrt(25))) |>
  mutate(diff = diff-mean(diff_data)/(sd(diff_data)/sqrt(25)))

#view(resamples)
#view(resamples_shifted)

#################
#  PART B
#################
pval_further = (mean(resamples_shifted$further <= mean(resamples$further)))/R
pval_closer = (mean(resamples_shifted$closer >= mean(resamples$diff)))/R
pval_diff = (mean(resamples_shifted$diff >= mean(resamples$diff)))/R

pval_further
pval_closer
pval_diff

#################
#  PART C
#################
resample_further_5th = quantile(resamples_shifted$further, c(0.05, 0.95))
resample_closer_5th = quantile(resamples_shifted$closer, c(0.05, 0.95))
resample_diff_5th = quantile(resamples_shifted$diff, c(0.05, 0.95))

resample_further_5th
resample_closer_5th
resample_diff_5th

#################
#  PART D
#################

# Confidence Interval
resample_further_ci = quantile(resamples$further, c(0.025, 0.975))
resample_closer_ci = quantile(resamples$closer, c(0.025, 0.975))
resample_diff_ci = quantile(resamples$diff, c(0.025, 0.975))

resample_further_ci
resample_closer_ci
resample_diff_ci




####################################
#####       QUESTION 3        ######
####################################

#################
#  PART A
#################
closer_data = finch_data$closer
diff_data = finch_data$diff

R <- 10000
rand <- tibble(further = rep(NA, R),
               closer = rep(NA, R),
               diff = rep(NA, R))

#Randomization for further
# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- further_data - 0
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift *
    sample(x = c(-1, 1),
           size = length(x.shift),
           replace = T)
  
  rand$further[i] <- mean(curr.rand)
}

#Randomization for closer
x.shift <- closer_data - 0
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift *
    sample(x = c(-1, 1),
           size = length(x.shift),
           replace = T)
  
  rand$closer[i] <- mean(curr.rand)
}


#Randomization for difference
x.shift <- diff_data - 0
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift *
    sample(x = c(-1, 1),
           size = length(x.shift),
           replace = T)
  
  rand$diff[i] <- mean(curr.rand)
}
#view(rand)



#################
#  PART B
#################
# p-value for further
(delta <- abs(mean(further_data) - 0))
(low <- 0 - delta) # mirror
(high<- 0 + delta)   # xbar

pval_further = mean(rand$further <= low) + mean(rand$further >= high)

# p-value for closer
(delta <- abs(mean(closer_data) - 0))
(low <- 0 - delta) # mirror
(high<- 0 + delta)   # xbar

pval_closer = mean(rand$closer <= low) + mean(rand$closer >= high)

# p-value for further
(delta <- abs(mean(diff_data) - 0))
(low <- 0 - delta) # mirror
(high<- 0 + delta)   # xbar

pval_diff = mean(rand$diff <= low) + mean(rand$diff >= high)


#################
#  PART C
#################
R <- 1000
mu0.iterate <- 0.001
starting.point <- mean(further_data)

#Calculating lower value
mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- further_data - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  (delta <- abs(mean(further_data) - mu.lower))
  (low <- mu.lower - delta) # mirror
  (high<- mu.lower + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

#Calculating upper value
mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- further_data - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  (delta <- abs(mean(further_data) - mu.upper))
  (low <- mu.upper - delta) # mirror
  (high<- mu.upper + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}

further_ci = c(mu.lower, mu.upper)


###################
#    CLOSER
###################
starting.point <- mean(closer_data)

#Calculating lower value
mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- closer_data - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  (delta <- abs(mean(closer_data) - mu.lower))
  (low <- mu.lower - delta) # mirror
  (high<- mu.lower + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

#Calculating upper value
mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- closer_data - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  (delta <- abs(mean(closer_data) - mu.upper))
  (low <- mu.upper - delta) # mirror
  (high<- mu.upper + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}

closer_ci = c(mu.lower, mu.upper)



###################
#    DIFFERENCE
###################
starting.point <- mean(diff_data)

#Calculating lower value
mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- diff_data - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  (delta <- abs(mean(diff_data) - mu.lower))
  (low <- mu.lower - delta) # mirror
  (high<- mu.lower + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

#Calculating upper value
mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- diff_data - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  (delta <- abs(mean(diff_data) - mu.upper))
  (low <- mu.upper - delta) # mirror
  (high<- mu.upper + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}

diff_ci = c(mu.lower, mu.upper)













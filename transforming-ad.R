library(mlbench)   #Version 2.1.3
library(Cubist)    #Version 0.4.2.1
library(tseriesChaos) #Version 0.1.13.1
library(TTR)       #Version 0.24.3
library(ifultools) #Version 2.0.26
library(wmtsa)     #Version 2.0.3
library(DescTools) #Version 0.99.49

library(scatterplot3d) #Version 0.3.44
library(forecast)   #Version 8.21
library(gtools)     #Version 3.9.4 
library(infotheo)   #Version 1.2.0.1

library(rEDM)       #Version 1.2.3
library(seasonal)#Version 1.9.0
library(ape)     #Version 5.7.1
library(rgee)    #Version 1.1.5
library(ggplot2) #Version 3.4.2
library(tidyr)   #Version 1.3.0
library(ggplot2) #Version 3.4.2
library(reticulate) #Version 1.28
library(earlywarnings) #Version 1.1.29



#######################################
# function 1: data-arranging function #
#######################################
# First, transfer the week to date "%Y-%m-%d", because the original CSV data didn't report date but week
#second, order incidence in date (from oldest date to newest date)

reorganise<- function (week, year, incidences){
  datasets=as.data.frame(cbind(week,year,incidences)) # merge
  colnames(datasets)=c("week", "year","incidences")
  datasets=datasets[datasets$week<=52,] #a year has 52 weeks, thus remove >53 week, otherwise as.Date function will give errors
  datasets$time<- as.Date(paste(datasets$year, datasets$week,1, sep="-"), "%Y-%U-%u")
  datasets$time<- as.Date(datasets$time, "%Y-%m-%d") #transfer time to standard format
  library(tidyverse)
  library(lubridate)
  library(dplyr)
  datasets=datasets[!is.na(datasets$time), ] #remove NA date
  datasets=datasets[!duplicated(datasets$time), ] #remove duplicate time
  datasets=datasets[order(as.Date(datasets$time, format="%Y-%m-%d")),] # sort incidences in date (from oldest date to newest date)
  
  return(as.data.frame(datasets))
}

###################################################################################################
#function 2: compute a Earlier warming signal (e.g. the absolute value of the maximum eigenvalue) #
###################################################################################################
abs_max_eigenvalue <- function (dataset=datasets_1, window_size=50, steps=1){
  
  library(rEDM)   #package version: 1.2.3
  library("sapa") #package version, 2.0.3
  raw_data=dataset$incidences
  # setting values
  tau <- 1                   # time lag is 1, i.e. x(T+1)=f(x(t))
  theta <- seq(0,2.5,by=0.5) # setting theta
  step_size <- 1             #moving forward 1 time-step
  window_indices <- seq(window_size, NROW(dataset), step_size)
  matrix_result <- as.data.frame( matrix(NaN, nrow = length(window_indices), ncol = 3) )
  index <- 0
  #
  for(j in window_indices) {
    index <- index + 1
    rolling_window <- raw_data[(j-window_size+1):j]   #select data for each rolling window
    norm_rolling_window <- (rolling_window - mean(rolling_window, na.rm=TRUE))/sd(rolling_window, na.rm=TRUE)
    #use the simplex projection to find best E (embedding dimension) for each window
    sim_r <- simplex(norm_rolling_window,lib=c(1,floor(length(norm_rolling_window)/2)),pred=c(floor(length(norm_rolling_window)/2)+1,length(norm_rolling_window)),E=c(1:11))
    E <-sim_r[which.min(sim_r$mae),"E"][1] 
    # find best theta
    smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta, silent=TRUE)
    best <- order(as.numeric(smap$mae))[1] # find location of best theta having lowest prediction error i.e. lowest mae 
    theta_best <- smap[best,]$theta        # find best theta having lowest prediction error i.e. lowest mae 
    # calculate interaction strength with best theta above and with the best E above
    smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta_best, silent=TRUE, save_smap_coefficients=TRUE)
    smap_co <- smap$smap_coefficients[[1]]
    # calculate maximum of eigenvalue
    matrix_eigen <- matrix(NA, nrow = NROW(smap_co), ncol = 3)
    for(k in 1:NROW(smap_co)){
      if(!is.na(smap_co[k,1]))
      {
        M <- rbind(as.numeric(smap_co[k, 1:E]), cbind(diag(E - 1), rep(0, E - 1)))
        M_eigen <- eigen(M)$values
        lambda1 <- M_eigen[order(abs(M_eigen))[E]]
        
        matrix_eigen[k,3] <- abs(lambda1) #the absolute value of the maximum eigenvalue
      }
    }
    # save results
    matrix_result[index,1] <- j+steps-1
    matrix_result[index,2] <- as.character( as.Date(dataset$time[c(j)],  "%Y-%m-%d") ) # add date back in to dataset
    matrix_result[index,3] <- mean(matrix_eigen[,3],na.rm=TRUE)
  }
  colnames(matrix_result)<-c("number","time", "abs.lambda")
  matrix_result$time<-as.Date(matrix_result$time, "%Y-%m-%d")
  matrix_result$number<-as.numeric(matrix_result$number)
  matrix_result$abs.lambda<-as.numeric(matrix_result$abs.lambda)
  return(as.data.frame(matrix_result))
  
}






# step 1: set working directory
setwd('C:/Users/QINGHUA ZHAO/Downloads/EWs-main (1)/EWs-main/EWs2023Sep24/EWs2023Sep24') # change it to your working routine

# step2: load data
original_data=read.csv("colombia.disease.Pertussis.csv") 
head(original_data)

#step3: reorganise the data. Here only 3 columns in original_data will be used to reorganise, i.e. "week", "year", and "Incidence" column 
datasets_1=reorganise(original_data$EpiWeek,original_data$EpiYear,original_data$Incidence )
head(datasets_1)

#function 4: compute a Earlier warming signal (e.g. the absolute value of the maximum eigenvalue) 
eigenvalues=abs_max_eigenvalue(dataset=datasets_1, window_size=50, steps=1) # 5o is the width of rolling window

# plot
ggplot(eigenvalues,aes(x = time, y = abs.lambda)) + 
  geom_line() + theme_bw() + xlab("Time") + ylab("Max(Eigen)") # larger than 1 means disease outbreak


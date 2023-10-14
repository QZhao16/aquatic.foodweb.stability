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




#### function: compute a Earlier warming signal (e.g. variance) 
the_variance <- function (dataset=original_data, window_size=50, steps=1){ #the width of rolling window is 50.
  inciddences=original_data$Incidence                      # only use the Incidence
  window_indices <- seq(window_size, NROW(dataset), steps) #cutting data into rolling windows
  #
  variances <-NULL
  for(j in window_indices) {
    variances <- c(variances, var(inciddences[(j-window_size+1):j]) )   #variance for each rolling window
  }
  return(variances)
}



# step 1: set working directory
setwd('C:/Users/QINGHUA ZHAO/Downloads/EWs-main (1)/EWs-main/EWs2023Sep24/EWs2023Sep24') # change it to your working routine

# step2: load data
original_data=read.csv("colombia.disease.Pertussis.csv") 
head(original_data)

#step 3, compute a Earlier warming signal (e.g. variance) 
results=the_variance(dataset=original_data, window_size=50, steps=1)

#plot results
plot(results, type="o", ylab="variance in each time window", col="blue")


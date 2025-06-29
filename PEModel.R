#Model of premature hive exiting behavior in honey bee colonies using parameters and 
#equations from perry et al. 2015 and Khourey et al. 2011

library(dde)
library(deSolve)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)

#Create delay differential equation function
#All parameters are inside function, except pe
bee_dde <- function(t, y, parms){
  
  #Brood break implementation
  #Either run this if statement, or assign L to 2000
  # if (t < 50){
  #   L = 2000
  # } else if (t< 70){
  #   L = 0
  # } else{
  L = 2000
  # }
  # 
  
  # if (t < 60){
  #    pe = parms[1]
  #  } else{
  #    pe = parms[2]
  #  }
  
  alphamin = 0.25
  sigma = 0.75
  v = 5000
  alphamax = 0.25
  phi = 1/9
  gammaA = 0.007
  gammaB = 0.018
  #cf is the added feeding parameter
  #cf = 0
  ct = 0.033
  b = 500
  a0 = 5
  pe = parms[1]
  mr = ((15*pe^2) + 1)
  delay = 12
  y1 = y[[1L]]
  y2 = y[[2L]]
  y3 = y[[3L]]
  y4 = y[[4L]]
  
  #tau = t - delay
  if(t < 12){
    y_lag = 0
  }
  
  else{
    y_lag = deSolve::lagvalue(t - 12)[2]
  }
  
  R = alphamin + (alphamax*(b^2/(b^2+y1^2))) - (sigma*(y4)/(y4+y3))
  a = 1/R
  
  if(a <= 10){
    Na = -0.02*(a-10)^2+3.5
  }
  else{
    Na = -0.015*(a-10)^2 +3.5
  }
  
  if(a <= 40/3){
    Ta = 0.06 * (a-5) +0.5
  }
  else{
    Ta = 1
  }
  
  Ma = (((a-a0)^4) + 3) / ((4.94 + 0.08*a) * (a-a0)^4)
  
  dfdt <- (ct*Na*y4) - (gammaA * (y4 +y3)) - (gammaB*y2) # + cf
  dBdt <- (((L*y1^2) * y3) / ((b^2 + y1^2) * (y3+v))) - (phi*y2)
  dHdt <- (phi*y_lag) - (R * y3) - pe*y3
  dFdt <- (Ta * R * y3) - (mr * Ma * y4)
  
  list(c(dfdt,dBdt,dHdt,dFdt))
}



#Label timesteps for simulation to run over
tt = seq(0,150)

#Label initial values for model
y0 = c(y1 = 1000, y2 = 0, y3 = 16000, y4 = 8000)

#Create a premature hive exiting parameter, values used in manuscript: 0.007, 0.04, 0.08, and 0.125
parms = c(pe = 0.007)


#Run model and time it
yout = dede(y = y0, times = tt, func = bee_dde,
                         parms = parms)

#Plot output, pe values used in manuscript: 0.007, 0.04, 0.08, and 0.125
matplot(yout[,1], yout[,-1], type = "l", lwd = 2, lty = 1,
        main = "Pe 0.007", xlab = "Days",
        ylab = "Stored food (g) or number of bees",
        ylim = c(0,25000))

#Prepare data to plot as a dataframe in ggplot
colnames(yout) = c("time","food", "brood", "hive", "foragers")

df <- as.data.frame(yout)
df_long <- pivot_longer(data = df, cols = !time, names_to = "cat")

ggplot()+
  geom_smooth(data = df_long, aes(x = time, y = value, col = cat), se = F)+
  theme(legend.position = "top")+
  labs(col = "")

##### PE VS PE PLOT ######

#Choose range of self removal rates
removal_rates <- seq(from = 0, to = 0.125, length.out = 10) 
removal_rates <- c(removal_rates, 0.119, 0.13)

#Create a function to repeatedly run simulations with different pe params using lapply later
run_many_sims <- function(parms){
  yout <- dede(y = c(y1 = 1000, y2 = 0, y3 = 16000, y4 = 8000), 
               times = seq(0,300), 
               func = bee_dde,
               parms = parms)
  colnames(yout) = c("time","food", "brood", "hive", "foragers")
  yout = data.frame(yout)
  yout$pe = parms
  yout
}

#Create function to get foraging age of list of simulations
get_foraging_age <- function(df){
  ageForaging <- 0.25 +(0.25*(500^2/(500^2+df$food^2))) - (0.75*(df$foragers/(df$foragers + df$hive)))
  ageForaging <- 1/ageForaging
  ages <- data.frame(time = df$time, age = ageForaging, pe = df$pe)
}

simResults <- lapply(X = removal_rates, FUN = run_many_sims)
ages <- lapply(X = simResults, FUN = get_foraging_age)

allAges <- data.frame(mapply(c, ages[[1]], ages[[2]], ages[[3]], ages[[4]], ages[[5]],
                             ages[[6]], ages[[7]], ages[[8]], ages[[9]], ages[[10]], ages[[11]], ages[[12]]))

allAges$pe <- round(allAges$pe, 4)
allAges$pe = factor(allAges$pe)

ggplot(allAges, aes(x = time, y = age, col = pe))+
  geom_line()+
  ylab("Age of First Foraging")+
  xlab("Time")

#Mite PE Trendline
mites <- read_csv("./MiteSRtrendline.csv")
colnames(mites) <- c("MiteCount", "SRCount", "TotalBees", "Rate")
ggplot(mites, aes(x = MiteCount, y = Rate))+
  geom_smooth(method = "lm", se = F)+
  geom_point()+
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15))+
  ylim(0, 0.15)


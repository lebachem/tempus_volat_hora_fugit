# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ------                         Survey Paper                           ------ #
# ------                        (Fitting a REM)                         ------ #
# ------                  Dataset 1: International Arms Trade           ------ #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
rm(list=ls())
load("load_data.RData")
# Download the package online: http://www.sn.ethz.ch/research/goldfish.html
library(goldfish)
library(data.table)

source("../help_functions.R")


beginning = 2016
ending = 2017
# Which countries are present in all 2 years? 
# Save the subgraph of those networks 
tmp_networks = list()
countries_ind =which(rowSums(EX[,beginning:ending -1949]) == 2)

tmp_networks[[1]] = amk[[beginning - 1949]][countries_ind, countries_ind]
tmp_networks[[2]] = amk[[ending - 1949]][countries_ind, countries_ind]

# Breslow Approximation----

# Save the events 
event_stream = which((tmp_networks[[2]]>0 ), arr.ind=T)

# Set the time of all events to 1 (equal to the Breslow approximation)
set.seed(1236)
event_stream = data.frame("time" = runif(nrow(event_stream)),
                          "sender" =as.character(names(countries_ind)[event_stream[,1]]), 
                          "receiver"= as.character(names(countries_ind)[event_stream[,2]]), 
                          "increment" = rep(1,nrow(event_stream)))
event_stream$receiver = as.character(event_stream$receiver)
event_stream$sender = as.character(event_stream$sender)
event_stream = event_stream[order(event_stream$time),]

# Add exogenous information 
polity_countries = polity[countries_ind,beginning - 1949]
countries <- data.frame(label = as.character(names(countries_ind)))
countries$gdp =  log(real_gdp[countries_ind,beginning - 1949])
countries$polity = polity_countries
countries$label = as.character(countries$label) 

# Define the network 
armsNetwork <- defineNetwork(nodes = countries,directed = TRUE, matrix = apply(tmp_networks[[1]]>0, 2, as.numeric))
# Link the network with the observed events 
armsNetwork <- linkEvents(x =armsNetwork , changeEvent = event_stream,
                          nodes = countries)
# Declare the network to be the target in the estimation procedure
armsDependent <- defineDependentEvents(events =  event_stream, 
                                       nodes = countries,
                                       defaultNetwork = armsNetwork)

# Estimation 
rem_arms_trade = estimate(armsDependent ~ inertia + outdeg_sender+  recip +
                            indeg +trans+ ego(countries$gdp) +
                            alter(countries$gdp) + diff(countries$polity),
               modelType = "REM")

rem_arms_trade_data = estimate(armsDependent ~ inertia + outdeg_sender+  recip + indeg +trans+ ego(countries$gdp) + alter(countries$gdp) + diff(countries$polity),
                          modelType = "REM",returnStatisticsOnly = T)

# Kalbfleisch-Prentence Approximation----


B = 100
p =8

coef = array(dim = c(p,B))
std_err = array(dim = c(p,B))
loglik = array(dim = c(B))
model = list()

for(b in 1:B){
  event_stream = which((tmp_networks[[2]]>0 ), arr.ind=T)
  set.seed(b)
  event_stream = data.frame("time" =runif(nrow(event_stream)),
                            "sender" =as.character(names(countries_ind)[event_stream[,1]]), 
                            "receiver"= as.character(names(countries_ind)[event_stream[,2]]), 
                            "increment" = 1)
  event_stream$receiver = as.character(event_stream$receiver)
  event_stream$sender = as.character(event_stream$sender)
  event_stream = event_stream[order(event_stream$time),]
  
  countries <- data.frame(label = as.character(names(countries_ind)))
  countries$gdp =  log(real_gdp[countries_ind,beginning - 1949])
  countries$polity = polity_countries
  countries$label = as.character(countries$label) 
  
  armsNetwork <- defineNetwork(nodes = countries,directed = TRUE, 
                               matrix = apply(tmp_networks[[1]]>0, 2, as.numeric))
 
  armsNetwork <- linkEvents(x =armsNetwork , changeEvent = event_stream, nodes = countries)
  armsDependent <- defineDependentEvents(events =  event_stream, 
                                         nodes = countries,
                                         defaultNetwork = armsNetwork)
  res = estimate(armsDependent ~ inertia + outdeg_sender+  recip + indeg +
                   trans+ ego(countries$gdp) + alter(countries$gdp) + diff(countries$polity),
                 modelType = "REM")
  
  coef[,b] = res$parameters
  std_err[,b] = res$standard.errors
  loglik[b] = res$log.likelihood
  model[[b]] = res
}
coefs = apply(coef,MARGIN = 1,FUN  = mean)
mean_t_statistics = apply(rbindlist(
  lapply(model, 
         FUN = function(x){data.frame(t(x$parameters/x$standard.errors))})),
  FUN =  mean,MARGIN = 2) 
1-pt(df = 8,mean_t_statistics) # p-values of the mean t-statistics

stadt_errs = apply(std_err,MARGIN = 1,FUN  = mean)

aics = mean(-2*loglik + 2*7)
save(aics,file = "aics_rem1.RData")

loglik_mean = mean(loglik)
save(loglik_mean,file = "loglik_mean_rem1.RData")
# Save the results ----
save.image("rem_1.RData")

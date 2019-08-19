# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ------                         Survey Paper                           ------ #
# ------                        (Fitting a REM)                         ------ #
# ------ Data Set 2: European Research Institution Email Correspondence ------ #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
 
# Preliminaries ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
library(goldfish)
source("../help_functions.R")

# Read the data -----
data_tmp = read.table("email-Eu-core-temporal-Dept3.txt")
names(data_tmp) = c("from", "to", "time")

# Use only the first two years of the correspondence data
# 31540000 seconds is approx a year
data_tmp_trunc = data_tmp[data_tmp$time< 31540000*2,]
# Delete possible self-loops
data_tmp_trunc = data_tmp_trunc[!(data_tmp_trunc$from == data_tmp_trunc$to),]
# Delete all group emails (emails with more than one receiver)
allowed_times = summary(factor(data_tmp_trunc$time),maxsum =8831) == 1
allowed_times = as.numeric(levels(factor(data_tmp_trunc$time))[allowed_times])
data_tmp_trunc = data_tmp_trunc[data_tmp_trunc$time %in% allowed_times, ]

# Save the actors of the relational data in a vector 
actors = unique(c(data_tmp_trunc$from, data_tmp_trunc$to))
# Change the actor labels to be 1, ..., n_actors
data_tmp_trunc$from_num = match(data_tmp_trunc$from, actors)
data_tmp_trunc$to_num = match(data_tmp_trunc$to, actors)

# The first year will be used as the default network (conditioned on)
# 31540000 seconds is approx a year
data_tmp_trunc1 = data_tmp_trunc[data_tmp_trunc$time <31540000,]

# The second year of observations is modeled with the REM
data_tmp_trunc2 =data_tmp_trunc[data_tmp_trunc$time >= 31540000,]

# Estimate the REM -----

event_stream = data.frame("time" = data_tmp_trunc2$time,
                          "sender" =data_tmp_trunc2$from_num, 
                          "receiver"=data_tmp_trunc2$to_num, 
                          "increment" = rep(1,nrow(data_tmp_trunc2)))
event_stream$receiver = as.character(event_stream$receiver)
event_stream$sender = as.character(event_stream$sender)
event_stream = event_stream[order(event_stream$time),]

build_stream = data.frame("time" = data_tmp_trunc1$time,
                          "sender" =data_tmp_trunc1$from_num, 
                          "receiver"=data_tmp_trunc1$to_num, 
                          "increment" = rep(1,nrow(data_tmp_trunc1)))
build_stream$receiver = as.character(build_stream$receiver)
build_stream$sender = as.character(build_stream$sender)
build_stream = build_stream[order(build_stream$time),]


actorsset <- data.frame(label = as.character(1:length(actors)))
actorsset$label = as.character(actorsset$label) 

network_1 = matrix(data = 0,nrow = length(actors), ncol = length(actors))
data_tmp_trunc1$dyad = factor(paste(data_tmp_trunc1$from_num,
                                    data_tmp_trunc1$to_num, sep =  "_"))
from_tmp = unlist(lapply(strsplit(as.character(levels(data_tmp_trunc1$dyad)),split = "_"),
                         FUN = function(x){
                           return(as.numeric(x[[1]]))
                         }))
to_tmp = unlist(lapply(strsplit(as.character(levels(data_tmp_trunc1$dyad)),split = "_"),
                       FUN = function(x){
                         return(as.numeric(x[[2]]))
                       }))

network_1[cbind(from_tmp, to_tmp)] = table(data_tmp_trunc1$dyad,exclude = F)

emailNetwork <- defineNetwork(nodes = actorsset,directed = TRUE,matrix = network_1)
emailNetwork <- linkEvents(x =emailNetwork , changeEvent = event_stream, 
                           nodes = actorsset)

emailDependent <- defineDependentEvents(events =  event_stream,
                                        nodes = actorsset,defaultNetwork = emailNetwork)
emailBuild <- defineDependentEvents(events =  build_stream,
                                        nodes = actorsset,defaultNetwork = emailNetwork)

rem_email = estimate(emailDependent ~ inertia + outdeg_sender+  recip + indeg +trans,
                     modelType = "REM")

rem_email_data = estimate(emailDependent ~ inertia + outdeg_sender+  recip + indeg +trans,
                     modelType = "REM",returnStatisticsOnly = T )

# Save GOF-Stats----
par(mfrow = c(1,1))

next_tie = c()
next_sender = c()
next_receiver = c()

for(i in 1:2536){
  tmp_mat = exp(rem_email_data$dep[[i]]$change.contributions[,,1]*rem_email$parameters[1] + 
                  rem_email_data$dep[[i]]$change.contributions[,,2]*rem_email$parameters[2] + 
                  rem_email_data$dep[[i]]$change.contributions[,,3]*rem_email$parameters[3]+ 
                  rem_email_data$dep[[i]]$change.contributions[,,4]*rem_email$parameters[4] +
                  rem_email_data$dep[[i]]$change.contributions[,,5]*rem_email$parameters[5])
  
  tmp_sender = tmp_mat[ as.numeric(event_stream[i+1,2]),]
  tmp_receiver =  tmp_mat[,as.numeric(event_stream[i+1,3])]

  next_tie[i] = sum(tmp_mat>tmp_mat[as.numeric(event_stream[i+1,c(2)]),as.numeric(event_stream[i+1,c(3)])])
  next_sender[i] = sum(tmp_sender>tmp_sender[as.numeric(event_stream[i+1,c(3)])])
  next_receiver[i] = sum(tmp_receiver>tmp_receiver[as.numeric(event_stream[i+1,c(2)])])
}

recall_tie = c()
recall_sender = c()
recall_receiver = c()

for(i in 1:400){
  recall_tie[i] = sum(next_tie<=i)/length(next_tie)
  recall_sender[i] = sum(next_sender<=i)/length(next_sender)
  recall_receiver[i] = sum(next_receiver<=i)/length(next_receiver)
  
}

plot_data = data.frame("cutoff" = 1:400, 
                       "recall_tie" = recall_tie, 
                       "recall_sender" = recall_sender, 
                       "recall_receiver" = recall_receiver)
line_df = data.frame(x = c(0,400),y = c(0,400*1/(89*88)))
pdf("recall_ties.pdf",width =8,height =8)
ggplot(data = plot_data ,mapping = aes(x = cutoff,y = recall_tie)) +
  geom_step(alpha = 0.8, size=1) +
  ylab("Recall")+
  xlab("Cutoff") +
  theme_pubr() +
  geom_line(data = line_df,aes(x,y), lty = 2 ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30)) 
dev.off()

line_df = data.frame(x = c(0,89),y = c(0,89*1/(89)))

pdf("recall_sender.pdf",width = 8,height =8)
ggplot(data = plot_data ,mapping = aes(x = cutoff,y = recall_sender)) +
  geom_step(alpha = 0.8, size=1) +
  ylab("Recall")+
  xlab("Cutoff") +
  theme_pubr() +
  geom_line(data = line_df,aes(x,y), lty = 2 ) +
  xlim(c(0,89)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30)) 
dev.off()

pdf("recall_receiver.pdf",width = 8,height =8)
ggplot(data = plot_data ,mapping = aes(x = cutoff,y = recall_receiver)) +
  geom_step(alpha = 0.8, size=1) +
  ylab("Recall")+
  xlab("Cutoff") +
  theme_pubr() +
  geom_line(data = line_df,aes(x,y), lty = 2 ) +
  xlim(c(0,89)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30)) 
dev.off()

# Save the results ----
save.image("rem_2.RData")

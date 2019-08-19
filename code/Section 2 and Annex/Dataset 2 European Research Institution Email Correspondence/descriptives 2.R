# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ------                         Survey Paper                           ------ #
# ------                        (Descriptives)                          ------ #
# ------ Data Set 2: European Research Institution Email Correspondence ------ #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# Preliminaries ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
rm(list=ls())

library(statnet)
library(texreg)
library(ggpubr)
library(reshape2)
# set seed
set.seed(123)

# Read the data -----
data_tmp = read.table("email-Eu-core-temporal-Dept3.txt")
names(data_tmp) = c("from", "to", "time")

# Use only the first tweo years of the email traffic
data_tmp_trunc = data_tmp[data_tmp$time< 31540000*2,]
data_tmp_trunc = data_tmp_trunc[!(data_tmp_trunc$from == data_tmp_trunc$to),]
# Delete all group emails 
allowed_times = summary(factor(data_tmp_trunc$time),maxsum =8831) == 1
allowed_times = as.numeric(levels(factor(data_tmp_trunc$time))[allowed_times])
data_tmp_trunc = data_tmp_trunc[data_tmp_trunc$time %in% allowed_times, ]

# Save the actors in a vector 
actors = unique(c(data_tmp_trunc$from, data_tmp_trunc$to))

# Change the actor labels to be 1, ..., n_actors
data_tmp_trunc$from_num = match(data_tmp_trunc$from, actors)
data_tmp_trunc$to_num = match(data_tmp_trunc$to, actors)

# The first aggregated network is the first year (yearly aggregation)
data_tmp_trunc1 = data_tmp_trunc[data_tmp_trunc$time <31540000,]
network_1 = matrix(data = 0,nrow = length(actors), ncol = length(actors))
network_1[cbind(data_tmp_trunc1$from_num, data_tmp_trunc1$to_num)] = 1
network_1 = network(network_1)

# The second aggregated network is the second half of observations (second year)
data_tmp_trunc2 =data_tmp_trunc[data_tmp_trunc$time >= 31540000,]
network_2 = matrix(data = 0,nrow = length(actors), ncol = length(actors))
network_2[cbind(data_tmp_trunc2$from_num, data_tmp_trunc2$to_num)] = 1
network_2 = network(network_2)
tmp_networks = list(network_1, network_2)

# Save those networks in an edge list for the use in gephi
matrix_1 = which(as.matrix(network_1)>0,arr.ind = TRUE)
matrix_1 = data.frame("Timestamp1" = 0,"Timestamp2" = 1,"Source" = matrix_1[,1],"Target"= matrix_1[,2])
matrix_2 = which(as.matrix(network_2)>0,arr.ind = TRUE)
matrix_2 = data.frame("Timestamp1" = 2,"Timestamp2" = 3,"Source" = matrix_2[,1],"Target"= matrix_2[,2])

matrix_3 = rbind(matrix_1, matrix_2)
write.csv(matrix_3,file = "matrix_3.csv",row.names = F)

# Convert the adjacency matrix into a igraph object for each year 
tmp_graph_1 = which(x =as.matrix(network_1)>0,arr.ind = T) # What trades were observed?
tmp_vector = as.numeric(tmp_graph_1) # Construct a vector for the edge list to generate the igraph object
tmp_vector[seq(from = 1, to = length(tmp_vector), by = 2)] = tmp_graph_1[,1] # The edge list is a vector in the form (i ->j, ...)
tmp_vector[seq(from = 2, to = length(tmp_vector), by = 2)] = tmp_graph_1[,2]
# Generate the graph with the edge vector
tmp_graph_1 = make_graph(tmp_vector, directed = TRUE)

# Same procedure for the second year
tmp_graph_2 = which(x = as.matrix(network_2)>0,arr.ind = T) # What trades were observed?
tmp_vector = as.numeric(tmp_graph_2) # Construct a vector for the edge list to generate the igraph object
tmp_vector[seq(from = 1, to = length(tmp_vector), by = 2)] = tmp_graph_2[,1] # The edge list is a vector in the form (i ->j, ...)
tmp_vector[seq(from = 2, to = length(tmp_vector), by = 2)] = tmp_graph_2[,2]
# Generate the graph with the edge vector
tmp_graph_2 = make_graph(tmp_vector, directed = TRUE)

transitivity(tmp_graph_1, type = c("global")) # ratio of the triangles and the connected triples in the graph
transitivity(tmp_graph_2, type = c("global"))

reciprocity(tmp_graph_1) # proportion of mutual connections, in a directed graph 
# -> probability that the opposite counterpart of a directed edge is also included in the graph
reciprocity(tmp_graph_2) 

in_graph_1 = degree_distribution(tmp_graph_1, mode = c("in"))
in_graph_2 = degree_distribution(tmp_graph_2, mode = c("in"))

in_graph_2[36:40] = 0
plot_data = data.frame("degree" =0:(length(in_graph_1)-1), 
                       "prop1" = in_graph_1, "prop2" = in_graph_2)

plot_data = melt(plot_data, id = "degree")
names(plot_data)[2] = "Year"
pdf("in_deg_2.pdf")
ggplot(data = plot_data,aes(x = degree, y = value, fill = Year))+
  geom_bar(stat="identity", position=position_dodge()) +
  theme_pubr() +
  scale_fill_manual(values=c('#000000','#8e8b8b'),labels = c("2016","2017")) +
  ylab("Proportion") +
  xlab("Degree") +
  guides(fill = F) +
  ylim(c(0,0.22)) +
  ggtitle("In-Degree")   +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30))  
  
dev.off()

out_graph_1 = degree_distribution(tmp_graph_1, mode = c("out"))
out_graph_2 = degree_distribution(tmp_graph_2, mode = c("out"))
out_graph_2[32:43] = 0


plot_data = data.frame("degree" =0:(length(out_graph_1)-1), 
                       "prop1" = out_graph_1, "prop2" = out_graph_2)

plot_data = melt(plot_data, id = "degree")
names(plot_data)[2] = "Year"
pdf("out_deg_2.pdf")

ggplot(data = plot_data,aes(x = degree, y = value, fill = Year))+
  geom_bar(stat="identity", position=position_dodge()) +
  theme_pubr() +
  scale_fill_manual(values=c('#000000','#8e8b8b'),labels = c("2016","2017")) +
  ylab("Proportion") +
  ylim(c(0,0.22)) +
  xlab("Degree") +
  guides(fill = F) +
  ggtitle("Out-Degree")   +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30))  

dev.off()


# Density 
sum(as.matrix(network_1) >0)/(89*88)
sum(as.matrix(network_2) >0)/(89*88)

# Repetition
sum(((as.matrix(network_1)) + (as.matrix(network_2))) == 2) /
  sum(as.matrix(network_2))

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ------                          Descriptives                          ------ #
# ------                    Fritz, Lebacher, Kauermann                  ------ #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

library(igraph)
library(ggpubr)
library(reshape2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
rm(list=ls())
load("load_data.RData")
beginning = 2016
ending = 2017

# Which countries are present in all 2 years? 
# Save the subgraph of those networks 
tmp_networks = list()
countries_ind =which(rowSums(EX[,beginning:ending -1949]) == 2)
amk_tmp =  amk[[beginning - 1949]][countries_ind, countries_ind]  

# Delete Isolates 

amk_tmp = amk_tmp[!((rowSums(amk_tmp) == 0) & (colSums(amk_tmp) == 0)), !((rowSums(amk_tmp) == 0) & (colSums(amk_tmp) == 0))]
dim(amk_tmp )
# Save those networks in the list tmp_networks
tmp_networks[[1]] = amk_tmp
matrix_1 = which(tmp_networks[[1]]>0,arr.ind = TRUE)
matrix_1 = data.frame("timestamp" = 1,"Source" = matrix_1[,1],"Target"= matrix_1[,2])

amk_tmp =  amk[[ending - 1949]][countries_ind, countries_ind]  
# Delete Isolates 

amk_tmp = amk_tmp[!((rowSums(amk_tmp) == 0) & (colSums(amk_tmp) == 0)), !((rowSums(amk_tmp) == 0) & (colSums(amk_tmp) == 0))]


# Save those networks in the list tmp_networks
tmp_networks[[2]] = amk_tmp
matrix_2 = which(tmp_networks[[2]]>0,arr.ind = TRUE)
matrix_2 = data.frame("timestamp" = 1,"Source" = matrix_2[,1],"Target"= matrix_2[,2])

tmp_networks[[1]] = amk[[beginning - 1949]][countries_ind, countries_ind]
tmp_networks[[2]] = amk[[ending - 1949]][countries_ind, countries_ind]

# Convert the adjacency matrix into a igraph object for each year 
tmp_graph_1 = which(x = tmp_networks[[1]]>0,arr.ind = T) # What trades were observed?
tmp_vector = as.numeric(tmp_graph_1) # Construct a vector for the edge list to generate the igraph object
tmp_vector[seq(from = 1, to = length(tmp_vector), by = 2)] = tmp_graph_1[,1] # The edge list is a vector in the form (i ->j, ...)
tmp_vector[seq(from = 2, to = length(tmp_vector), by = 2)] = tmp_graph_1[,2]
# Generate the graph with the edge vector
tmp_graph_1 = make_graph(tmp_vector, directed = TRUE)

# Same procedure for the second year
tmp_graph_2 = which(x = tmp_networks[[2]]>0,arr.ind = T) # What trades were observed?
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
in_graph_2[15:17] = 0
plot_data = data.frame("degree" =0:(length(in_graph_1)-1), 
                       "prop1" = in_graph_1, "prop2" = in_graph_2)

plot_data = melt(plot_data, id = "degree")
names(plot_data)[2] = "Year"
pdf("in_deg.pdf")
ggplot(data = plot_data,aes(x = degree, y = value, fill = Year))+
  geom_bar(stat="identity", position=position_dodge()) +
  theme_pubr() +
  scale_fill_manual(values=c('#000000','#8e8b8b'),labels = c("2016","2017")) +
  ylab("Proportion") +
  xlab("Degree") +
  ylim(c(0,0.75)) +
  guides(fill = F) +
  ggtitle("In-Degree") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(family = "Sans",size = 25)) +
  theme(axis.title.y = element_text(family = "Sans",size = 25)) +
  theme(axis.text.x = element_text(family = "Sans",size = 15))  +
  theme(axis.text.y = element_text(family = "Sans",size = 15)) +
  theme(plot.title = element_text(family = "Sans",size = 25)) 
  
dev.off()

out_graph_1 = degree_distribution(tmp_graph_1, mode = c("out"))
out_graph_1[61:64] = 0
out_graph_2 = degree_distribution(tmp_graph_2, mode = c("out"))
plot_data = data.frame("degree" =0:(length(out_graph_1)-1), 
                       "prop1" = out_graph_1, "prop2" = out_graph_2)

plot_data = melt(plot_data, id = "degree")
names(plot_data)[2] = "Year"
pdf("out_deg.png")

ggplot(data = plot_data,aes(x = degree, y = value, fill = Year))+
  geom_bar(stat="identity", position=position_dodge()) +
  theme_pubr() +
  ylim(c(0,0.75)) +
  scale_fill_manual(values=c('#000000','#8e8b8b'),labels = c("2016","2017")) +
  ylab("Proportion") +
  xlab("Degree") +
  guides(fill = F) +
  ggtitle("Out-Degree") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_text(family = "Sans",size = 25)) +
  theme(axis.title.y = element_text(family = "Sans",size = 25)) +
  theme(axis.text.x = element_text(family = "Sans",size = 15))  +
  theme(axis.text.y = element_text(family = "Sans",size = 15)) +
  theme(plot.title = element_text(family = "Sans",size = 25)) 

dev.off()




# Density 
sum(tmp_networks[[1]] >0)/(148*147)
sum(tmp_networks[[2]] >0)/(148*147)

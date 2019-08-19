# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ------                         Survey Paper                           ------ #
# ------                      (Fitting a (S)TERGM)                      ------ #
# ------ Data Set 2: European Research Institution Email Correspondence ------ #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# Preliminaries ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
rm(list=ls())

library(statnet)
library(texreg)
source("../help_functions.R")

# set seed
set.seed(123)

# Read the data -----
data_tmp = read.table("email-Eu-core-temporal-Dept3.txt")
names(data_tmp) = c("from", "to", "time")

data_tmp_trunc = data_tmp[data_tmp$time< 31540000*2,]
# Exclude Self-loops
data_tmp_trunc = data_tmp_trunc[!(data_tmp_trunc$from == data_tmp_trunc$to),]
# Delete all group emails 
allowed_times = summary(factor(data_tmp_trunc$time),maxsum =8831) == 1
allowed_times = as.numeric(levels(factor(data_tmp_trunc$time))[allowed_times])
data_tmp_trunc = data_tmp_trunc[data_tmp_trunc$time %in% allowed_times, ]

# Save the actors of the relational data in a vector 
actors = unique(c(data_tmp_trunc$from, data_tmp_trunc$to))

# Change the actor labels to be 1, ..., n_actors
data_tmp_trunc$from_num = match(data_tmp_trunc$from, actors)
data_tmp_trunc$to_num = match(data_tmp_trunc$to, actors)

# The first aggregated network is the first half of observations  (Day 1 - 365)
data_tmp_trunc1 = data_tmp_trunc[data_tmp_trunc$time <31540000,]
network_1 = matrix(data = 0,nrow = length(actors), ncol = length(actors))
network_1[cbind(data_tmp_trunc1$from_num, data_tmp_trunc1$to_num)] = 1
network_1 = network(network_1)

# The second aggregated network is the second half of observations (Day 366 - 730)
data_tmp_trunc2 =data_tmp_trunc[data_tmp_trunc$time >= 31540000,]
network_2 = matrix(data = 0,nrow = length(actors), ncol = length(actors))
network_2[cbind(data_tmp_trunc2$from_num, data_tmp_trunc2$to_num)] = 1
network_2 = network(network_2)

# Estimation of the TERGM ----
lrecip<-t(as.matrix(network_1))

tergm_email = ergm(network_2~  edgecov(network_1) + 
                     edges+
                     edgecov(lrecip)+
                     gwidegree(log(2), fixed = TRUE)+
                     gwodegree(log(2), fixed = TRUE)+
                     gwesp(log(2), fixed = TRUE),control = control.ergm(MCMC.samplesize = 5000))
load("results_tergm_email.Rdata")
save(tergm_email, file = "results_tergm_email.Rdata")

gof_tergm_email = gof(tergm_email)
pdf("gof_email_tergm.pdf")
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
my_gof(gof_tergm_email, plotlogodds = T, new_names = c("Repetition","Edges","Reciprocity","GWID","GWOD",
                                                 "GWESP"))
dev.off()



# Estimation of the STERGM ----

# Combine them in a list
nets<-network.list(list(network_1,network_2))

# run the STERGM ----
stergm_email<-stergm(nets,formation = ~edges+
                 mutual+
                 gwidegree(log(2),fixed=T)+
                 gwodegree(log(2),fixed=T)+
                 gwesp(log(2),fixed=T),
               dissolution  = ~edges+
                 mutual+
                 gwidegree(log(2),fixed=T)+
                 gwodegree(log(2),fixed=T)+
                 gwesp(log(2),fixed=T),
               times=1:2,estimate = "CMLE",control = control.stergm(CMLE.MCMC.interval = 5000))
gof_stergm_email<-gof(stergm_email)
summary(gof_stergm_email)

pdf("gof_email_tergm.pdf")
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
my_gof(gof_tergm_email, plotlogodds = T, 
       new_names = c("Repetition","Edges","Reciprocity","GWID","GWOD","GWESP"))
                                                                  
dev.off()
pdf("gof_email_stergm.pdf")
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
my_gof(gof_stergm_email$formation, plotlogodds = T, new_names = c("Edges","Reciprocity","GWID","GWOD",
                                                       "GWESP"))
my_gof(gof_stergm_email$dissolution, plotlogodds = T, new_names = c("Edges","Reciprocity","GWID","GWOD",
                                                                  "GWESP"))
dev.off()



# Save the results ----
save.image("stergm_tergm_2.RData")
#load("stergm_tergm_2.RData")

# Extract the coefficients
theta_stergm<-coef(stergm_email)
theta_tergm<-coef(tergm_email)
Mat<-cbind(theta_tergm,c(0,theta_stergm$formation),c(0,theta_stergm$dissolution))

# Save the P-values ----
p_val<-list()
p_val[[2]]<-c(0,summary(stergm_email)$formation$coefs[,4])
p_val[[3]]<-c(0,summary(stergm_email)$dissolution$coefs[,4])
p_val[[1]]<-summary(tergm_email)$coefs[,4]

# Save GOF-Stats----
lgof<-list()
lgof[[2]] <- c(summary(stergm_email)$formation$aic,summary(stergm_email)$formation$mle.lik)
lgof[[3]] <- c(summary(stergm_email)$dissolution$aic,summary(stergm_email)$dissolution$mle.lik)
lgof[[1]] <- c(AIC(tergm_email),logLik(tergm_email))
gof.names <- c("AIC", "Log Likelihood","Sample Size") #names of GOFs

est <- Mat
se <- cbind(summary(tergm_email)$coefs[,2],c(0,summary(stergm_email)$formation$coefs[,2]),c(0,summary(stergm_email)$dissolution$coefs[,2]))
colnames(est) <- c("TERGM","STERGM, Formation","STERGM, Dissolution")

# add row labels:
rownames(est) <- c("lagged Network",rownames(Mat)[-1])

# create a texreg object:
tr <- list()
for (j in 1:ncol(est)) {
  tr[[j]] <- createTexreg(
    coef.names = rownames(est), 
    coef = est[, j], 
    se = se[, j],
    pvalues = p_val[[j]],
    ci.low = numeric(0),
    ci.up = numeric(0),
    gof.names <- c("AIC",  "Log Likelihood"),
    gof =lgof[[j]]
    
  )
}

# for text output:
screenreg(tr, custom.model.names = colnames(est), 
          custom.note = "",stars = c(0.001,  0.01, 0.05),gof=gof,digits = 3)

# for LaTeX output:
texreg(tr, custom.model.names = colnames(est),
       custom.note = "",digits=3,stars = c(0.001,  0.01, 0.05))

# ROC- and PR-Curves 

# First get the change statistics 
change_info_tergm = ergmMPLE(network_2~  edgecov(network_1) + 
                               edges+
                               edgecov(lrecip)+
                               gwidegree(log(2), fixed = TRUE)+
                               gwodegree(log(2), fixed = TRUE)+
                               gwesp(log(2)),output = "array")
# Then construct the linear predictor for calculating the probability of 
# ties conditional on the rest of the network 
eta = lapply(1:length(tergm_email$coef), function(x){tergm_email$coef[x] * change_info_tergm$predictor[,,x]})
eta = Reduce("+",eta)
eta = as.numeric(eta)
eta = eta[-is.na(eta)]
P_tergm = exp(eta)/(1 + exp(eta))


# STERGM
# Save a matrix of the network at t where the diagonals are NA not 0 
tmp_mat = (as.matrix(nets[[1]]))
diag(tmp_mat) = NA
# Use this matrix to indicate which ties are governed 
# by the formation or dissolution model
formation_network = which(tmp_mat == 0, arr.ind = T)
dissolution_network = which(tmp_mat == 1, arr.ind = T)

# The covaraites slightly differ -> bew change statistics 
change_info_stergm = ergmMPLE(network_2~ edges+
                                mutual+
                                gwidegree(log(2),fixed=T)+
                                gwodegree(log(2),fixed=T)+
                                gwesp(log(2),fixed=T),output = "array")


# This part is unchanged
eta_formation = lapply(1:length(stergm_email$formation.fit$coef), 
                       function(x){stergm_email$formation.fit$coef[x] * 
                           change_info_stergm$predictor[cbind(formation_network,x)]})
eta_formation = Reduce("+",eta_formation)
eta_dissolution = lapply(1:length(stergm_email$dissolution.fit$coef), 
                         function(x){stergm_email$dissolution.fit$coef[x] * 
                             change_info_stergm$predictor[cbind(dissolution_network,x)]})
eta_dissolution = Reduce("+",eta_dissolution)
eta_stergm = c(eta_formation, eta_dissolution)
P_stergm = exp(eta_stergm)/(1 + exp(eta_stergm))
real_stergm = c(stergm_email$formation.fit$network[formation_network],
                stergm_email$dissolution.fit$network[dissolution_network])
# Calculate the ROC values
roc_tergm = roc_my(P_tergm,nets[[2]])
roc_stergm = roc_my(P_stergm,real_stergm)

df_roc_tergm = data.frame("specificity" =1- roc_tergm$curve[,1],
                          "sensitivity" = roc_tergm$curve[,2], "Model" = "TERGM")
df_roc_stergm = data.frame("specificity" =1- roc_stergm$curve[,1],
                           "sensitivity" = roc_stergm$curve[,2], "Model" = "STERGM")
df_roc= rbind(df_roc_tergm, df_roc_stergm)

pdf("roc_comparison_2.pdf",width = 8,height =8)

ggplot(data = df_roc ,mapping = aes(x = specificity,y = sensitivity, lty = Model)) +
  geom_step(alpha = 0.8, size=1) +
  scale_x_reverse() +
  geom_abline(slope = 1, intercept = 1, lty = 2, color = "grey40") +
  ylab("Sensitivity")+
  xlab("Specificity") +
  theme_pubr() +
  scale_linetype_manual(name = "", 
                        labels = c(paste0("TERGM (",round(roc_tergm$auc,digits = 3) ,")"),
                                   paste0("STERGM (",round(roc_stergm$auc,digits = 3) ,")")), 
                        values =  c(1,4,2),guide = guide_legend(keywidth = 4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30))

dev.off()

# Calculate the pr values
pr_tergm = pr_my(P_tergm,nets[[2]])
pr_stergm = pr_my(P_stergm,real_stergm)

df_pr_tergm = data.frame("Precision" =pr_tergm$curve[,1],
                         "Recall" = pr_tergm$curve[,2], "Model" = "TERGM")
df_pr_stergm = data.frame("Precision" = pr_stergm$curve[,1],
                          "Recall" = pr_stergm$curve[,2], "Model" = "STERGM")

df_pr= rbind(df_pr_tergm, df_pr_stergm)

pdf("pr_comparison_2.pdf",width = 8,height =8)
ggplot(data = df_pr ,mapping = aes(x = Precision,y = Recall, lty = Model)) +
  geom_step(alpha = 0.8, size=1) +
  ylab("Precision")+
  xlab("Recall") +
  theme_pubr() +
  scale_linetype_manual(name = "", 
                        labels = c(paste0("TERGM (",round(pr_tergm$auc.integral,digits = 3) ,")"),
                                   paste0("STERGM (",round(pr_stergm$auc.integral,digits = 3) ,")")), 
                        values = c(1,4,2),guide = guide_legend(keywidth = 3)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30)) 

dev.off()

